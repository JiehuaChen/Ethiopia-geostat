########
# 
########
library(proj4)
library(rgdal)
library(raster)

fielddata_folder <- "../../../Fielddata"


fdat <- read.table(paste(fielddata_folder, "/ET_field_data.csv", sep=""), header=T, sep=",")

field <- read.table(paste(fielddata_folder, "/ethiopia_usaid.csv", sep=""), header=T, sep=",")# field data
field <- as.data.frame(field)

# delete points falling out of ethiopia boundary
field <- field[field$X>30, ]
field <- aggregate(CULTIVATION~X+Y, field, mean)
field$CULTIVATION <- ifelse(field$CULTIVATION>0, 1, 0)
fdat <- rbind(field, cbind(X=fdat$Lon, Y=fdat$Lat, CULTIVATION = fdat$CMA))

# project Lat/Lon profile coordinates in to the LAEA CRS of "etgrid"
fdat.laea <- as.data.frame(project(cbind(fdat$X, fdat$Y), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(fdat.laea) <- c("x","y")
fdat <- cbind(fdat, fdat.laea)

### Specify grid cell ID's (GID's) and center point coordinates
# Define pixel resolution (res.pixel, in m)
# res.pixel <- 1000

# Gtif data
# get all the neighoring grid locations
locations_gid <- matrix(NA, dim(fdat.laea)[1], 18)

gtiffolder <- "../../../GEOdata/ET_1k_Gtif_CMA/"
# covariates interested   
grid.list <- list.files(gtiffolder, pattern = "\\.tif$")[list.files(gtiffolder, pattern = "\\.tif$")!= "pred_grid_1K.tif"]
grid.list.loc <- paste(gtiffolder, grid.list, sep="") 

# read in the prediction grids, and attach covariates in it
nm <-  grid.list.loc[1]
#raster fileinfo
tif.info <- GDALinfo(nm, silent=TRUE)
totalrows <- tif.info[1]
totalcols <- tif.info[2]
bandsnum <- tif.info[3]
origin.x <- tif.info[4]
res.x <- tif.info[6]
res.y <- tif.info[7]
if(attr(tif.info, "ysign")==-1){
	#origin of data array should be the left upper corner
	origin.y <- tif.info[5] + totalrows*res.y
}else{
	origin.y <- tif.info[5]	
}
	
pixeln.x <- ceiling((fdat.laea$x-origin.x)/res.x)-1
pixeln.y <- ceiling((origin.y-fdat.laea$y)/res.y)-1
neighbor.x <- c(-1, 0, 1, -1,0,1, -1, 0, 1)
neighbor.y <- c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
locations_gid  <-  matrix(NA, dim(fdat.laea)[1], 18)

for(i in 1:9){
	locations_gid[, (2*i-1)] <- (pixeln.x+neighbor.x[i])*res.x + origin.x + res.x/2
	locations_gid[, (2*i)] <- origin.y - ((pixeln.y+neighbor.y[i])*res.y + res.y/2)
}

colnames(locations_gid) <- c("2nb1x","2nb1y" ,"1nb1x", "1nb1y", "2nb2x", "2nb2y", "1nb2x", "1nb2y","locx", "locy", "1nb3x", "1nb3y", "2nb3x", "2nb3y", "1nb4x", "1nb4y", "2nb4x", "2nb4y")

for(j in 1:9){
	temp <- cbind(fdat[,1], x=locations_gid[, (2*j-1)], y = locations_gid[, (2*j)])
	temp <- as.data.frame(temp)
	coordinates(temp) = ~ x+y
	for (i in 1:length(grid.list)) {
		print(paste("extracting", grid.list.loc[i]))
		rmap_bndry_new <- raster(grid.list.loc[i]) # raster is producing small file size in the memory as compared to readGDAL
		temp@data[paste(strsplit(grid.list[i], split=".tif")[[1]],substring(colnames(locations_gid)[(2*j-1)],1, 4), sep="_")] <- extract (
		  	x = rmap_bndry_new, # evi data 
		  	y = temp, # original data
		  	method = "simple"
		 )
	}
	fdat <- cbind(fdat, temp@data[,-1])
}

# predict CMA
cma_data <- cbind(CMA= fdat$CMA, fdat[, (dim(fdat)[2]+2):(dim(fdat)[2])])
cma_data <- na.omit(cma_data)
# BART prediction
library(BayesTree)
load("LC.RData")

x <- cma_data[,-1]
y <- cma_data[,1]

bart.est <- bart(x, y,x.test = as.data.frame(predict_grid_1k_values.narm),   ndpost=500, nskip=2000, keepevery=10)

save.image("bart.est.cma.RData")

predict.bart.mean <-  apply(pnorm(bart.est$yhat.test),2,mean)
predict.bart.sd <- apply(pnorm(bart.est$yhat.test), 2, sd)

predict_grid_1k <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
 predict.mean =predict.bart.mean,
 predict.se = predict.bart.sd)
)

gridded(predict_grid_1k) <- TRUE

writeGDAL(
  	dataset = predict_grid_1k["predict.mean"],
  	fname ="/data6/EthiopiaData/cmapredict.tif",
  	drivername = "GTiff",
  	type = "Float32",
  	Overwrite<- TRUE)
  	

  	  	  	
writeGDAL(
  	dataset = predict_grid_1k["predict.se"],
  	fname ="/data6/EthiopiaData/cmapredict_se.tif",
  	drivername = "GTiff",
  	type = "Float32",
  	Overwrite<- TRUE)
  	
 


