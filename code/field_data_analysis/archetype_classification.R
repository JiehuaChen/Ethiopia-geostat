################
# post-stratification estimator of CMA based on archetype classification
#################

# load in libraries

library(rgdal)
library(proj4)
library(shapefiles)
library(maptools)
library(arm)
library(archetypes)
library(compositions)
library(raster)
require(downloader)

#############
# data downloading
###############

### Load field, grid covariate data and woreda shape files
download("https://www.dropbox.com/s/cioqau3tndqhupg/ET_field_data.csv.zip", "ET_field_data.csv.zip", mode="wb")
unzip("ET_field_data.csv.zip", overwrite=T)
download("https://www.dropbox.com/s/kr4ow5y66kj2itv/ET_1K_Gtif.zip", "ET_1K_Gtif.zip", mode="wb")
unzip("ET_1K_Gtif.zip", overwrite=T)
download("https://www.dropbox.com/s/rysda8349lbpho0/ETH_adm.zip", "ETH_adm.zip", mode = "wb")
unzip("ETH_adm.zip", overwrite=T)


#################
# read in data
##################

# read in field data
fdat <- read.table("ET_field_data.csv", header=T, sep=",")

### Define grid "coordinate reference system" (CRS)
coordinates(fdat) = ~Lon+Lat
proj4string(fdat) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
fdat.laea <- spTransform(fdat, CRS=CRS("+proj=laea +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs"))
coordnames(fdat.laea)<- c("x","y")

### Grid overlay loop
grid.list <- c("BSAv.tif","BSAn.tif","RED.tif","BLUE.tif","NIR.tif","MIR.tif","LSTd.tif","LSTn.tif","EVI.tif","ELEV.tif","RELIEF.tif","MAP.tif","MAT.tif")
for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grid.cov <- raster(paste("ET_1K_Gtif/"+grid.list[i], sep="")) 
  fdat.laea@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grid.cov, 
    y = fdat.laea,
    method = "simple")
}
fdat <- as.data.frame(fdat.laea)

# read in the grids for the whole map
# only keep the locations which is not cliped
predict_grid_1k_tif <- readGDAL("pred_grid_1K.tif")
predict_grid_1k_coords <- coordinates(predict_grid_1k_tif)[c(predict_grid_1k_tif@data$band1)==1&!is.na(predict_grid_1k_tif@data$band1), ]

predict_grid_1k <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
    mask = rep(1, dim(predict_grid_1k_coords)[1]))
)
 
# use 'raster' to read a geotif from disk
for(i in 1:length(grid.list)){
	print(paste("extracting", grid.list[i]))
	rmap_bndry_new <- raster(grid.list[i]) # raster is producing small file size in the memory as compared to readGDAL
	predict_grid_1k@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract (
  	x = rmap_bndry_new, # covariates data 
  	y = predict_grid_1k, # original data
  	method = "simple"
	)
}

# the woreda shape file is in "/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/GEOdata.git/ETH_adm
woreda_shp <-  readShapePoly("ETH_adm/ETH_adm3.shp")

######################
# archetype classficiation
########################

gtifdata <- fdat[, do.call("rbind", strsplit(grid.list, split=".tif"))]
cdata <- as.data.frame(clr(gtifdata))

### Identification of compositional archetypes / endmembers
set.seed(20513)
cdata.arc <- stepArchetypes(data=cdata, k=1:10, nrep=5, verbose=T)
screeplot(cdata.arc)

# Select no. of archetypes, screeplot suggests k=6
cdata.arc6 <- bestModel(cdata.arc[[6]])
predicted_prop_archetype <- predict(cdata.arc6, cdata)
predict_archetype <- as.character(apply(predicted_prop_archetype, 1, which.max))

#######################
# multilevel modeling
#######################

# run multilevel on the predicted_archetypes
M1 <- glmer(CMA~(1|predict_archetype), family=binomial(link="logit"), data=fdat)


###############################
# post-stratification estimator
###############################

# post-stratification estimator, wrapped as a function to be used later
# input: a glmer or lmer object, and ratios of each stratum

post_strat_estimator <- function(M1,table_archetypes){
	CMA_post_mean <- sum(coef(M1)[[1]]*table_archetypes)
	# no error term variance here, since we assume the logistic regreesion
	CMA_post_se <- sqrt(sum((t((se.coef(M1)$fixef)^2 + (se.coef(M1)[[2]])^2))*t(table_archetypes^2)))
	return(list(post_mean = CMA_post_mean, post_se = CMA_post_se))
}

predict.cdata <- as.data.frame(clr(predict_grid_1k@data[, do.call("rbind", strsplit(grid.list, split=".tif"))]))
predict.archetypes <- predict(cdata.arc6, predict.cdata)
predict.archetypes.max <- apply(predict.archetypes, 1, which.max)
archetype.predict.map <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
 predict.archetype =predict.archetypes.max)
)

# the overall post-stratification estimator
table_archetypes <- table(predict.archetypes.max)/length(predict.archetypes.max)
post_strat_estimator(M1, table_archetypes)

# for each woreda
woreda_name <- "Beyeda"
# the index of the woreda in the shapefile
woreda_index <- (1:length(woreda_shp))[woreda_shp@data$NAME_3==woreda_name]
woreda_cliped_shp <-  woreda_shp

woreda_cliped_shp@polygons <-  woreda_cliped_shp@polygons[woreda_index]

# reproject coordinates
k <- length(woreda_cliped_shp@polygons[[1]]@Polygons)
for(i in 1:k){
	projected.coords <- project(woreda_cliped_shp@polygons[[1]]@Polygons[[i]]@coords, "+proj=laea +datum=WGS84 +lat_0=5 +lon_0=20")
  slot(slot(woreda_cliped_shp@polygons[[1]], "Polygons")[[i]], "coords") <- projected.coords
}

# overlay the archetype map with the woreda polygon
predict_archetype_woreda <- archetype.predict.map@data$predict.archetype[!is.na(overlay(archetype.predict.map, woreda_cliped_shp))]
table_archetypes_woreda <- mapply(mean, lapply(1:6, "==",  predict_archetype_woreda))

# the post-stratification estimate for woreda "Beyeda"
post_strat_estimator(M1, table_archetypes_woreda)

