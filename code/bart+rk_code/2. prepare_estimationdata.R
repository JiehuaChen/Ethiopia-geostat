labdata_folder <- "../../../LABdata"

# read in the lab data
lab <- read.table(paste(labdata_folder,"/Samples.csv", sep=""), header=T, sep=",") # lab data
field <- read.table(paste(labdata_folder, "/Profiles.csv", sep=""), header=T, sep=",")# field data
lab_field <- merge(lab, field, by="PID")

# project Lat/Lon profile coordinates in to the LAEA CRS of "etgrid"
coordinates(lab_field) = ~ Lon+Lat # assign coordinates; #proj4string(top.soil) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(lab_field) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") ## +init=epsg:32637" #for transformation, we use the funcrion
# reproject it into lambert
lab_field.laea <- spTransform(lab_field, CRS=CRS("+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
coordnames(lab_field.laea)<- c("x", "y")

#remove duplicates
ids<-zerodist(lab_field.laea)
if(dim(ids)[1]>0){
	lab_field.laea<-lab_field.laea[-c(ids[,1], ids[, 2]),]
}
# save the data
#save(lab_field.laea, file="lab_field.laea.RData")

# extract covariates for lab data locations
grid.list <- c("BLUE.tif", "BSAn.tif", "BSAs.tif", "BSAv.tif", "CTI.tif", "ELEV.tif", "EVI.tif", "FPAR.tif", "LAI.tif", "LSTd.tif", "LSTn.tif", "MAP.tif", "MAT.tif", "MIR.tif", "NDVI.tif", "NIR.tif", "RED.tif", "RELIEF.tif", "WSAn.tif", "WSAs.tif", "WSAv.tif")
 

for (i in 1:length(grid.list)) {
	cat(paste("extracting", grid.list[i], "\n"))
	rmap_bndry_new <- raster(paste(gtiffolder, "/", grid.list[i], sep=""))
	lab_field.laea@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract (
  	x = rmap_bndry_new, # grid data
  	y = lab_field.laea, # original data
  	method = "simple"
	)
}


