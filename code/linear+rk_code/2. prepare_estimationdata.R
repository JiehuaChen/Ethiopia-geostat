labdata_folder <- "~/ethiosis/LABdata"
gtifdata_folder <- "~/ethiosis/GEOdata"

# read in the lab data
setwd(labdata_folder)

lab <- read.table("Samples.csv", header=T, sep=",") # lab data
field <- read.table("Profiles.csv", header=T, sep=",")# field data
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
setwd(gtifdata_folder)


grid.list <- c("BLUE.tif","CTI.tif", "ELEV.tif", "EVI.tif", "LAI.tif", "LCOV.tif", "LSTd.tif", "LSTn.tif","MAP.tif", "MAT.tif", "MFI.tif", "MIR.tif", "NDVI.tif", "NIR.tif", "NPP.tif","RED.tif", "RELIEF.tif")
 
for (i in 1:length(grid.list)) {
	print(paste("extracting", grid.list[i]))
	rmap_bndry_new <- raster(grid.list[i]) # raster is producing small file size in the memory as compared to readGDAL
	lab_field.laea@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract (
  	x = rmap_bndry_new, # evi data 
  	y = lab_field.laea, # original data
  	method = "simple"
	)
}
lab_field.laea@data["x"] <- lab_field.laea@coords[, "x"]
lab_field.laea@data["y"] <- lab_field.laea@coords[, "y"]

