# preparing prediction grids with covariates attached

# attach packages
suppressPackageStartupMessages(require(rgdal, quietly=TRUE, warn.conflicts=FALSE)) # to import and export spatial data
suppressPackageStartupMessages(require(raster, quietly =TRUE, warn.conflicts=FALSE)) # for handling raster maps

# folder contraining all the covariates tif files
gtiffolder <- "../../../GEOdata/ET_1k_Gtif"
# covariates interested   
grid.list <- c("BLUE.tif", "BSAn.tif", "BSAs.tif", "BSAv.tif", "CTI.tif", "ELEV.tif", "EVI.tif", "FPAR.tif", "LAI.tif", "LSTd.tif", "LSTn.tif", "MAP.tif", "MAT.tif", "MIR.tif", "NDVI.tif", "NIR.tif", "RED.tif", "RELIEF.tif", "WSAn.tif", "WSAs.tif", "WSAv.tif")
 
# read in the prediction grids, and attach covariates in it
# only keep the locations which is not cliped

predict_grid_1k_tif <- readGDAL(paste(gtiffolder, "/", "pred_grid_1K.tif", sep=""), silent=TRUE)

predict_grid_1k_coords <- coordinates(predict_grid_1k_tif)[c(predict_grid_1k_tif@data$band1)==1&!is.na(predict_grid_1k_tif@data$band1), ]

predict_grid_1k <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
    mask = rep(1, dim(predict_grid_1k_coords)[1]))
)
  
# use 'raster' to read a geotif from disk
for(i in 1:length(grid.list)){
	cat(paste("extracting", grid.list[i], "\n"))
	rmap_bndry_new <- raster(paste(gtiffolder, "/", grid.list[i], sep="")) # raster is producing small file size in the memory as compared to readGDAL
	predict_grid_1k@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract (
  	x = rmap_bndry_new, # covariates data 
  	y = predict_grid_1k, # original data
  	method = "simple"
	)
}

predict_grid_1k_values <- predict_grid_1k@data[, -1]
predict_grid_1k_values.narm <- predict_grid_1k_values[!is.na(rowMeans(predict_grid_1k_values)), ]
predict_grid_1k_values.narm <- as.matrix(predict_grid_1k_values.narm)
predict_grid_1k_coords <- predict_grid_1k_coords[!is.na(rowMeans(predict_grid_1k_values)), ]

write.table(predict_grid_1k_values.narm, "predcov.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)



