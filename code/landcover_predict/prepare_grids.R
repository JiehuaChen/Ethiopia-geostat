# preparing prediction grids with covariates attached

# attach packages
library(rgdal) # to import and export spatial data
library(raster) # for handling raster maps


# covariates interested   

gtiffolder <- "../../../GEOdata/ET_1k_Gtif_CMA/"
# covariates interested   
grid.list <- list.files(gtiffolder, pattern = "\\.tif$")[list.files(gtiffolder, pattern = "\\.tif$")!= "pred_grid_1K.tif"]
grid.list.loc <- paste(gtiffolder, grid.list, sep="")



predict_grid_1k_tif <- readGDAL(paste(gtiffolder, "/", "pred_grid_1K.tif", sep=""), silent=TRUE)

predict_grid_1k_coords <- coordinates(predict_grid_1k_tif)[c(predict_grid_1k_tif@data$band1)==1&!is.na(predict_grid_1k_tif@data$band1), ]

predict_grid_1k <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
    mask = rep(1, dim(predict_grid_1k_coords)[1]))
)

nm <- grid.list.loc[1]

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
	
pixeln.x <- ceiling((predict_grid_1k@coords[,1]-origin.x)/res.x)-1
pixeln.y <- ceiling((origin.y-predict_grid_1k@coords[,2])/res.y)-1
neighbor.x <- c(-1, 0, 1, -1,0,1, -1, 0, 1)
neighbor.y <- c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
pred_locations_gid  <-  matrix(NA, dim(predict_grid_1k@coords)[1], 18)

for(i in 1:9){
	pred_locations_gid[, (2*i-1)] <- (pixeln.x+neighbor.x[i])*res.x + origin.x + res.x/2
	pred_locations_gid[, (2*i)] <- origin.y - ((pixeln.y+neighbor.y[i])*res.y + res.y/2)
}
colnames(pred_locations_gid) <- c("2nb1x","2nb1y" ,"1nb1x", "1nb1y", "2nb2x", "2nb2y", "1nb2x", "1nb2y","locx", "locy", "1nb3x", "1nb3y", "2nb3x", "2nb3y", "1nb4x", "1nb4y", "2nb4x", "2nb4y")

for(j in 1:9){
	temp <- cbind(predict_grid_1k$mask, x=pred_locations_gid[, (2*j-1)], y = pred_locations_gid[, (2*j)])
	temp <- as.data.frame(temp)
	coordinates(temp) = ~ x+y 
	for (i in 1:length(grid.list)) {
		print(paste("extracting", grid.list.loc[i]))
		predict_grid_1k_new <- raster(grid.list.loc[i]) # raster is producing small file size in the memory as compared to readGDAL
		temp@data[paste(strsplit(grid.list[i], split=".tif")[[1]],substring(colnames(pred_locations_gid)[(2*j-1)],1, 4), sep="_")] <- extract (
		  	x = predict_grid_1k_new, # evi data 
		  	y = temp, # original data
		  	method = "simple"
		 )
	}
	predict_grid_1k <- cbind(predict_grid_1k, temp@data[,-1])
}

predict_grid_1k_values <- predict_grid_1k[, -(1:3)]
predict_grid_1k_values.narm <- predict_grid_1k_values[!is.na(rowMeans(predict_grid_1k_values)), ]
predict_grid_1k_values.narm <- as.matrix(predict_grid_1k_values.narm)
predict_grid_1k_coords <- predict_grid_1k_coords[!is.na(rowMeans(predict_grid_1k_values)), ]

save.image("GEOdata.RData")

