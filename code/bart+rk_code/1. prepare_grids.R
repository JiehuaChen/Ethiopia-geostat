# preparing prediction grids with covariates attached

# attach packages
suppressPackageStartupMessages(require(rgdal, quietly=TRUE, warn.conflicts=FALSE)) # to import and export spatial data
suppressPackageStartupMessages(require(raster, quietly =TRUE, warn.conflicts=FALSE)) # for handling raster maps

GID <- function(location){
    res.pixel <- 1000
    xgid <- floor(location[,1]/res.pixel)
    ygid <- floor(location[,2]/res.pixel)
    gidx <- ifelse(location[,1]<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
    gidy <- ifelse(location[,2]<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
    GID <- paste(gidx, gidy, sep="-")
    return(list(xgid=xgid, ygid=ygid, GID=GID))
}


# folder contraining all the covariates tif files
gtiffolder <- "../../../../../remotesensing_data/ET_grids/"
tiffiles <- list.files(gtiffolder, ".tif", full.names=T)
grids <- stack(tiffiles)
# covariates interested   
 
# read in the prediction grids, and attach covariates in it
# only keep the locations which is not cliped

predict_grid_1k_tif <- readGDAL(paste("../../../GEOdata/ET_1K_Gtif/pred_grid_1K.tif", sep=""), silent=TRUE)

predict_grid_1k_coords <- coordinates(predict_grid_1k_tif)[c(predict_grid_1k_tif@data$band1)==1&!is.na(predict_grid_1k_tif@data$band1), ]

predict_grid_1k <- SpatialPointsDataFrame(
  coords = predict_grid_1k_coords,
  data = data.frame(
  mask = rep(1, dim(predict_grid_1k_coords)[1]))
)
  
predict_grid_1k_values <- extract(grids, predict_grid_1k)
predict_grid_1k_values.narm <- predict_grid_1k_values[!is.na(rowMeans(predict_grid_1k_values)), ]
predict_grid_1k_values.narm <- as.matrix(predict_grid_1k_values.narm)
predict_grid_1k_coords <- predict_grid_1k_coords[!is.na(rowMeans(predict_grid_1k_values)), ]

predict_grid_1k_GID <- GID(predict_grid_1k_coords)
predict_grid_1k_values_withGID  <- data.frame(GID = predict_grid_1k_GID[[3]], predict_grid_1k_values.narm)



