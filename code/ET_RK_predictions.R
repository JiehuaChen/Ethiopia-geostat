### RK prediction workflow for topsoil properties of agricultural areas in Ethiopia
# Soil survey data provided by: Crop Nutrition Laboratory Services, Nairobi, http://www.cropnuts.com
# J. Chen & M. Walsh, June 2013

### Required packages
# set.repositories(Ind=1:2)
# install.packages(c("downloader","rgdal","raster","gstat","MASS","rpart","randomForest"), dependencies=T)
library(downloader)
library(rgdal)
library(raster)
library(gstat)
library(MASS)
library(randomForest)

### Set local working directories, e.g.
dat_dir <- "/Users/markuswalsh/Documents/Ethiopia/CNLS/Analyses"
cov_dir <- "/Users/markuswalsh/Documents/Ethiopia/ET_GeoDat"
map_dir <- "/Users/markuswalsh/Documents/Ethiopia/ET_Soil_Pred"

### Load soil survey data
setwd(dat_dir)
# download("https://www.dropbox.com/s/x7il9g863dwg1h8/CNLS_dat.zip, "CNLS_dat.zip", mode="wb")
# unzip("CNLS_dat.zip", overwrite=T)
samp <- read.table("Samples.csv", header=T, sep=",")
prof <- read.table("Profiles.csv", header=T, sep=",")
cnlsdat <- merge(prof, samp, by="PID")

### Project Lat/Lon survey coordinates to LAEA CRS
coordinates(cnlsdat) = ~Lon+Lat
proj4string(cnlsdat) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
cnlsdat.laea <- spTransform(cnlsdat, CRS=CRS("+proj=laea +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs"))
coordnames(cnlsdat.laea)<- c("x","y")

### Remove geo-duplicates
ids <- zerodist(cnlsdat.laea)
if (dim(ids)[1]>0){
	cnlsdat.laea<-cnlsdat.laea[-ids[,2],]
}

### Overlay gridded covariates with soil survey locations
# Specify covariate grids to include in subsequent analyses
setwd(cov_dir)
grid.list <- c("BLUE.tif","CTI.tif", "ELEV.tif", "EVI.tif", "LAI.tif", "LSTd.tif", "LSTn.tif","MAP.tif", "MAT.tif", "MIR.tif", "NIR.tif","RED.tif", "RELIEF.tif")

# Grid overlay loop
for (i in 1:length(grid.list)){
	print(paste("extracting", grid.list[i]))
	grid.cov <- raster(grid.list[i]) 
	cnlsdat.laea@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
  	x = grid.cov, 
  	y = cnlsdat.laea,
  	method = "simple")
}
cnlsdat.laea@data["x"] <- cnlsdat.laea@coords[,"x"]
cnlsdat.laea@data["y"] <- cnlsdat.laea@coords[,"y"]

### Alternative regression models
# Regression formula set-up: specify soil property of interest as "prop", e.g.,
prop <- "SOC"
yvar <- log(cnlsdat.laea[prop]@data[[1]])

xvar <- paste0(strsplit(grid.list, split=".tif"))
fmla <- as.formula(paste("yvar~", paste(xvar, collapse="+")))

# Stepwise linear regression (using the MASS package)
yvar.lm <- lm(fmla, data=cnlsdat.laea)
yvar.st <- stepAIC(yvar.lm, direction="both")
summary(yvar.st)

# Random forests (using randomForest)
set.seed(230613)
yvar.rf <- randomForest(fmla, data=cnlsdat.laea, mtry=3, importance=T, proximity=T, na.action=na.omit)
print(yvar.rf)
importance(yvar.rf)

# BART model
xtrain <- as.matrix(cnlsdat.laea@data[,xvar])
# delete missing rows
na.index <- (is.na(rowMeans(cbind(xtrain, yvar))))*(1:length(yvar))
yvar <- yvar[-na.index]
xtrain <- xtrain[-na.index,]
bart.est <- bart(xtrain, yvar, ntree=100, sigquant=0.9, k=2, power=2, base=0.95, nskip=1000)

# residuals from bart, which can be used for RK
residuals.bart <- yvar - bart.est$yhat.train.mean


### Set up the Region Of Interest (ROI) for prediction, and merge with covariates.
# Note that this step takes about 10-15 min to complete once the Gtif files have been dowloaded.
setwd(cov_dir)
# download("https://www.dropbox.com/s/1j8r1zdf7wkmvgx/ET_Gtif.zip, "ET_Gtif.zip", mode="wb")
# unzip("ET_Gtif.zip", overwrite=T)

# Read prediction grid
predict_grid_1k <- readGDAL("predgrid_1k.tif")
predict_grid_1k_coords <- coordinates(predict_grid_1k)[as.numeric(predict_grid_1k@data$band1)==1&!is.na(predict_grid_1k@data),]

# Convert predition grid to SpatialPointsDataFrame 
predict_grid_1k <- SpatialPointsDataFrame(
	coords = predict_grid_1k_coords,
  	data = data.frame(
  	mask = rep(1,dim(predict_grid_1k_coords)[1])))
  
# Attach covariates included in "grid.list"
for(i in 1:length(grid.list)){
	print(paste("extracting", grid.list[i]))
	pred_map <- raster(grid.list[i])
	predict_grid_1k@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
  	x = pred_map, 
  	y = predict_grid_1k,
  	method = "simple")
}
