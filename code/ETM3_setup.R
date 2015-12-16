#' Mehlich-3 nutrient and nutrient mass balance data setup
#' December 2015
# Markus Walsh

# install.packages(c("downloader","compositions","raster","rgdal"), dependencies=T)
require(downloader)
require(compositions)
require(rgdal)
require(raster)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("ETM3_data", showWarnings=F)
setwd("./ETM3_data")

# Download grids
download("https://www.dropbox.com/s/d8gkwuwpb28l0ty/ET_grids.zip?dl=0", "ET_grids.zip", mode="wb")
unzip("ET_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Download survey data
unzip("ET_Mehlich3.zip", overwrite=T)
prof <- read.table("Profiles.csv", header=T, sep=",") ## survey locations and Woreda names
samp <- read.table("Samples.csv", header=T, sep=",") ## sample data
etm3 <- merge(prof, samp, by="PID")

# Compositional analysis setup --------------------------------------------
vars <- c("PID","Woreda","Lat","Lon","P","K","S","Ca","Mg")
etnb <- na.omit(etm3[vars])
fpart <- c("P","K","S","Ca","Mg") ## all values in mg/kg
etnb$Fv <- 1000000-rowSums(etnb[fpart]) ## calculates "fill value" (Fv), in mg/kg soil
cpart <- c("P","K","S","Ca","Mg","Fv")

# Sequential binary partion & integrated log ratio (ilr) transform
cdata <- acomp(etnb[cpart])
bpart <- t(matrix(c( 1, 1, 1, 1, 1,-1,
                     1,-1, 1,-1,-1, 0,
                     0, 1, 0,-1,-1, 0,
                     1, 0,-1, 0, 0, 0,
                     0, 0, 0, 1,-1, 0), ncol=6, nrow=5, byrow=T))
CoDaDendrogram(X=acomp(cdata), signary=bpart) ## mass balance mobile graph  			
idata <- as.data.frame(ilr(cdata, V=bpart))
colnames(idata) <- c("V0","V3","V4","V5","V6")
etnb <- cbind(etnb, idata)

# Overlay with gridded covariates -----------------------------------------
# Project survey coords to grid CRS
etnb.proj <- as.data.frame(project(cbind(etnb$Lon, etnb$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(etnb.proj) <- c("x","y")
etnb <- cbind(etnb, etnb.proj)
coordinates(etnb) <- ~x+y
projection(etnb) <- projection(grids)

# Extract gridded variables at survey locations
etnbgrid <- extract(grids, etnb)
etm3 <- as.data.frame(etnb)
etm3 <- cbind.data.frame(etm3, etnbgrid)
etm3 <- na.omit(etm3)

# Write data file ---------------------------------------------------------
write.csv(etm3, "etm3.csv", row.names=F)
