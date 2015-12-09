### Post-stratification of cultivated area surveys of Ethiopia
# M. Walsh & J. Chen, Sep. 2013 

### Set local working directory e.g.,
dat_dir <- "/Users/markuswalsh/Documents/Ethiopia/Field data/Analysis"
setwd(dat_dir)

### Load packages
require(downloader)
require(arm)
require(rgdal)
require(raster)

### Load field data and prediction grid(s)
download("https://www.dropbox.com/s/cioqau3tndqhupg/ET_field_data.csv.zip", "ET_field_data.csv.zip", mode="wb")
unzip("ET_field_data.csv.zip", overwrite=T)
download("https://www.dropbox.com/s/ly0e7jqwndzw26i/ET_CMA_pred.zip", "ET_CMA_pred_grid.zip", mode="wb")
unzip("ET_CMA_pred.zip", overwrite=T)
fdat <- read.table("ET_field_data.csv", header=T, sep=",")

### Define grid "coordinate reference system" (CRS)
coordinates(fdat) = ~Lon+Lat
proj4string(fdat) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
fdat.laea <- spTransform(fdat, CRS=CRS("+proj=laea +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs"))
coordnames(fdat.laea)<- c("x","y")

### Grid overlay loop
grid.list <- c("CMA_prob.tif", "CMA_se.tif")
for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grid.cov <- raster(grid.list[i]) 
  fdat.laea@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grid.cov, 
    y = fdat.laea,
    method = "simple")
}
fdat <- as.data.frame(fdat.laea)

### Stratum definitions (str)
str <- c(0.1,0.25,0.5,0.75,0.9)
attach(fdat)
fdat$str[CMA_prob <= str[1]] <- 1
fdat$str[CMA_prob > str[1] & CMA_prob <= str[2]] <- 2
fdat$str[CMA_prob > str[2] & CMA_prob <= str[3]] <- 3
fdat$str[CMA_prob > str[3] & CMA_prob <= str[4]] <- 4
fdat$str[CMA_prob > str[4] & CMA_prob <= str[5]] <- 5
fdat$str[CMA_prob > str[5]] <- 6
detach(fdat)

### Post-stratification glmm model
M <- lmer(CMA~1+(1|str), family=binomial(link="logit"), data=fdat)
display(M)
m.coef <- coef(M)
m.se <- se.coef(M)
coefplot(m.coef$str[,1], m.se$str[,1], varnames=rownames(m.coef$str), xlim=c(-6,6), CI=2, cex.var=1, cex.pts=1.2, main="")

### Delineate strata
grids <- readGDAL(grid.list[1])
names(grids)[1] <- sub(".tif", "", grid.list[1])
for (i in grid.list[-1]) {
  grids@data[sub(".tif", "", i[1])] <- readGDAL(paste(i))$band1
}
etgrid <- as.data.frame(grids)

attach(etgrid)
etgrid$str[CMA_prob <= str[1]] <- 1
etgrid$str[CMA_prob > str[1] & CMA_prob <= str[2]] <- 2
etgrid$str[CMA_prob > str[2] & CMA_prob <= str[3]] <- 3
etgrid$str[CMA_prob > str[3] & CMA_prob <= str[4]] <- 4
etgrid$str[CMA_prob > str[4] & CMA_prob <= str[5]] <- 5
etgrid$str[CMA_prob > str[5]] <- 6
detach(etgrid)

# Strata map
str_map <- SpatialPixelsDataFrame(points=etgrid[c("x","y")], data=etgrid, proj4string=CRS("+proj=laea +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs"))
writeGDAL(data=str_map["str"], fname="CMA_str.tif", drivername="GTiff", type="Byte", mvFlag=0)

### Post-stratification estimates
# Function inputs: lmer model object (M), and stratum weights (i.e., "proportional areas" = str_wgt)
post_str_est <- function(M, str_wgt){
  post_mean <- sum(coef(M)[[1]]*str_wgt)
  post_se <- sqrt(sum((t((se.coef(M)$fixef)^2 + (se.coef(M)[[2]])^2))*t(str_wgt^2)))
  post_pred <- exp(post_mean)/(1+exp(post_mean))
  post_lower <- exp(post_mean-2*post_se)/(1+exp(post_mean-2*post_se))
  post_upper <- exp(post_mean+2*post_se)/(1+exp(post_mean+2*post_se))
  return(list(post_pred = post_pred, lower = post_lower, upper = post_upper))
}
str_wgt <- table(etgrid$str)/length(etgrid$str)
post_str_est(M, str_wgt)

