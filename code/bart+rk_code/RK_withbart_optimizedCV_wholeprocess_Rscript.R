#!/usr/local/lib64/R/bin/Rscript --vanilla


# clean R workspace
rm(list=ls())

cmd.args <- commandArgs(trailingOnly=TRUE)

soil_property <-  cmd.args[1]

suppressPackageStartupMessages(require(BayesTree, quietly=T, warn.conflicts=FALSE))
suppressPackageStartupMessages(require(gstat, quietly=T, warn.conflicts=FALSE))

# run 1. prepare_grids: prepare grid data for the whole prediction grid
cat("Step 1: Prepare prediction grids with covariates attached \n")
source("1. prepare_grids.R")

# run 2. prepare_estimationdata
cat("Step 2: Prepare estimation data with covariates attached \n")
source("2. prepare_estimationdata.R")
ls()
# get cv results
cat("Step 3: Load tuned hyperparameters from Cross-Validation for BART \n")
source("getcvresults.R")


# if the soil_property is not calibrated, we use the default hyperparameter setting of bart funciton

if(sum(results25.sse[,1]==soil_property)==1){
	cv_results <- results25.sse[results25.sse[,1]==soil_property, ]
	ntree.est <- 25
	sigdf.est <- as.numeric(cv_results[2])
	sigquant.est <- as.numeric(cv_results[3])
	k.est <- as.numeric(cv_results[4])
} else{
	ntree.est <- 25
	sigdf.est <- 3 
	sigquant.est <- 0.9
	k.est <- 5
}

# estimate BART model
cat(paste("Step 4.1: Estimating the BART model for the soil property", soil_property, "\n"))
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <-  lab_field.laea[soil_property]@data[[1]]
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
Y_lm <- lab_field.laea[soil_property]@data[[1]]
Y_lm <- log(Y_lm)


# prepare prediction covariates with no missing data
source("bartfunc.R")
source("makeind.R")
dyn.load("src/mbart.so")


setwd("MCMC_results")

bart.est <- bart_saveresults(X_lm, Y_lm, sigdf=sigdf.est, sigquant=sigquant.est, k=k.est, ntree=ntree.est, ndpost=500, nskip=10000, keepevery=10)

cat("Step 4.1: Prediction from the estimated BART model \n")
# bart prediction
xdat.dir <- paste("../../../../GEOdata/ET_1k_Gtif/", "predcov.txt", sep="")
MCMCresults.dir <- paste("", "MCMC*.txt", sep="")
rgy.dir <- paste("", "rgy.txt", sep="")
pipe_run <- (pipe(paste("../src_prediction/predict -x ", xdat.dir, " -f ", "\"", MCMCresults.dir, "\"", " -s \" \" -r ", rgy.dir," -o predictedY.txt", sep="")))
readLines(pipe_run)
close(pipe_run)

bart.predict <- read.table("predictedY.txt", header=TRUE)
bart.predict.mean <- bart.predict[,1]
bart.predict.sd <- bart.predict[,2]

# add the residuals of BART into the dataframe

lab_field.laea@data["bart.resi"] <- Y_lm - bart.est$yhat.train.mean


cat("Step 4.2: Ordinary Kriging on the BART residuals \n")

## STEP TWO: variogram for residuals ____________________________________
#compute exponential variogram of model residuals
vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=10E3)
vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 10E3, nugget = 0.01))


#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.
predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
coordinates(predict_grid_1k_coords) =~x+y
predict_grid_1k_coords@proj4string@projargs <- " +proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.exp, nmax = 100, beta=0, na.action=na.omit)


# copy the kriging results to the SpatialPixelsDataFrame map_prediction
krige.P.resi$pred.resi<-krige.P.resi$var1.pred
krige.P.resi$variance.resi<- krige.P.resi$var1.var

# compute the final prediction by adding the regression model prediction and the interpoalted residual
predlogP <- bart.predict.mean +krige.P.resi$pred.resi

# change to original value
krige.P.resi$RKpred_log <- predlogP
krige.P.resi$RKpred<-exp(predlogP)

# to understand the variance of prediction, we can see lower and upper limits of prediction as alternative of variance map (var1.var)
alpha<-0.05 # this is 95% confidence level.
loglower<-predlogP - qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)
logupper<-predlogP + qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)

# compute the variance of the regression+kriging prediction error by adding the regression prediction error variance and the kriging variace of the residuals
krige.P.resi$RKlower<-exp(loglower)
krige.P.resi$RKupper<-exp(logupper)
krige.P.resi$se_log <- sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)

cat("Step 5: Write the estimated maps in disk, and those maps will be in map_results folder \n")

map_folder <- "../map_results/"
setwd(map_folder)
#linear+rk_N_upper_ppm_20130626.tif
predlogP.new <- krige.P.resi
gridded(predlogP.new)<-TRUE

writeGDAL (
  dataset=predlogP.new["RKupper"],
  fname=paste("bart+rk_",soil_property,"_upper_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKlower"],
  fname=paste("bart+rk_",soil_property,"_lower_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKpred"],
  fname=paste("bart+rk_",soil_property,"_mean_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKpred_log"],
  fname=paste("bart+rk_",soil_property,"_mean_log_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["se_log"],
  fname=paste("bart+rk_",soil_property,"_se_log_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


