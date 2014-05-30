#!/usr/bin/Rscript --vanilla

# clean R workspace
rm(list=ls())

cmd.args <- commandArgs(trailingOnly=TRUE)

soil_property <-  cmd.args[1]
suppressPackageStartupMessages(require(BayesTree, quietly=T, warn.conflicts=FALSE))
suppressPackageStartupMessages(require(gstat, quietly=T, warn.conflicts=FALSE))

# run 1. prepare_grids: prepare grid data for the whole prediction grid
cat("Step 1: Prepare prediction grids with covariates attached \n")
source("1. prepare_grids_withshp.R")

# run 2. prepare_estimationdata
cat("Step 2: Prepare estimation data with covariates attached \n")
source("2. prepare_estimationdata.R")

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
library(BayesTree)

bart.est <- bart(X_lm, Y_lm, x.test=predict_grid_1k_values.narm, sigdf=sigdf.est, sigquant=sigquant.est, k=k.est, ntree=ntree.est, ndpost=500, nskip=1000, keepevery=10)

bart.predict.mean <-  apply(bart.est$yhat.test,2,mean)
bart.predict.sd <- apply(bart.est$yhat.test, 2, sd)
bart.predict <- bart.est$yhat.test

# add the residuals of BART into the dataframe
trained_est <- bart.est$yhat.train
predict_map <- matrix(NA, length(bart.predict.mean), dim(trained_est)[1])

cat("Step 4.2: Ordinary Kriging on the BART residuals \n")

for(k in 1:dim(trained_est)[1]){
    cat(k, "Simulation \n")
	
	lab_field.laea@data["bart.resi"] <- Y_lm - trained_est[k,]
	## STEP TWO: variogram for residuals ____________________________________
	#compute exponential variogram of model residuals
	vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=10E3)
	vgmf.P.resi.exp <- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 10E3, nugget = 0.01))
	
	#Now interpolate model residuals by simple kriging to the nodes of the grid.
	#Simple kriging means that we know the model mean (mu), which is in this case 0.
	predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
	coordinates(predict_grid_1k_coords) =~x+y
	predict_grid_1k_coords@proj4string@projargs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
	krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.exp, nmax = 100, beta=0, na.action=na.omit, nsim=1)

	# compute the final prediction by adding the regression model prediction and the interpoalted residual
	predlogP <- bart.predict[k,] +krige.P.resi@data$sim1
	predict_map[,k] <- exp(predlogP)
}

cat("Step 5: Write the estimated maps in disk, and those maps will be in map_results folder \n")
krige.P.resi$mean = rowMeans(predict_map)
krige.P.resi$sd = apply(predict_map, 1, sd)

map_folder <- "../map_results/"
setwd(map_folder)
#linear+rk_N_upper_ppm_20130626.tif
predlogP.new <- krige.P.resi
gridded(predlogP.new)<-TRUE


writeGDAL (
  dataset=predlogP.new["mean"],
  fname=paste("bart+rk_",soil_property,"_mean_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["sd"],
  fname=paste("bart+rk_",soil_property,"_sd_ppm_20130726.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


