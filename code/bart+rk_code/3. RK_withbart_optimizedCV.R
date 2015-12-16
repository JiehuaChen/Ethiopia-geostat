# preexisted data: 
library(gstat)
library(rbart)


soil_property <- "Ca"

# get cv results
source("getcvresults.R")
# prepare prediction covariates with no missing data



num_draw <- 50

print(paste("estimating", soil_property))
cv_results <- results25.sse[results25.sse[,1]==soil_property, ]
ntree.est <- 25
sigdf.est <- as.numeric(cv_results[2])
sigquant.est <- as.numeric(cv_results[3])
k.est <- as.numeric(cv_results[4])

# estimate BART model
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <-  lab_field.laea[soil_property]@data[[1]]
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
Y_lm <- lab_field.laea[soil_property]@data[[1]]
Y_lm <- log(Y_lm)

mcmcresults_dir <- paste("MCMCresults_", soil_property, sep="")
if(file.exists(mcmcresults_dir)){
    setwd(mcmcresults_dir)
}else{
    dir.create(mcmcresults_dir)
    setwd(mcmcresults_dir)
}

bart.est <- bart_saveresults(X_lm, Y_lm, sigdf=sigdf.est, sigquant=sigquant.est, k=k.est, ntree=ntree.est, ndpost=500, nskip=10000, keepevery=10, outformat="json")
setwd("..")


# make predictions

rgy.dir <- file.path(paste("MCMCresults_", soil_property, sep=""), "rgy.txt")

forestpath_dir <- paste("MCMCresults_", soil_property, sep="")

bart_prediction <- predict_func(rangefile=rgy.dir, tiffilename="", forestpath=forestpath_dir, rframe=predict_grid_1k_values_withGID)


bart_prediction_draws <- log(bart_prediction[,2:51])
bart_prediction_mean <- rowMeans(bart_prediction_draws)
bart_prediction_sd <- apply(bart_prediction_draws, 1, sd) 



# add the residuals of BART into the dataframe

lab_field.laea@data["bart.resi"] <- Y_lm - bart.est$yhat.train.mean
## STEP TWO: variogram for residuals ___________________________________
#compute exprimental variogram of model residuals
vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=10E3)
vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 10E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.exp, main="exp")
str(vgmf.P.resi.exp) #"SSErr"= num 1.25e-06  --> is the better fit model


#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.
predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
coordinates(predict_grid_1k_coords) =~x+y
proj4string(predict_grid_1k_coords) = CRS(proj4string(lab_field.laea)) 
krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.exp, nmax = 100, beta=0, na.action=na.omit)


#Copy the kriging results to the SpatialPixelsDataFrame map_prediction
krige.P.resi$pred.resi<-krige.P.resi$var1.pred
krige.P.resi$variance.resi<- krige.P.resi$var1.var

#________________________________________________


# STEP FOUR:prediction by adding regressional model prediction and interpolated residual __________________________

#Compute the final prediction by adding the regression model prediction and the interpoalted residual

predlogP <- bart_prediction_mean +krige.P.resi$pred.resi
# change to original value
krige.P.resi$RKpred_log <- predlogP
krige.P.resi$RKpred<-exp(predlogP)

# to understand the variance of prediction, we can see lower and upper limits of prediction as alternative of variance map (var1.var)
alpha<-0.05 # this is 95% confidence level.
loglower<-predlogP - qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart_prediction_sd^2+ krige.P.resi$var1.var)
logupper<-predlogP + qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart_prediction_sd^2+ krige.P.resi$var1.var)

#compute the variance of the regression+kriging prediction error by adding the regression prediction error variance and the kriging variace of the residuals
krige.P.resi$RKlower<-exp(loglower)
krige.P.resi$RKupper<-exp(logupper)
krige.P.resi$se_log <- sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)

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


