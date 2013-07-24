# preexisted data: rmap_bndry, lab_field.laea
library(BayesTree)

soil_property <- "pH"
cv_results <- results25.sse[results25.sse[,1]==soil_property, ]
ntree.est <- 25
sigdf.est <- as.numeric(cv_results[2])
sigquant.est <- as.numeric(cv_results[3])
k.est <- as.numeric(cv_results[4])

map_folder <- "~/ethiosis/map_results/"
map_folder <- "/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/spatial.git/code/bart+rk_code/map_results/"
# estimate linear model
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <-  lab_field.laea[soil_property]@data[[1]]
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm)&Y_lm>0, ]
Y_lm <- lab_field.laea[soil_property]@data[[1]]
Y_lm <- log(Y_lm)


# prepare prediction covariates with no missing data
source("/Users/jiehuachen/Documents/research/afsis/Spectrum/bart_code/bartfunc.R")
source("/Users/jiehuachen/Documents/research/afsis/Spectrum/bart_code/makeind.R")

dyn.load("/Users/jiehuachen/Documents/research/afsis/Spectrum/bart_code/bart/src/mbart.so")

setwd("/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/spatial.git/code/bart+rk_code/MCMC_results")


bart.est <- bart_saveresults(X_lm, Y_lm, sigdf=sigdf.est, sigquant=sigquant.est, k=k.est, ntree=ntree.est, ndpost=500, nskip=10000, keepevery=10)

# bart prediction
xdat.dir <- paste("../../../../GEOdata.git/ET_1k_Gtif/", "predcov.txt", sep="")
MCMCresults.dir <- paste("", "MCMC*.txt", sep="")
rgy.dir <- paste("", "rgy.txt", sep="")
pipe_run <- (pipe(paste("/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/spatial.git/code/bart+rk_code/src_prediction/predict -x ", xdat.dir, " -f ", "\"", MCMCresults.dir, "\"", " -s \" \" -r ", rgy.dir," -o predictedY.txt", sep="")))
readLines(pipe_run)
close(pipe_run)

bart.predict <- read.table("predictedY.txt", header=TRUE)
bart.predict.mean <- bart.predict[,1]
bart.predict.sd <- bart.predict[,2]

# add the residuals of BART into the dataframe

lab_field.laea@data["bart.resi"] <- Y_lm - bart.est$yhat.train.mean

## STEP TWO: variogram for residuals ____________________________________
#compute exprimental variogram of model residuals
vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=10E3)
vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 10E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.exp, main="exp")
str(vgmf.P.resi.exp) #"SSErr"= num 1.25e-06  --> is the better fit model


#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.
predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
coordinates(predict_grid_1k_coords) =~x+y
predict_grid_1k_coords@proj4string@projargs <- " +proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.exp, nmax = 100, beta=0, na.action=na.omit)


#Copy the kriging results to the SpatialPixelsDataFrame map_prediction
krige.P.resi$pred.resi<-krige.P.resi$var1.pred
krige.P.resi$variance.resi<- krige.P.resi$var1.var

#________________________________________________


# STEP FOUR:prediction by adding regressional model prediction and interpolated residual __________________________

#Compute the final prediction by adding the regression model prediction and the interpoalted residual

predlogP <- bart.predict.mean +krige.P.resi$pred.resi
# change to original value
krige.P.resi$RKpred<-exp(predlogP)

# to understand the variance of prediction, we can see lower and upper limits of prediction as alternative of variance map (var1.var)
alpha<-0.05 # this is 95% confidence level.
loglower<-predlogP - qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)
logupper<-predlogP + qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(bart.predict.sd^2+ krige.P.resi$var1.var)

#compute the variance of the regression+kriging prediction error by adding the regression prediction error variance and the kriging variace of the residuals
krige.P.resi$RKlower<-exp(loglower)
krige.P.resi$RKupper<-exp(logupper)

setwd(map_folder)
#linear+rk_N_upper_ppm_20130626.tif
predlogP.new <- krige.P.resi
gridded(predlogP.new)<-TRUE

writeGDAL (
  dataset=predlogP.new["RKupper"],
  fname=paste("bart+rk_",soil_property,"_upper_ppm_20130614.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKlower"],
  fname=paste("bart+rk_",soil_property,"_lower_ppm_20130614.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKpred"],
  fname=paste("bart+rk_",soil_property,"_mean_ppm_20130614.tif", sep=""),
  drivername= "GTiff",
  type="Float32")



