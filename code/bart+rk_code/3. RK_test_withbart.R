# preexisted data: rmap_bndry, lab_field.laea
library(BayesTree)

soil_property <- "P"
map_folder <- "~/ethiosis/map_results/"

# estimate linear model
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <-  lab_field.laea[soil_property]@data[[1]]
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm), ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm), ]
Y_lm <- lab_field.laea[soil_property]@data[[1]]
Y_lm <- ifelse(Y_lm==0, 0.5, Y_lm)
Y_lm <- log(Y_lm)

# prepare prediction covariates with no missing data
source("/Users/jiehuachen/Documents/research/afsis/Spectrum/bart_code/bartfunc.R")
source("/Users/jiehuachen/Documents/research/afsis/Spectrum/bart_code/makeind.R")

bart.est <- bart_saveresults(X_lm, Y_lm, sigquant=0.9, k=2, ntree=50, ndpost=2000, nskip=1000, keepevery=10)

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
vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=20E3)
plot(vg.P.resi)

#fit variogram model for residuals
vgmf.P.resi.sph<- fit.variogram(vg.P.resi, model = vgm(psill = 0.3, model = "Sph", range = 20E3, nugget = 1))
plot(vg.P.resi,vgmf.P.resi.sph, main="sph")
str(vgmf.P.resi.sph) #"SSErr"= num  1.26e-06

vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 20E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.sph, main="exp")
str(vgmf.P.resi.exp) #"SSErr"= num 1.25e-06  --> is the better fit model

vgmf.P.resi.gau<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Gau", range = 20E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.gau, main="gau")
str(vgmf.P.resi.gau) #"SSErr"= num 1.94e-06



#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.
predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
coordinates(predict_grid_1k_coords) =~x+y
predict_grid_1k_coords@proj4string@projargs <- " +proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.sph, nmax = 100, beta=0, na.action=na.omit)


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


#predlogP$RK.var.P<-exp(predlogP$variance.resi)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=x, y=y, fill=RKpred.P))+ coord_equal() # fill mean give the value to the cell with color attached to it

#compute the variance of the regression+kriging prediction error by adding the regression prediction error variance and the kriging variace of the residuals
krige.P.resi$RKlower<-exp(loglower)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=s1, y=s2, fill=RKlower))+ coord_equal()

krige.P.resi$RKupper<-exp(logupper)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=s1, y=s2, fill=RKupper))+ coord_equal()

setwd(map_folder)

predlogP.new <- krige.P.resi
gridded(predlogP.new)<-TRUE
writeGDAL (
  dataset=predlogP.new["RKupper"],
  fname=paste("rkupper.",soil_property,"_new.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKlower"],
  fname=paste("rklower.",soil_property,"_new.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


writeGDAL (
  dataset=predlogP.new["RKpred.P"],
  fname=paste("RKpred.",soil_property,"_new.tif", sep=""),
  drivername= "GTiff",
  type="Float32")



