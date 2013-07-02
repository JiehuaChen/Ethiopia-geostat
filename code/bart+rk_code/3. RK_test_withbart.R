# preexisted data: rmap_bndry, lab_field.laea
library(BayesTree)

soil_property <- "P"
map_folder <- "/Users/jiehuachen/Documents/research/afsis/Ethiopia/datasets/predicted_maps/"

# estimate linear model
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <- log(lab_field.laea[soil_property]@data[[1]]+1)
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm), ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm), ]
Y_lm <- log(lab_field.laea[soil_property]@data[[1]]+1)



# prepare predictions covariates

load("Gtif_1k.RData")

# prepare prediction covariates with no missing data

predict_grid_1k_values <- predict_grid_1k@data[, -1]
predict_grid_1k_coords <- coordinates(predict_grid_1k_tif)[c(predict_grid_1k_tif@data$band1)==1&!is.na(predict_grid_1k_tif@data$band1), ]
predict_grid_1k_values.narm <- predict_grid_1k_values[!is.na(rowMeans(predict_grid_1k_values)), ]
predict_grid_1k_values.narm <- predict_grid_1k_values.narm[, colnames(X_lm)]
predict_grid_1k_values.narm <- as.matrix(predict_grid_1k_values.narm)
predict_grid_1k_coords <- predict_grid_1k_coords[!is.na(rowMeans(predict_grid_1k_values)), ]

bart.est <- bart(X_lm, Y_lm, x.test = predict_grid_1k_values.narm, sigquant=0.9, k=2, ntree=50, ndpost=100, nskip=100)

predict.bart.mean <-  bart.est$yhat.test.mean
predict.bart.sd <- apply(bart.est$yhat.test, 2, sd)


# add the residuals of BART into the dataframe

lab_field.laea@data["bart.resi"] <- Y_lm - bart.est$yhat.train.mean

## STEP TWO: variogram for residuals ____________________________________
#compute exprimental variogram of model residuals
vg.P.resi <- variogram(bart.resi ~ 1, data = lab_field.laea, cutoff=30E3)
plot(vg.P.resi)

#fit variogram model for residuals
vgmf.P.resi.sph<- fit.variogram(vg.P.resi, model = vgm(psill = 0.3, model = "Sph", range = 30E3, nugget = 1))
plot(vg.P.resi,vgmf.P.resi.sph, main="sph")
str(vgmf.P.resi.sph) #"SSErr"= num  1.26e-06

vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 30E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.sph, main="exp")
str(vgmf.P.resi.exp) #"SSErr"= num 1.25e-06  --> is the better fit model

vgmf.P.resi.gau<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Gau", range = 30E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.gau, main="gau")
str(vgmf.P.resi.gau) #"SSErr"= num 1.94e-06



#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.
predict_grid_1k_coords <- as.data.frame(predict_grid_1k_coords)
coordinates(predict_grid_1k_coords) =~x+y
predict_grid_1k_coords@proj4string@projargs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
krige.P.resi<- krige(bart.resi ~1, lab_field.laea, newdata=predict_grid_1k_coords, model = vgmf.P.resi.sph, nmax = 100, beta=0, na.action=na.omit)


#Copy the kriging results to the SpatialPixelsDataFrame map_prediction
krige.P.resi$pred.resi<-krige.P.resi$var1.pred
krige.P.resi$variance.resi<- krige.P.resi$var1.var
#________________________________________________


# STEP FOUR:prediction by adding regressional model prediction and interpolated residual __________________________

#Compute the final prediction by adding the regression model prediction and the interpoalted residual

predlogP <- predict.bart.mean +krige.P.resi$pred.resi
# change to original value
krige.P.resi$RKpred.P<-exp(predlogP)

# to understand the variance of prediction, we can see lower and upper limits of prediction as alternative of variance map (var1.var)
alpha<-0.05 # this is 95% confidence level.
loglower<-predlogP - qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(predict.bart.sd^2+ krige.P.resi$var1.var)
logupper<-predlogP + qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(predict.bart.sd^2+ krige.P.resi$var1.var)


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



