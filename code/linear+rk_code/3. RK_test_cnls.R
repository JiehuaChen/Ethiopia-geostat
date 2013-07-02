library(gstat)

# preexisted data: predict_grid_1k, lab_field.laea

soil_property <- "P"
map_folder <- "~/ethiosis/map_results/"

# estimate linear model
covariates.names <- do.call("rbind",strsplit(grid.list, split=".tif"))
X_lm <- t(do.call("rbind", lab_field.laea@data[, names(lab_field.laea@data)%in%covariates.names]))
Y_lm <- log(lab_field.laea[soil_property]@data[[1]])
lab_field.laea <- lab_field.laea[!is.na(rowMeans(X_lm)+Y_lm), ]
X_lm <- X_lm[!is.na(rowMeans(X_lm)+Y_lm), ]
Y_lm <- log(lab_field.laea[soil_property]@data[[1]])

lm.P<-lm(Y_lm~ X_lm + x + y, data=lab_field.laea)
summary(lm.P)
stepAIC(lm.P)

# stepAIC selected covariates for Boron are then
variable_included <- c("bio1", "bio12")

lm.P2 <-lm(Y_lm ~ bio1 + bio12+ x + y, data=lab_field.laea) # stepAIC drops  rfl_blue_1
bestmodel<-lm.P2  # excluding rfl_blue_1_new1
summary(bestmodel)


#make data frame with candidate predictors (full model), and soil variable of interest
d <- subset(x=lab_field.laea, select = c(variable_included, "x", "y", soil_property))
d <- na.omit(d) # this can drop the NA values
#Add model residuals to d
d$P.fit <- predict(bestmodel, data=d)
d$P.resi <- Y_lm-d$P.fit
lab_field.laea$P.resi<-d$P.resi

#add projection attributes to d (saved in beginning of this script)

proj4string(d) <- (lab_field.laea)@proj4string@projargs
proj4string(rmap_bndry) <- (lab_field.laea)@proj4string@projargs

#remove points with exactly same coordinates
ids<-zerodist(d)
if(dim(ids)[1]>0){
	d<-d[--c(ids[,1], ids[, 2]),]
}

## STEP TWO: variogram for residuals ____________________________________
#compute exprimental variogram of model residuals
vg.P.resi <- variogram(P.resi ~ 1, data = d, cutoff=100E3)
plot(vg.P.resi)

#fit variogram model for residuals
vgmf.P.resi.sph<- fit.variogram(vg.P.resi, model = vgm(psill = 0.02, model = "Sph", range = 30E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.sph, main="sph")
str(vgmf.P.resi.sph) #"SSErr"= num  1.26e-06

vgmf.P.resi.exp<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Exp", range = 30E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.sph, main="exp")
str(vgmf.P.resi.exp) #"SSErr"= num 1.25e-06  --> is the better fit model

vgmf.P.resi.gau<- fit.variogram(vg.P.resi, model = vgm(psill =0.02, model = "Gau", range = 30E3, nugget = 0.01))
plot(vg.P.resi,vgmf.P.resi.gau, main="gau")
str(vgmf.P.resi.gau) #"SSErr"= num 1.94e-06

#____________________________________________________


# STEP THREE: predict P by using fitted the model (the bestmodel) selected from the 
# stepAIC _____________________


lmresult.P<-predict.lm(bestmodel, newdata =rmap_bndry, se.fit = TRUE, data=d) # se.fit is the confidence interval 
str(lmresult.P)
lmresult.P <- na.omit(lmresult.P)

#sidat1.ov$lmpred.P<-lmresult.P$fit
#spplot(lmpred.P, main="lmPred")

#head(sidat1.ov)
#View(sidat1.ov)

#Now interpolate model residuals by simple kriging to the nodes of the grid.
#Simple kriging means that we know the model mean (mu), which is in this case 0.


krige.P.resi<- krige(P.resi ~1, d, newdata=rmap_bndry, model = vgmf.P.resi.exp, nmax = 100, beta=0, na.action=na.omit)
# beta=0 means simple kriging since we know the mean, we donot assume mean is constant for residuals.
# spplot(krige.P.resi, main= "krige.P.resi")


#Copy the kriging results to the SpatialPixelsDataFrame map_prediction
krige.P.resi$pred.resi<-krige.P.resi$var1.pred
krige.P.resi$variance.resi<- krige.P.resi$var1.var
#________________________________________________


# STEP FOUR:prediction by adding regressional model prediction and interpolated residual __________________________

#Compute the final prediction by adding the regression model prediction and the interpoalted residual

predlogP <- lmresult.P$fit +krige.P.resi$pred.resi


# to plot
# ggplot()+  geom_raster(data=sidat1.ov, mapping=aes(x=x, y=y, fill=predlogP))+  coord_equal()

# change to original value
krige.P.resi$RKpred.P<-exp(predlogP)
#sidat1.ov$RKpred.P<-exp(predlogP)

# to understand the variance of prediction, we can see lower and upper limits of prediction as alternative of variance map (var1.var)
alpha<-0.05 # this is 95% confidence level.
loglower<-predlogP - qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(lmresult.P$se.fit^2+ krige.P.resi$var1.var)
logupper<-predlogP + qnorm(p=1-alpha/2,mean=0,sd=1)*sqrt(lmresult.P$se.fit^2+ krige.P.resi$var1.var)


#predlogP$RK.var.P<-exp(predlogP$variance.resi)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=x, y=y, fill=RKpred.P))+ coord_equal() # fill mean give the value to the cell with color attached to it

#compute the variance of the regression+kriging prediction error by adding the regression prediction error variance and the kriging variace of the residuals
krige.P.resi$RKlower<-exp(loglower)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=s1, y=s2, fill=RKlower))+ coord_equal()

krige.P.resi$RKupper<-exp(logupper)
#ggplot()+ geom_raster(data=as.data.frame(predlogP), mapping=aes(x=s1, y=s2, fill=RKupper))+ coord_equal()

setwd(map_folder)

# name convensions: (prediction method)_(property predicted)_(statistics predicted)_(unit of the predictions)_(version of data which the predictions are based on).tif

predlogP.new <- krige.P.resi
gridded(predlogP.new)<-TRUE
writeGDAL (
  dataset=predlogP.new["RKupper"],
  fname=paste("linear+rk.",soil_property,"_new.tif", sep=""),
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



