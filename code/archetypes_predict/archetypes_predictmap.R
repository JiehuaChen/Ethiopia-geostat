### Ethiopian topsoil macro-nutrient composition & archetypal analyses
# Data provided by: Crop Nutrition Laboratory Services, Nairobi, http://www.cropnuts.com
# M. Walsh & J. Chen, May 2013

# Set local working directory, e.g.
# setwd("/Users/markuswalsh/Documents/Ethiopia/CNLS/Analyses")

### Load packages
# install.packages(c("downloader","arm","archetypes","compositions"), dependencies=T)
library(arm)
library(archetypes)
library(compositions)

### Load data
# download("https://www.dropbox.com/s/x7il9g863dwg1h8/CNLS_dat.zip", "CNLS_dat.zip", mode="wb")
# unzip("CNLS_dat.zip", overwrite=T)
labdata_folder <- "/Users/jiehuachen/Documents/research/afsis/Ethiopia/datasets/CNLS_dat/20130626"
setwd(labdata_folder)
samp <- read.csv("Samples.csv",  sep=",")
prof <- read.table("Profiles.csv", header=T, sep=",")
cnlsdat <- merge(prof, samp, by="Sample_ID")

### Macro-nutrient compositions and log ratio transforms
macronut <- c("N","P","K","S","Ca","Mg")
cdata <- as.data.frame(clr(cnlsdat[macronut]))
idata <- as.data.frame(clrInv(cdata))
pcplot(idata)

### Identification of compositional archetypes / endmembers
set.seed(20513)
cdata.arc <- stepArchetypes(data=cdata, k=1:6, nrep=5, verbose=T)
screeplot(cdata.arc)

# Select no. of archetypes, screeplot suggests k=7
cdata.arc7 <- bestModel(cdata.arc[[5]])
arcpar <- clrInv(parameters(cdata.arc7))
barplot(arcpar, names.arg=paste("Archetype", 1:5), legend.text=TRUE, space=0.5, horiz=T, axisnames=TRUE, args.legend=list(x = "topright"), xlim=c(0,1.2))

# predict archetypes (DA's)
setwd("/Users/jiehuachen/Documents/research/afsis/Ethiopia/map_results/map_v1/linear+rk_macronuts_mean_ppm_4717")
predict.N <- readGDAL("RKpred_N_linear_ppm_4717.tif")@data$band1
predict.P <- readGDAL("RKpred_P_linear_ppm_4717.tif")@data$band1
predict.K <- readGDAL("RKpred_K_linear_ppm_4717.tif")@data$band1
predict.S <- readGDAL("RKpred_S_linear_ppm_4717.tif")@data$band1
predict.Ca <- readGDAL("RKpred_Ca_linear_ppm_4717.tif")@data$band1
predict.Mg <- readGDAL("RKpred_Mg_linear_ppm_4717.tif")@data$band1


predict.macronut <- cbind(predict.N, predict.P, predict.K, predict.S, predict.Ca, predict.Mg)
predict.coords <- coordinates(readGDAL("RKpred_Mg_4717_ppm.tif"))
predict.coords <- predict.coords[!is.na(rowMeans(predict.macronut)), ]
predict.macronut <- predict.macronut[!is.na(rowMeans(predict.macronut)), ]


predict.cdata <- as.data.frame(clr(predict.macronut))
predict.archetypes <- predict(cdata.arc7, predict.cdata)

# get the archetype index for the largest archetype coefficients
predict.archetypes.max <- apply(predict.archetypes, 1, which.max)

# get the decimal representation of the archetype coefficients
predict.archetypes.weightmean <- (round(predict.archetypes*10)/10)%*%10^c(1:5)

archetype.predict.map <- SpatialPointsDataFrame(
  coords = predict.coords,
  data = data.frame(
 predict.archetype =predict.archetypes.max, 
 predict.archetype.mean =predict.archetypes.weightmean
	)
)

gridded(archetype.predict.map) <- TRUE
writeGDAL (
  dataset=archetype.predict.map["predict.archetype"],
  fname=paste("archtype_predicted_4717.tif", sep=""),
  drivername= "GTiff",
  type="Float32")

writeGDAL (
  dataset=archetype.predict.map["predict.archetype.mean"],
  fname=paste("archtype_mean_predicted_4717.tif", sep=""),
  drivername= "GTiff",
  type="Float32")


