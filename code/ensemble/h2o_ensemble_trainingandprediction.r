# H2O Ensemble Learning 
# Jiehua Chen 


## Model Training and Prediction

####library(devtools)
####install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")

### Load in h2oEnsemble library, and customized funcitons

library(h2oEnsemble) 
h2o.ensemble_predict <- function(object, newdata) {
  
  if (class(object) != "h2o.ensemble") {
    stop("object must be of class, h2o.ensemble")
  }
  if (!grepl("H2O", class(object$metafit))) {
    stop("H2O performance metrics are not supported for SuperLearner-based metalearners.")
  }
  
  # Training_frame may be a key or an H2OFrame object
  if ((!inherits(newdata, "Frame") && !inherits(newdata, "H2OFrame")))
    tryCatch(newdata <- h2o.getFrame(newdata),
             error = function(err) {
               stop("argument \"newdata\" must be a valid H2OFrame or id")
             })
  
  if (object$family == "binomial") {
    newdata_levelone <- h2o.cbind(sapply(object$basefits, function(ll) h2o.predict(object = ll, newdata = newdata)[,3]))
  } else {
    newdata_levelone <- h2o.cbind(sapply(object$basefits, function(ll) h2o.predict(object = ll, newdata = newdata)[,1]))
  }
  names(newdata_levelone) <- names(object$basefits)
  newdata_predictions <- as.data.frame(h2o.predict(object$metafit, newdata = newdata_levelone))

  return(newdata_predictions)
}


h2o.init(nthreads = -1,max_mem_size="5g")  # Start an H2O cluster with nthreads = num cores on your machine
h2o.removeAll() 

### load in all data
train <- h2o.importFile(path = normalizePath("./EthM3_Cal.csv"))
test <- h2o.importFile(path = normalizePath("./EthM3_Val.csv"))
data_whole <- h2o.rbind(train, test)


cov_names <- names(train)[32:length(names(train))]

#### y is the column name of the property for prediction, x is the array of covariate names
y <- "V0"
x <- cov_names

### Define base learners and meta learner

h2o.randomForest.1 <- function(..., ntrees = 50, nbins = 100, seed = 1) {
  h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
}

h2o.randomForest.2 <- function(..., ntrees = 20, nbins = 20, seed = 1) {
  h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
}

h2o.deeplearning.1 <- function(..., hidden = c(100,100), activation = "Rectifier", seed = 1) {
  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
}

h2o.deeplearning.2 <- function(..., hidden = c(200,200), activation = "Rectifier", seed = 1) {
  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
}

h2o.gbm.1 <- function(..., ntrees = 50, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
h2o.gbm.2 <- function(..., ntrees = 50, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)

learner <- c("h2o.randomForest.wrapper","h2o.gbm.wrapper", "h2o.deeplearning.wrapper", "h2o.randomForest.1", "h2o.randomForest.2","h2o.deeplearning.1", "h2o.deeplearning.2", "h2o.gbm.1", "h2o.gbm.2")
h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative,  lambda_search=TRUE)
metalearner <- "h2o.glm_nn"

### Model training

ensemble_model <- h2o.ensemble_training(data_whole, y, x, learner, metalearner)


### Load prediction covariates into h2o

library(raster)
library(rgdal)
grids <- stack(list.files("../../../../../remotesensing_data/ET_grids/", full.names = TRUE))

predictcov <-  getValues(grids)
coordniates_cov <- coordinates(grids)

predictcov_narm <- na.omit(predictcov)
coordinates_cov_narm <- coordniates_cov[!is.na(rowSums(predictcov)), ]

predictcov_h2o <- as.h2o(predictcov_narm, "newdata")


### Prediction

predictions_h2o <- h2o.ensemble_predict(ensemble_model$ensemblemodel, predictcov_h2o)

### Shut down h2o server when the computation is finished

h2o.shutdown()

### Write out tif file

predictions_results <- as.data.frame(cbind(coordinates_cov_narm, predictions=as.data.frame(predictions_h2o)))
coordinates(predictions_results) <- c("x", "y")
gridded(predictions_results) <- TRUE
writeGDAL (
           dataset=predictions_results["predict"],
           fname=paste("h2o_ensemble_",y,".tif", sep=""),
           drivername= "GTiff",
           type="Float32")

