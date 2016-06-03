# H2O Ensemble Learning 
# Jiehua Chen 



## Model Calibration and Selection

### install h2oEnsemble package if not installed

####library(devtools)
####install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")


### Load in h2oEnsemble library, and customized funcitons
library(h2oEnsemble) 
h2o_ensemble_predictionerror <- function(train, test, y, x, learner, metalearner){

    if(is.factor(train[,y])){
        family="binomial"
    }else{
        family="gaussian"
    }

    fit <- h2o.ensemble(x = x, y = y, 
                    training_frame = train, 
                    family = family, 
                    learner = learner, 
                    metalearner = metalearner,
                    cvControl = list(V = 5))
    
    newperf <- h2o.ensemble_performance(fit, newdata = test)
    return(newperf)
}

### Start h2o server
#### max_mem_size sets the maximal memory allocated for h2o


h2o.init(nthreads = -1, max_mem_size="5g")  # Start an H2O cluster with nthreads = num cores on your machine, -1 means using all cores
h2o.removeAll() 

### Load in data
train <- h2o.importFile(path = normalizePath("./EthM3_Cal.csv"))
test <- h2o.importFile(path = normalizePath("./EthM3_Val.csv"))

cov_names <- names(train)[32:length(names(train))]

#### y is the column name of the property for prediction, x is the array of covariate names
##### logarithm transformation if needed

y <- "P"
x <- cov_names

train[,y] <- log(train[,y])
test[,y] <- log(test[,y])


### Define learner and meta learner. 
#### . There are many options for setting up learner and metalearner, so we should use this script try different options and get the best one;
#### . You can see that there are a lot of flexibilities for defining those learners;

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


### Ensemble learning wrapper to implement the diagram procedure
#### check out ensembling_onelevel.pdf

h2o_ensemble_predictionerror(train, test, y, x, learner, metalearner)


### Shut down h2o server when the computation is finished

h2o.shutdown()


