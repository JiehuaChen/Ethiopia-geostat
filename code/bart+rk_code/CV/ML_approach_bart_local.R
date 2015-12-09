### CV


# load library
library(BayesTree, quietly=TRUE)

# read in data
labdata <- read.csv("labdata.csv")
covdata <- read.csv("covdata.csv")
nfolder <- 5 # 5 folder CV

#if(line <- readLines(f,n=1, warn=FALSE)>0){

cv_parameters <- read.table("cv_parameters.txt", sep=",", stringsAsFactors=FALSE)

prediction.error <- rep(NA, dim(cv_parameters)[1])

for(k in 1:dim(cv_parameters[1:3,])[1]){
    print(k)
    # soil property to calibrate
    line <- cv_parameters[k,]
    # delete missing data
    Y <- labdata[, line[1,1]]
    X <- covdata
    Y.temp <- Y[!is.na(rowMeans(cbind(Y, X)))&Y>0]
    X <- X[!is.na(rowMeans(cbind(Y, X)))&Y>0, ]
    Y <- log(Y.temp)

    # cv parameters
    ntree.cv = as.numeric(line[1,5])
    nfolder <- 5 #5 fold cv 
    sigdf.quant= c(as.numeric(line[1,2]), as.numeric(line[1,3]))
    k.cv = as.numeric(line[1,4])

    # set random seed
    set.seed(line[1,6])
    sizetest <- ceiling(length(Y)/nfolder)
    testindex <- sort(sample(1:length(Y), size=sizetest))
    x.train <- X[-testindex, ]
    y.train <- Y[-testindex]
    x.test <- X[testindex, ]
    y.test <- Y[testindex]		

    # run bart	   
    temp.bart <- bart(x.train, y.train, x.test, sigdf = sigdf.quant[1], sigquant=sigdf.quant[2], k = k.cv, verbose=FALSE, usequant=TRUE, ntree = ntree.cv, ndpost=5000, keepevery=10, nskip=10000)
    prediction.error[k] <- sum((y.test - temp.bart$yhat.test.mean)^2)

}

write.table(cbind(cv_parameters, prediction.error), "prediction_errors_cv.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
