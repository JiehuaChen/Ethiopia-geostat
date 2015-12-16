# get cv results
source("findmin.R")
source("relative_error.R")
tree25 <- read.table("ethiopia_bart_2015_sp_10-100_computed_messages.txt", sep=",")
totaldata <- read.csv("etm3.csv",  stringsAsFactors=FALSE)

tree25_cv <- relative_error(totaldata, tree25)
results25.sse <- findmin(tree25_cv, 8)
