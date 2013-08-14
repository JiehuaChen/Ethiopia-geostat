# get cv results
source("/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/spatial.git/code/bart+rk_code/findmin.R")
source("/Users/jiehuachen/Documents/research/afsis/Ethiopia/git/spatial.git/code/bart+rk_code/relative_error.R")
tree25 <- read.table("ethiopia_bart_25_computed_messages.txt", sep=",")
labdata <- read.csv("data/labdata.csv", stringsAsFactors=FALSE)
covdata <- read.csv("data/covdata.csv",  stringsAsFactors=FALSE)

tree25_cv <- relative_error(labdata, covdata, tree25)
results25.sse <- findmin(tree25_cv, 7)
