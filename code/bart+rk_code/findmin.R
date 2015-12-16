findmin <- function(tree.cv, error.index){
	mean.sse.cv <- aggregate(tree.cv[,error.index], by=list(tree.cv[,1], tree.cv[,2], tree.cv[,3], tree.cv[,4], tree.cv[,5]), mean)
 
 min.mean.sse.cv.soil <- matrix(NA, length(unique(mean.sse.cv[,1])), dim(mean.sse.cv)[2])
 rownum <- 1
 for(i in unique(mean.sse.cv[,1])){
 	mean.sse.cv.soil <- mean.sse.cv[mean.sse.cv[,1]==i, ]
 	min.mean.sse.cv.soil[rownum, ] <-  as.matrix(mean.sse.cv.soil[which.min(mean.sse.cv.soil[,6]), ])
  rownum <- rownum + 1
 }
 if(error.index==7){
 	colnames(min.mean.sse.cv.soil) <-c("soil_prop", "sigdf", "sigquant", "k", "ntree", "mean_sse")
 }
 if(error.index==7){
 	colnames(min.mean.sse.cv.soil) <-c("soil_prop", "sigdf", "sigquant", "k", "ntree", "relative_error")
 } 
 return(min.mean.sse.cv.soil)
}

