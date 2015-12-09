soil.property = c("N","P","S","K","Ca","Mg","SOC", "pH")
sigdf.quant.cv= matrix(c(3, 3, 10, 0.90, 0.99, 0.75), 3, 2)
k.cv = c(1, 2, 3, 5)
ntree.cv = seq(10, 100, by=20)
ncv <- 50
random.seed = ceiling(runif(ncv, 0, 100000))

parameters.cv.list <- list(as.matrix(soil.property), as.matrix(sigdf.quant.cv), as.matrix(k.cv), as.matrix(ntree.cv), as.matrix(random.seed))
# combine all the paramters: length(soil.property)*dim(sigdf.quant.cv)[1]*length(k.cv)*length(ntree.cv)*length(random.seed)
source("colMeans_new.R")
total.dim <-  length(soil.property)*dim(sigdf.quant.cv)[1]*length(k.cv)*length(ntree.cv)*length(random.seed)

meshgrid <- function(parameters.cv.list){
	temp <- parameters.cv.list[[1]]
	for(i in 2:length(parameters.cv.list)){
		temp <- apply(temp,2,rep, dim(parameters.cv.list[[i]])[1])

		unique.values <- as.matrix(unique(temp))
		temp2 <- matrix(NA, dim(temp)[1], dim(parameters.cv.list[[i]])[2])
		for(j in 1:dim(unique.values)[1]){
			temp2[apply(matrix(apply(temp, 1, "==", unique.values[j, ]), dim(temp)[1], dim(temp)[2], byrow=TRUE), 1, mean)==1,] <- parameters.cv.list[[i]]
		}
		temp <- cbind(temp, temp2)
		print(i)
	}		
	return(temp)	
}

cv_parameters <- meshgrid(parameters.cv.list)
write.table(cv_parameters, "cv_parameters.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)

