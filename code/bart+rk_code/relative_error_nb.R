relative_error<- function(totaldata, tree.cv){
	unique_randomseed <- unique(tree.cv[,6])
	unique_soil_property <- unique(as.character(tree.cv[,1]))
	nfolder <- 5
	var_y <- rep(NA, 3)
	for(i in 1:length(unique_soil_property)){
		Y <- totaldata[, unique_soil_property[i]]
		X <- totaldata[, 18:45]
		Y.temp <- Y[!is.na(rowMeans(cbind(Y, X)))]
		X <- X[!is.na(rowMeans(cbind(Y, X))), ]
		Y <- Y.temp
		for(j in 1:length(unique_randomseed)){
			set.seed(unique_randomseed[j])
			sizetest <- ceiling(length(Y)/nfolder)
			testindex <- sort(sample(1:length(Y), size=sizetest))
			x.test <- X[testindex, ]
			y.test <- Y[testindex]		
			var_y <- rbind(var_y, c(unique_soil_property[i], unique_randomseed[j], var(y.test)*sizetest))
		}
	}
	var_y <- var_y[-1, ]
	relative_error <- rep(NA, dim(tree.cv)[1])
	for(i in 1:dim(tree.cv)[1]){
		var_y_one <- var_y[as.character(tree.cv[i,1])==var_y[,1]&tree.cv[i,6]==var_y[,2],3]
		relative_error[i] <- tree.cv[i,7]/as.numeric(var_y_one)
	}
	return(cbind(tree.cv, relative_error))	
}
