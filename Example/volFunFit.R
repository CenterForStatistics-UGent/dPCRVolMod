########################
#                      #
# Author  M. Vynck     #
# Last update 03JAN-18 #
#                      #
########################

rejection.sampling <- function(density,n.obs){
	cutoff <- array(NA,n.obs)
	d.approx <- approxfun(density)
	unif.sample <- runif(n.obs,min(density$x),max(density$x))
	cutoff <- d.approx(unif.sample)
	max.y <- max(density$y)
	keep <- runif(n.obs,0,max.y)<=cutoff
	final.sample <- unif.sample[keep]
	return(final.sample)
}


estVolsAbs <- function(density, pos, n, estimate = TRUE, iter = 20){
		# density: a marginal volume distribution
		# pos: a vector of length one, specifying the number of positives
		#      in the target
		#
		# n: a vector of length one, specifying the number of partitions
		#      in the target
		#
		# estimate: whether to return estimates for the absolute quantity
		#           (TRUE) or to return the partition status and their
		#           associated volumes (FALSE)
		#
		# iter: number of iterations for obtaining unbiased estimates
		#       (default set to 20)		
		
		V<-rejection.sampling(density=density, n.obs=1000000)
		V<-V[V>0]
		v.mean <- mean(V)
		
		#we do not know the true concentration, so calculate an initial estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.init <- unname(-log(1-pos/n))
		C.init <- lambda.init/v.mean
		#cat(C.init,"\n")
		#loop until convergence

		#with this initial estimate, calculate the conditional distributions
		copies<-array(NA,length(V))	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]*C.init))
		}		
		copies.binary<-(copies>=1)*1
		df<-data.frame(copies.binary,V)
		colnames(df)<-c("Y","vol")

		n.pos <- round(mean(copies.binary)*length(V))
		n.neg <- round((1-mean(copies.binary))*length(V))
		vol.pos<-sample(df$vol[df$Y==1],n.pos,replace=TRUE)
		vol.neg<-sample(df$vol[df$Y==0],n.neg,replace=TRUE)
		
		##assign a volume to each positive/negative partition
		df.pos<-cbind(1,vol.pos[1:n.pos])
		df.neg<-cbind(0,vol.neg[1:n.neg])
		df<-rbind(df.pos,df.neg)
		df<-data.frame(df)
		colnames(df)<-c("Y","vol")
		
		#fit model using volume
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))
		C.new <- exp(coef(fit))
		df$vol<-v.mean
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))		
		bias <- exp(coef(fit)) - C.init
		C.new <- C.init - bias

		#simulate data with given volume distribution
		copies<-array(NA,length(V))
		new.lambda <- C.new*v.mean	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
		}
		#we do not know the true concentration, so calculate an estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.new <- unname(-log(table(copies)[1]/length(copies)))
		C.new.biased <- lambda.new/v.mean
		
		estimates <- C.init-bias+C.init-C.new.biased
		C.new.temp <- C.init-bias
		max.iter = 20
		iter<-1
		bias.count <- 1e07
		#while the initial biased concentration is not close to the estimated biased concentration
		#adjust C.new until convergence
		#when converged, C.new will be an unbiased estimator for C
		while(abs(C.init-C.new.biased)>1e-4&&iter<max.iter){
			iter<-iter+1
			C.new.temp <- C.new.temp+(C.init-C.new.biased)
			copies<-array(NA,length(V))
			new.lambda <- C.new.temp*v.mean	
			for(i in 1:length(V)){
				copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
			}
			#we do not know the true concentration, so calculate an estimate ignoring volume variability
			#this will be an underestimate of the true concentration
			lambda.new <- unname(-log(table(copies)[1]/length(copies)))
			C.new.biased <- lambda.new/v.mean
			estimates <- c(estimates,C.new.temp+C.init-C.new.biased)
			bias.count <- c(bias.count,C.init-C.new.biased)
			#cat("The difference is:",C.init-C.new.biased,"\t The estimated concentration is:",C.new.temp+C.init-C.new.biased,"\n")
		}
		if(iter<max.iter){
			#cat("Converged at concentration of:",C.new.temp+C.init-C.new.biased,"\n")
			C.best <- estimates[length(estimates)]
		} else {
			#cat("Not converged. Estimated concentration during optimization:",median(estimates[order(abs(bias.count))[1:10]]),"\n")	
			C.best <- 	median(estimates[order(abs(bias.count))[1:10]])	
		}


		#obtain final volumes and final estimate
		#with this initial estimate, calculate the conditional distributions
		copies<-array(NA,length(V))	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]*C.best))
		}		
		copies.binary<-(copies>=1)*1
		df<-data.frame(copies.binary,V)
		colnames(df)<-c("Y","vol")
		
		cut <- n
		n.pos <- round(mean(copies.binary)*cut)
		n.neg <- round((1-mean(copies.binary))*cut)
		vol.pos<-sample(df$vol[df$Y==1],n.pos,replace=TRUE)
		vol.neg<-sample(df$vol[df$Y==0],n.neg,replace=TRUE)
		
		##assign a volume to each positive/negative partition
		df.pos<-cbind(1,vol.pos[1:n.pos])
		df.neg<-cbind(0,vol.neg[1:n.neg])
		df<-rbind(df.pos,df.neg)
		df<-data.frame(df)
		colnames(df)<-c("Y","vol")
		
		#fit model using volume
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))
		C.final <- unname(exp(coef(fit)))
		C.final.var <- unname(as.numeric(vcov(fit)) * C.final^2)
		
		if(estimate){
			return(list(C.final,C.final.var))
		} else {
			return(df)
		}
}


estVolsRel <- function(density, pos, n, estimate = TRUE, iter = 20){
		# density: a marginal volume distribution
		# pos: a vector of length two, specifying the number of positives
		#      in the target and reference respectively
		#
		# n: a vector of length two, specifying the number of partitions
		#      in the target and reference respectively
		#
		# estimate: whether to return estimates for the relative quantity
		#           (TRUE) or to return the partition status and their
		#           associated volumes (FALSE)
		#
		# iter: number of iterations for obtaining unbiased estimates
		#       (default set to 20)		
		
		V<-rejection.sampling(density=density, n.obs=1000000)
		V<-V[V>0]

		# calculate average volume
		v.mean <- mean(V)
		
		#we do not know the true concentration, so calculate an initial estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.init <- unname(-log(1-pos[1]/n[1]))
		C.init <- lambda.init/v.mean
		
		lambda.init.norm <- unname(-log(1-pos[2]/n[2]))
		C.init.norm <- lambda.init.norm/v.mean

		#### RETRIEVE ORIGINAL CONCENTRATIONS
		#### FIRST FOR TARGET
		
		#with this initial estimate, calculate the conditional distributions
		copies<-array(NA,length(V))	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]*C.init))
		}		
		copies.binary<-(copies>=1)*1
		df<-data.frame(copies.binary,V)
		colnames(df)<-c("Y","vol")

		n.pos <- round(mean(copies.binary)*length(V))
		n.neg <- round((1-mean(copies.binary))*length(V))
		vol.pos<-sample(df$vol[df$Y==1],n.pos,replace=TRUE)
		vol.neg<-sample(df$vol[df$Y==0],n.neg,replace=TRUE)
		
		##assign a volume to each positive/negative partition
		df.pos<-cbind(1,vol.pos[1:n.pos])
		df.neg<-cbind(0,vol.neg[1:n.neg])
		df<-rbind(df.pos,df.neg)
		df<-data.frame(df)
		colnames(df)<-c("Y","vol")
		
		#fit model using volume
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))
		C.new <- exp(coef(fit))
		df$vol<-v.mean
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))		
		bias <- exp(coef(fit)) - C.init
		C.new <- C.init - bias

		#simulate data with given volume distribution
		copies<-array(NA,length(V))
		new.lambda <- C.new*v.mean	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
		}
		#we do not know the true concentration, so calculate an estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.new <- unname(-log(table(copies)[1]/length(copies)))
		C.new.biased <- lambda.new/v.mean
		
		estimates <- C.init-bias+C.init-C.new.biased
		C.new.temp <- C.init-bias
		max.iter = iter
		iter<-1
		bias.count <- 1e07
		#while the initial biased concentration is not close to the estimated biased concentration
		#adjust C.new until convergence
		#when converged, C.new will be an unbiased estimator for C
		while(abs(C.init-C.new.biased)>1e-4&&iter<max.iter){
			iter<-iter+1
			C.new.temp <- C.new.temp+(C.init-C.new.biased)
			copies<-array(NA,length(V))
			new.lambda <- C.new.temp*v.mean	
			for(i in 1:length(V)){
				copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
			}
			#we do not know the true concentration, so calculate an estimate ignoring volume variability
			#this will be an underestimate of the true concentration
			lambda.new <- unname(-log(table(copies)[1]/length(copies)))
			C.new.biased <- lambda.new/v.mean
			estimates <- c(estimates,C.new.temp+C.init-C.new.biased)
			bias.count <- c(bias.count,C.init-C.new.biased)
			#cat("The difference is:",C.init-C.new.biased,"\t The estimated concentration is:",C.new.temp+C.init-C.new.biased,"\n")
		}
		if(iter<max.iter){
			#cat("Converged at concentration of:",C.new.temp+C.init-C.new.biased,"\n")
			C.best <- estimates[length(estimates)]
		} else {
			#cat("Not converged. Estimated concentration during optimization:",median(estimates[order(abs(bias.count))[1:10]]),"\n")	
			C.best <- 	median(estimates[order(abs(bias.count))[1:10]])	
		}


		#### FOR REFERENCE

		#with this initial estimate, calculate the conditional distributions
		copies<-array(NA,length(V))	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]*C.init.norm))
		}		
		copies.binary<-(copies>=1)*1
		df<-data.frame(copies.binary,V)
		colnames(df)<-c("Y","vol")

		n.pos <- round(mean(copies.binary)*length(V))
		n.neg <- round((1-mean(copies.binary))*length(V))
		vol.pos<-sample(df$vol[df$Y==1],n.pos,replace=TRUE)
		vol.neg<-sample(df$vol[df$Y==0],n.neg,replace=TRUE)
		
		##assign a volume to each positive/negative partition
		df.pos<-cbind(1,vol.pos[1:n.pos])
		df.neg<-cbind(0,vol.neg[1:n.neg])
		df<-rbind(df.pos,df.neg)
		df<-data.frame(df)
		colnames(df)<-c("Y","vol")
		
		#fit model using volume
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))
		C.new <- exp(coef(fit))
		df$vol<-v.mean
		fit<-glm(Y~1+offset(log(vol)),data=df,family=binomial(cloglog))		
		bias <- exp(coef(fit)) - C.init.norm
		C.new <- C.init.norm - bias

		#simulate data with given volume distribution
		copies<-array(NA,length(V))
		new.lambda <- C.new*v.mean	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
		}
		#we do not know the true concentration, so calculate an estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.new <- unname(-log(table(copies)[1]/length(copies)))
		C.new.biased <- lambda.new/v.mean
		
		estimates <- C.init.norm-bias+C.init.norm-C.new.biased
		C.new.temp <- C.init.norm-bias
		max.iter = 20
		iter<-1
		bias.count <- 1e07
		#while the initial biased concentration is not close to the estimated biased concentration
		#adjust C.new until convergence
		#when converged, C.new will be an unbiased estimator for C
		while(abs(C.init.norm-C.new.biased)>1e-4&&iter<max.iter){
			iter<-iter+1
			C.new.temp <- C.new.temp+(C.init.norm-C.new.biased)
			copies<-array(NA,length(V))
			new.lambda <- C.new.temp*v.mean	
			for(i in 1:length(V)){
				copies[i]<-rpois(1,(V[i]/v.mean*new.lambda))
			}
			#we do not know the true concentration, so calculate an estimate ignoring volume variability
			#this will be an underestimate of the true concentration
			lambda.new <- unname(-log(table(copies)[1]/length(copies)))
			C.new.biased <- lambda.new/v.mean
			estimates <- c(estimates,C.new.temp+C.init.norm-C.new.biased)
			bias.count <- c(bias.count,C.init.norm-C.new.biased)
			#cat("The difference is:",C.init-C.new.biased,"\t The estimated concentration is:",C.new.temp+C.init-C.new.biased,"\n")
		}
		if(iter<max.iter){
			#cat("Converged at concentration of:",C.new.temp+C.init-C.new.biased,"\n")
			C.best.norm <- estimates[length(estimates)]
		} else {
			#cat("Not converged. Estimated concentration during optimization:",median(estimates[order(abs(bias.count))[1:10]]),"\n")	
			C.best.norm <- 	median(estimates[order(abs(bias.count))[1:10]])	
		}



		# calculate the conditional distributions
		copies<-array(NA,length(V))	
		copies.norm<-array(NA,length(V))	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]*C.best))
			copies.norm[i]<-rpois(1,(V[i]*C.best.norm))
		}		
		copies.binary<-(copies>=1)*1
		copies.binary.norm<-(copies.norm>=1)*1
		df<-data.frame(copies.binary,copies.binary.norm,V)
		colnames(df)<-c("Y.target","Y.norm","vol")

		cut <- max(n)
		n.pos.pos <- round(mean(((df$Y.target==1)+(df$Y.norm==1))>1)*cut)
		n.neg.neg <- round(mean(((df$Y.target==0)+(df$Y.norm==0))>1)*cut)
		n.pos.neg <- round(mean(((df$Y.target==1)+(df$Y.norm==0))>1)*cut)
		n.neg.pos <- round(mean(((df$Y.target==0)+(df$Y.norm==1))>1)*cut)
		n.total <- n.pos.pos+n.neg.neg+n.pos.neg+n.neg.pos
		vol.pos.pos<-sample(df$vol[((df$Y.target==1)+(df$Y.norm==1))>1],n.pos.pos,replace=TRUE)
		vol.neg.neg<-sample(df$vol[((df$Y.target==0)+(df$Y.norm==0))>1],n.neg.neg,replace=TRUE)
		vol.pos.neg<-sample(df$vol[((df$Y.target==1)+(df$Y.norm==0))>1],n.pos.neg,replace=TRUE)
		vol.neg.pos<-sample(df$vol[((df$Y.target==0)+(df$Y.norm==1))>1],n.neg.pos,replace=TRUE)

		
		##assign a volume to each positive/negative partition
		df.pos.pos<-cbind(1,1,vol.pos.pos)
		df.neg.neg<-cbind(0,0,vol.neg.neg)
		df.pos.neg<-cbind(1,0,vol.pos.neg)
		df.neg.pos<-cbind(0,1,vol.neg.pos)
		

		df<-rbind(df.pos.pos,df.pos.neg,df.neg.neg,df.neg.pos)
		df<-data.frame(df)
		colnames(df)<-c("Y.target","Y.norm","vol")
		df<-data.frame(Y=c(df$Y.target,df$Y.norm),vol=c(df$vol,df$vol),target=c(rep(TRUE,n.total),rep(FALSE,n.total)))
		
		#fit model using volume
		fit<-glm(Y~1+factor(target)+offset(log(vol)),data=df,family=binomial(cloglog))
		C.final <- unname(exp(coef(fit)[2]))
		C.final.var <- unname(as.numeric(vcov(fit)[2,2]) * C.final^2)
		
		if(estimate){
			return(list(C.final,C.final.var))
		} else {
			return(df)
		}
}


reduceVarAbs <- function(density, pos, n, estimate = TRUE, iter = 20, runs = 10){
		# density: a marginal volume distribution
		# pos: a vector of length one, specifying the number of positives
		#      in the target
		#
		# n: a vector of length one, specifying the number of partitions
		#      in the target
		#
		# estimate: whether to return estimates for the absolute quantity
		#           (TRUE) or to return the partition status and their
		#           associated volumes (FALSE)
		#
		# iter: number of iterations for obtaining unbiased estimates
		#       (default set to 20)		
		#
		# runs: number of times to perform the estimation
		
		conc <- c()
		conc.var <- c()
		for(i in 1:runs){
			est <- estVolsAbs(density=density, pos=pos, n=n, estimate = estimate, iter = iter)
			conc[i] <- est[[1]]
			conc.var[i] <- est[[2]]
		}
		index <- which(conc==sort(conc)[ceiling(runs/2)])
		if((runs/2)%%1==0){index<-c(index, which(conc==sort(conc)[ceiling(runs/2)+1]))}
		return(list(mean(conc[index]),mean(conc.var[index]),conc,conc.var))
}


reduceVarRel <- function(density, pos, n, estimate = TRUE, iter = 20, runs = 10){
		# density: a marginal volume distribution
		# pos: a vector of length one, specifying the number of positives
		#      in the target
		#
		# n: a vector of length one, specifying the number of partitions
		#      in the target
		#
		# estimate: whether to return estimates for the absolute quantity
		#           (TRUE) or to return the partition status and their
		#           associated volumes (FALSE)
		#
		# iter: number of iterations for obtaining unbiased estimates
		#       (default set to 20)		
		#
		# runs: number of times to perform the estimation
		
		conc <- c()
		conc.var <- c()
		for(i in 1:runs){
			cat(i)
			est <- estVolsRel(density=density, pos=pos, n=n, estimate = estimate, iter = iter)
			conc[i] <- est[[1]]
			conc.var[i] <- est[[2]]
		}
		index <- which(conc==sort(conc)[ceiling(runs/2)])
		if((runs/2)%%1==0){index<-c(index, which(conc==sort(conc)[ceiling(runs/2)+1]))}
		return(list(mean(conc[index]),mean(conc.var[index]),conc,conc.var))
}