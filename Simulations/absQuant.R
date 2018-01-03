########################
#                      #
# Author  M. Vynck     #
# Last update 03JAN-18 #
#                      #
########################

source('volumefunctions.R')
dvol<-read.csv("dropletvolumeDong.csv")

#set variability to (x*100)%
var.set<-0.20
dvol[,1]<-((dvol[,1]-mean(dvol[,1],na.rm=T))*var.set/sd(dvol[,1]/mean(dvol[,1],na.rm=T),na.rm=T)+mean(dvol[,1],na.rm=T))
#dvol[,2]<-((dvol[,2]-mean(dvol[,2],na.rm=T))*var.set/sd(dvol[,2]/mean(dvol[,2],na.rm=T),na.rm=T)+mean(dvol[,2],na.rm=T))
#dvol[,3]<-((dvol[,3]-mean(dvol[,3],na.rm=T))*var.set/sd(dvol[,3]/mean(dvol[,3],na.rm=T),na.rm=T)+mean(dvol[,3],na.rm=T))
#dvol[,4]<-((dvol[,4]-mean(dvol[,4],na.rm=T))*var.set/sd(dvol[,4]/mean(dvol[,4],na.rm=T),na.rm=T)+mean(dvol[,4],na.rm=T))

#d4<-density(dvol[,4],na.rm=T,bw="SJ",n=10000)
#d3<-density(dvol[,3],na.rm=T,bw="SJ",n=10000)
#d2<-density(dvol[,2],na.rm=T,bw="SJ",n=10000)
d1<-density(dvol[,1],na.rm=T,bw="SJ",n=10000)

V<-rejection.sampling(density=d1, n.obs=1000000)
V<-V[V>0]
cut<-15000
C<-5	
nsim <- 500
bias.sims <- matrix(100,nsim,3)
var.sims <- matrix(100,nsim,3)

#fit model using gamma distribution
fit.gamma<-fit.gamma(V)
shape<-fit.gamma[1]

for(z in 1:nsim){	
		cat("Iteration ",z," of ",nsim,".\n",sep="")
		#simulate data with given volume distribution
		v.mean <- mean(V)
		copies<-array(NA,length(V))
		true.lambda <- C*v.mean	
		for(i in 1:length(V)){
			copies[i]<-rpois(1,(V[i]/v.mean*true.lambda))
		}
		copies.binary<-(copies>=1)*1
		
		df.start <- data.frame(Y=copies.binary[1:cut],V=V[1:cut])
		
		fit<-glm(Y~1+offset(log(V)),data=df.start,family=binomial(cloglog))
		C.real <- exp(coef(fit))
		C.real.var <- as.numeric(vcov(fit)) * C.real^2
		df.start$V<-mean(df.start$V)
		fit<-glm(Y~1+offset(log(V)),data=df.start,family=binomial(cloglog))
		C.naive <- exp(coef(fit)) 	
		C.naive.var <- as.numeric(vcov(fit))* C.naive^2	
		
		lambda.gamma <- shape*(1/((1-mean(copies.binary[1:cut]))^(1/shape))-1)
		C.gamma<-lambda.gamma/mean(V[1:cut])
		C.gamma.var<- 1/(cut-1) * (1+lambda.gamma/shape)^2*((1+lambda.gamma/shape)^shape-1) /mean(V[1:cut])^2		
		C.gamma.var<-mean(df.start$Y)/(((1-mean(df.start$Y))^(2*(1/shape)+1))*(cut-1))/mean(df.start$V)^2



		
		#we do not know the true concentration, so calculate an initial estimate ignoring volume variability
		#this will be an underestimate of the true concentration
		lambda.init <- unname(-log(table(copies)[1]/length(copies)))
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
		C.final <- exp(coef(fit))
		C.final.var <- as.numeric(vcov(fit)) * C.final^2

		bias.sims[z,]<-(c(C.naive,C.gamma,C.final)-C.real)/C.real*100
		var.sims[z,]<-(c(C.naive.var,C.gamma.var,C.final.var)-C.real.var)/C.real.var*100
}	
	
save.image(paste("d1_",C,"_",var.set,".RData",sep=""))