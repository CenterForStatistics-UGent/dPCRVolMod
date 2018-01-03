########################
#                      #
# Author  M. Vynck     #
# Last update 03JAN-18 #
#                      #
########################

dvol<-read.csv("dropletvolumeDong.csv")
source('volumefunctions.R')
library(MASS)
library(actuar)

########################
# EXPLORATORY ANALYSIS #
########################
d4<-density(dvol[,4],na.rm=T,bw=0.02,n=10000)
d3<-density(dvol[,3],na.rm=T,bw=0.02,n=10000)
d2<-density(dvol[,2],na.rm=T,bw=0.02,n=10000)
d1<-density(dvol[,1],na.rm=T,bw=0.02,n=10000)
plot(d1,xlim=c(0.5,1.1),ylim=c(0,14), main="Densities of droplet volumes for 4 channels",xlab="Volume (nl)",ylab="Density")
lines(d2,lty=2,col="blue")
lines(d3,lty=3,col="red")
lines(d4,lty=4,col="green")
legend("topleft",c("Sample 1","Sample 2","Sample 3","Sample 4"),col=c("black","blue","red","green"),lty=c(1:4))

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(5,4,2,1))
for(i in 1:4){
	quantiles<-qgamma(seq(0,1,0.01),shape=shape[i],scale=scale[i])
	qqplot(quantiles,dvol[,i][complete.cases(dvol[,i])],xlab="Theoertical quantiles",ylab="Sample quantiles",main=paste("Sample",i,sep=" "))
	qqline(dvol[,i][complete.cases(dvol[,i])], distribution= function(p) qgamma(p,shape=shape[i],scale=scale[i]))
}