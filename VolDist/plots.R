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
library(cowplot)
library(ggplot2)
########################
# EXPLORATORY ANALYSIS #
########################
d4<-density(dvol[,4],na.rm=T,bw=0.02,n=10000)
d3<-density(dvol[,3],na.rm=T,bw=0.02,n=10000)
d2<-density(dvol[,2],na.rm=T,bw=0.02,n=10000)
d1<-density(dvol[,1],na.rm=T,bw=0.02,n=10000)
pdf("density.pdf",width=3.25,height=2.5)
par(mai=c(0.6,0.55,0.03,0.03))
par(oma=c(0,0,0,0))
plot(d1,xlim=c(0.5,1.1),ylim=c(0,12), main="",xlab="",ylab="",cex.lab=.7, cex.axis=.7, cex.main=.7, cex.sub=.7)
mtext(side = 1, text = "Volume", line = 2, cex=.7)
mtext(side = 2, text = "Density", line = 2, cex=.7)
lines(d2,lty=2,col="blue")
lines(d3,lty=3,col="red")
lines(d4,lty=4,col="green")
legend("topleft",c("Sample 1","Sample 2","Sample 3","Sample 4"),col=c("black","blue","red","green"),lty=c(1:4),cex=0.7)
dev.off()

shape<-array(NA,4)
scale<-array(NA,4)
for(i in 1:4){
	fit <- unlist(fitdistr(dvol[,i][complete.cases(dvol[,i])],"gamma"))
    shape[i] <- fit[1]
    scale[i] <- 1/fit[2]
}
pdf("QQ.pdf",width=3.25,height=3.25)
par(mai=c(0.5,0.1,0.2,0.2))
par(omi=c(0,0.4,0.0,0))
par(mfrow=c(2,2))
for(i in 1:4){
	quantiles<-qgamma(seq(0,1,0.01),shape=shape[i],scale=scale[i])
	qqplot(quantiles,dvol[,i][complete.cases(dvol[,i])],xlab="Theoretical quantiles",ylab="Sample quantiles",main=paste("Sample",i,sep=" "),cex.lab=.7, cex.axis=.7, cex.main=.7, cex.sub=.7, cex=.7)
	mtext(side = 1, text = "Theoretical quantiles", line = 2, cex=.7)
	if(i%in%c(1,3)){mtext(side = 2, text = "Sample quantiles", line = 2, cex=.7)}
	qqline(dvol[,i][complete.cases(dvol[,i])], distribution= function(p) qgamma(p,shape=shape[i],scale=scale[i]))
}
dev.off()
plotlist <- list()

for(i in 1:4){
	sample <-dvol[,i][complete.cases(dvol[,i])]
 y <- quantile(sample, c(0.25, 0.75))
  x <- qgamma(c(0.25, 0.75),shape=shape[i],scale=scale[i])
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
	df <- data.frame(sample)
	plotlist[[i]]<-ggplot(data=df,aes(sample=sample))+stat_qq(distribution=qgamma,dparams=list(shape=shape[i],scale=scale[i]),size=0.5)+xlab("Theoretical quantiles")+ylab("Sample quantiles")+theme_bw()+ geom_abline(slope = slope, intercept = int)+theme(axis.title=element_text(size=8),axis.text=element_text(size=7))

}
pdf("QQ.pdf",width=3.25,height=3.25)
plot_grid(plotlist=plotlist,labels=LETTERS[1:4])
dev.off()