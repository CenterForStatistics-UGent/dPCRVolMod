########################
#                      #
# Author  M. Vynck     #
# Last update 03JAN-18 #
#                      #
########################

vols <- read.csv("dropletvolumeDong.csv")

ks.matrix <- matrix(NA,3,3)
for(i in 1:3){
	for(j in (i+1):4){
		ks.matrix[i,j-1]<-ks.test(vols[,i][complete.cases(vols[,i])],vols[,j][complete.cases(vols[,j])])$p.value
	}
}
p.adjust(as.vector(ks.matrix)[c(1,4,5,7,8,9)],method="holm")