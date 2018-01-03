########################
#                      #
# Author  M. Vynck     #
# Last update 03JAN-18 #
#                      #
########################

# source needed functions
source("volFunFit.R")

                      
###################################
# example absolute quantification #
###################################

# set seed for reproducibility of examples
set.seed(1)

# suppose we have volumes from a normal distribution
# with mean 0.85nl and sd 0.17nl (i.e. an RSD of 20%)
vols <- rnorm(200, 0.85, 0.17)

# we have obtained the following number of positives:
pos.count <- 15890
n.count <- 16500

# without volume correction this corresponds to
# a concentration of
-log(1-pos.count/n.count)/0.85

# i.e. 3.88 copies per nl

# with volume correction this corresponds to
# a concentration (and associated standard error) of
estVolsAbs(density(vols,bw="SJ"), pos.count, n.count)

# i.e. 4.11 copies per nl
# indeed, ignoring volume variability leads
# to an underestimation!



###################################
# example relative quantification #
###################################

# set seed for reproducibility of examples
set.seed(1)

# we have obtained the following number of positives:
pos.count <- c(15890,10574)
n.count <- c(16500,16241)

# without volume correction this corresponds to
# a relative quantity of
(-log(1-pos.count[1]/n.count[1]))/(-log(1-pos.count[2]/n.count[2]))

# i.e. 3.13

# with volume correction this corresponds to
# a relative quantity (and associated standard error) of
estVolsRel(density(vols,bw="SJ"), pos.count, n.count)

# i.e. 3.31
# indeed, ignoring volume variability leads
# to an underestimation!


#################################
# reducing sampling variability #
#################################

# set seed for reproducibility of examples
set.seed(1)

# the estimator is sampling based and introduces some
# sampling variability, which may be significantly
# reduced by running multiple iterations of the estimation
# procedure, at the expense of increased computation time

# we have obtained the following number of positives:
pos.count <- 15890
n.count <- 16500

redEst <- reduceVarAbs(density(vols,bw="SJ"), pos.count, n.count, runs=10)
redEst[[1]] # concentration estimate
redEst[[2]] # std error of concentration estimate

# similar use of the function reduceVarRel for
# relative quantification