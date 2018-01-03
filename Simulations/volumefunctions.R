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

fit.gamma <- function(V){
	#install.packages("nloptr")
	library(nloptr)
	library(actuar)
	GammaNLL <- function(pars, data){
		alpha <- pars[[1]]
		theta <- pars[[2]]
		return (-sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)))
	}

	shape<-NA
	scale<-NA
	fit <- nloptr(x0 = c(1, 4), eval_f = GammaNLL, lb = c(0,0), data = V,
		opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e5))
	shape <- fit$solution[1]
	scale<- fit$solution[2]
	return(c(shape,scale))
}
