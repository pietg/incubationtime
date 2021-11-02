	library(Rcpp)
	library(pracma)
	library(nloptr)
	sourceCpp("nonparMLE.cpp")
	source("Weibull_nloptr.R")
	source("lognormal_nloptr.R")

	NumIt = 1000
	n = 500
	p = 0.95
	
# data vectors
	S <- vector("numeric", n)
	E <- vector("numeric", n)
 
 #upper bound for exposure time 
	M <- 30
	
#parameters Weibull distribution
	a <- 3.035140901
	b <- 7.0897556

# input parameters for parametric distributions
	x=vector("numeric",2)
	alpha=vector("numeric",2)
  
	data <- matrix(0, nrow= n, ncol= 2, byrow = FALSE)
	MLEMat <- matrix(0, nrow= NumIt, ncol= 3, byrow = FALSE)
	colnames(MLEMat) <- c("MLE","Weibull","log-normal")


for (iter in 1: NumIt)
{
  	sim = iter
  	
  	print(iter)
   
	set.seed(sim)

	E <- runif(n,1,M)
	 
# generate data
	for (i in 1: n)
	{ 	
		y <- runif(1,0,E[i])
		#u <- runif(1,0,1)
  		#v <- b^(-1/a)*(log(1/(1-c*u)))^(1/a)
  		v <- rweibull(1,shape=a,scale=b)
		S[i] <- y+v
		data[i,1] <- max(0,S[i]-E[i])
		data[i,2] <- S[i]
	}

# Compute NPMLE	
	output <- NPMLE(n,data,bandwidth=6*n^(-1/5),percentile=0.95)
	percentile_NPMLE <- output$quantile

# provide options for parametric methods	
	opts <- list( "algorithm"= "NLOPT_LN_COBYLA","xtol_rel"= 1.0e-10, 
					"maxeval"= 10000, "tol_constraints_ineq" = rep(1.0e-10,2))

# x0 is initial value for Weibull fitting
	x0 <- c(3,7)
	
# x1 is initial value for Log-normal fitting
	x1 <- c(1,0.5)

# res1 is result for Weibull fitting	
	res1 <- nloptr(x0 = x0, eval_f = eval_f0,lb = lb0, ub = ub0, eval_g_ineq = eval_g_ineq0, opts = opts,data = data)

# res2 is result for Log-normal fitting        				
    res2 <- nloptr(x0 = x1, eval_f = eval_f1,lb = lb1, ub = ub1, eval_g_ineq = eval_g_ineq1, opts = opts,data = data)
 
 # parameter estimates for Weibull       
    a1 <- res1$solution[1]
    b1 <- res1$solution[2]

 # parameter estimates for log-normal     
    a2 <- res2$solution[1]
    b2 <- res2$solution[2]
	
	
 #	medians
		
	percentile_Weibull <- b1*(-log(1-p))^(1/a1)
	percentile_lognormal <- exp(a2+sqrt(2*b2^2)*erfinv(2*p-1))

# if one wants to have the results in a file, "decomment" the following lines
	
	#write(percentile_NPMLE,file = "median_MPMLE.txt", ncol =1, append = TRUE)
	#write(percentile_Weibull,file = "median_Weibull.txt", ncol =1, append = TRUE)
	#write(percentile_lognormal,file = "median_lognormal.txt", ncol =1, append = TRUE)

#make matrix for box plot	
	MLEMat[iter,] = c(percentile_NPMLE,percentile_Weibull,percentile_lognormal)

}

pdf("BoxPlot_95th_percentile.pdf")
boxplot(MLEMat,las=1)
abline(h=b*(-log(1-p))^(1/a),lwd=2,lty=1,col = "red")
dev.off()