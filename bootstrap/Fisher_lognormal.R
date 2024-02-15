	library(pracma)
	library(nloptr)
	library(Deriv)
	library(MASS)
	source("lognormal_nloptr.R")
	source("lognormal_Fisher.R")
	
	NumIt =1
	n=500
	ngrid =100
	
	#data vectors
	S <- vector("numeric", n)
	E <- vector("numeric", n)

	coverage <- vector("numeric", ngrid)
	length_interval <- vector("numeric", ngrid)
	sigma <- vector("numeric", ngrid)
	
	f2 <- vector("numeric", ngrid)
	upperbound <- vector("numeric", ngrid)
	lowerbound <- vector("numeric", ngrid)

 #upper bound for exposure time 
	M <- 30
 #upper bound for incubation time 
	M1 <- 15
	
	grid  <- vector("numeric", ngrid)
	for (i in 1: ngrid)
		grid[i] = M1*i/ngrid
		
	#parameters Weibull distribution
	
	a <- 3.035140901
	b <- 7.0897556
	
	#parameters log-normal distribution
	c <- 1.763329
	d <- 0.3728888

	# provide options for parametric methods	
	opts <- list( "algorithm"= "NLOPT_LN_COBYLA","xtol_rel"= 1.0e-10, 
					"maxeval"= 10000, "tol_constraints_ineq" = rep(1.0e-10,2))

	Fisher_matrix <- matrix(0, nrow= 2, ncol= 2, byrow = FALSE)
	Inverse_Fisher_matrix <- matrix(0, nrow= 2, ncol= 2, byrow = FALSE)
	data <- matrix(0, nrow= n, ncol= 2, byrow = FALSE)
	
	# u1 is initial value for lognormal fitting
	u1 <- c(1,0.5)
	
	for (j in 1:ngrid)
	  coverage[j] <- 0
	  
	f <- function(x) {plnorm(x, meanlog = c, sdlog = d, lower.tail = TRUE, log.p = FALSE)}
	#f <- function(x) {1-exp(-(x/b)^a)}	
	u0 <-seq(0,M1,by=0.01)
	v0<-f(u0)
	
for (iter in 1: NumIt)
{
  	sim = NumIt+iter
  	
  	print(iter)
  	
  	for (i in 1:2)
  	{
  	  for (j in 1:2)
  	    Fisher_matrix[i,j]=0
  	}

   
	set.seed(sim)

	E <- runif(n,1,M)
	# generate data
	for (i in 1: n)
	{ 	
		y <- runif(1,0,E[i])
		#u <- runif(1,0,1)
  		#v <- b^(-1/a)*(log(1/(1-c*u)))^(1/a)
  		#v <- rweibull(1,shape=a,scale=b)
  		v <- rlnorm(1, meanlog =c, sdlog = d)
		S[i] <- y+v
		data[i,1] <- max(0,S[i]-E[i])
		data[i,2] <- S[i]
		
	}
	
	# res2 is result for lognormal fitting	
	res2 <- nloptr(x0 = u1, eval_f = eval_f1,lb = lb1, ub = ub1, 
	               eval_g_ineq = eval_g_ineq1, opts = opts,data = data)
	
	#parameter estimates for Weibull       
	c1 <- res2$solution[1]
	c2 <- res2$solution[2]
    
    for (i in 1: n)
	{
		if (S[i]<=E[i])
		{
			Fisher_matrix[1,1] <- Fisher_matrix[1,1] + g11(S[i],c1,c2)
			Fisher_matrix[2,2] <- Fisher_matrix[2,2] + g22(S[i],c1,c2)
			Fisher_matrix[1,2] <- Fisher_matrix[1,2] + g12(S[i],c1,c2)
			Fisher_matrix[2,1] <- Fisher_matrix[1,2]
		}
		else
		{
			Fisher_matrix[1,1] <- Fisher_matrix[1,1] + h11(E[i],S[i],c1,c2)
			Fisher_matrix[2,2] <- Fisher_matrix[2,2] + h22(E[i],S[i],c1,c2)
			Fisher_matrix[1,2] <- Fisher_matrix[1,2] + h12(E[i],S[i],c1,c2)
			Fisher_matrix[2,1] <- Fisher_matrix[1,2]
		}
	}
	
	Inverse_Fisher_matrix <- n*inv(Fisher_matrix)
	
	a1 <- Inverse_Fisher_matrix[1,1]
	a2 <- Inverse_Fisher_matrix[2,2]
	a3 <- Inverse_Fisher_matrix[1,2]
	
	for (i in 1: ngrid)
	{
		x <- grid[i]
		z1 = -exp(-(log(x)-c1)^2/(2*c2^2))/(c2*sqrt(2*pi))
		z2 = (c1-log(x))*exp(-(log(x)-c1)^2/(2*c2^2))/(c2^2*sqrt(2*pi))
		sigma[i] <- sqrt(a1*z1^2+a2*z2^2+2*a3*z1*z2)
	}

	
	for (j in 1:ngrid)
	{
	  x <- grid[j]
	  f2[j] <- plnorm(x, meanlog = c1, sdlog = c2, lower.tail = TRUE, log.p = FALSE)
	  upperbound[j] <- f2[j] + 1.96*sigma[j]/sqrt(n)
	  lowerbound[j] <- f2[j] - 1.96*sigma[j]/sqrt(n)
	  if (f(x) < lowerbound[j]|| f(x) >  upperbound[j])
	    coverage[j] = coverage[j] +1
	  
	  length_interval[j] <- upperbound[j]-lowerbound[j]
	}
	write(length_interval,file = "length_interval_parametric_lognormal_lognormal.txt", ncol =ngrid, append = TRUE)	
}
	coverage <- 1-coverage/NumIt
	
	write(coverage,file = "coverage_lognormal_lognormal.txt", ncol =ngrid, append = TRUE)
	
	pdf("coverages.pdf")	
	plot(c(-10000,-10000),xlim=c(3,10), ylim=c(0.0,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(grid,coverage,lty=1,lwd=2)
	lines(c(min(grid),max(grid)),c(0.95,0.95),col="red",lty=2,lwd=2)
	dev.off()
	
	x1 <- grid
	y1 <- f2
	y2 <- lowerbound
	y3 <- upperbound
	
	plot(c(-100,-100),xlim=c(0,M1), ylim=c(0,1), main= "",ylab="",xlab="",bty="n",las=1)
	lines(u0,v0,lwd=2,col="red")
	lines(x1,y1,lwd=2,col="blue")
	#lines(x1,y2,lwd=2)
	#lines(x1,y3,lwd=2)
	segments(x1,y2,x1,y3)
	
	