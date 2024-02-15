library(Rcpp)
sourceCpp("bootstrap_SMLE.cpp")

	NumIt =100
	n=500
	ngrid =100
	
	#data vectors
	S <- vector("numeric", n)
	E <- vector("numeric", n)

	percentage  <- vector("numeric", ngrid)
	coverage <- vector("numeric", ngrid)
	length_interval <- vector("numeric", ngrid)
	lowerbound <- vector("numeric", ngrid)
	upperbound <- vector("numeric", ngrid)

 #upper bound for exposure time 
	M <- 30
 #upper bound for incubation time 
	M1 <- 15
	
	grid  <- vector("numeric", ngrid)
	for (i in 1: ngrid)
		grid[i] = M1*i/ngrid

#parameters Weibull distribution
	#a<-3.035140901
	#b<-0.002619475
	#c<- 1-exp(-b*M1^a)

	#parameters Weibull distribution	
	a <- 3.035140901
	b <- 7.0897556
	
	#parameters log-normal distribution
	c <- 1.763329
	d <- 0.3728888
	
	data <- matrix(0, nrow= n, ncol= 2, byrow = FALSE)
	
	for (j in 1:ngrid)
	  coverage[j] <- 0
	  
	 f <- function(x) {plnorm(x, meanlog = c, sdlog = d, lower.tail = TRUE, log.p = FALSE)}
	 #f <- function(x){pweibull(x, shape = a, scale = b, lower.tail = TRUE, log.p = FALSE)}

	
for (iter in 1: NumIt)
{
  sim = iter
   
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

# Compute NPMLE	
	output <- ComputeIntervals_df(n,M1,ngrid,data)
	length_interval = output$length_interval
	
	write(length_interval,file = "length_intervals_SMLE.txt", ncol =ngrid, append = TRUE)
	
	for (j in 1:ngrid)
	{
		x <- grid[j]
		lowerbound[j] <- output$CI_df[j,3]
		upperbound[j] <- output$CI_df[j,4]
		if (f(x) < lowerbound[j] || f(x) > upperbound[j])
			coverage[j] = coverage[j] +1
	}	
	
	print(c(iter,coverage[ngrid/2]))
}

coverage <- 1-coverage/NumIt

write(coverage,file = "coverage_SMLE.txt", ncol =ngrid, append = TRUE)

B <- output$CI_df
   x1<-B[,1]
   y1<-B[,2]
   y2<-B[,3]
   y3<-B[,4]

	x0 <-seq(0,M1,by=0.01)
	y0<-f(x0)

	plot(c(-100,-100),xlim=c(0,M1), ylim=c(0,1), main= "",ylab="",xlab="",bty="n",las=1)
	lines(x0,y0,lwd=2,col="red")
	lines(x1,y1,lwd=2,col="blue")
	#lines(x1,y2,lwd=2)
	#lines(x1,y3,lwd=2)
	segments(x1,y2,x1,y3)
	
	pdf("coverages.pdf")
	
	plot(c(-10000,-10000),xlim=c(3,10), ylim=c(0.0,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(grid,coverage,lty=1,lwd=2)
	lines(c(min(grid),max(grid)),c(0.95,0.95),col="red",lty=2,lwd=2)
	
	dev.off()
