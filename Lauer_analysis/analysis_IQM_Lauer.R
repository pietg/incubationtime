	library(Rcpp)
	library(coarseDataTools)
	sourceCpp("incubationIQM.cpp")
	
	## data set test_data.csv, consisting of 181 triples. 
	## The triples consist of Exit time E and S_L and S_R, the left and right 		## boundaries of the interval for the time S of becoming symptomatic
	## We should have S_R > S_L
	 
	A<-read.table("data_L2.txt")
	
	## Nonparametric MLE, computed by incubationIQM.cpp for the 181 data triples
	## Chosen bandwidth for density estimate  is h=3.5
	
	output <- MLE(A,3.5)

	## List containing output
	B <- output$MLE
	C <- output$masses
	D <- output$SMLE
	E <- output$density
	
	## We only consider the density part
	x<-E[,1]
	y<-E[,2]
	
	## Lauer_data.txt is data file in the right format 
	## for using the function dic.fit in coarseDataTools
	
	A<-read.table("Lauer_data.txt")
	colnames(A) = c("EL","ER","SL","SR","type")
	dic.fit(A,start.par2=10,par1.int=c(1,4),par2.int=c(1,12),dist="W")
	dic.fit(A,dist="L")	
	
	## parameter estimates for Weibull density, coming from dic.fit
	a1 <- 2.453
	b1 <- 6.258
	
	## parameter estimates for log-normal density, coming from dic.fit
	a <- 1.621
	b <- 0.418
	
	## corresponding Weibull density
 	f1 <- function(x) {(a1/b1)*((x/b1)^(a1-1))*exp(-(x/b1)^a1)}
	x1 <-seq(0,max(x),by=0.1)
	y1<-f1(x1)
	
	## corresponding log-normal density
 	f2 <- function(x) {(1/(x*b*sqrt(2*pi)))*exp(-(log(x)-a)^2/(2*b^2))}
	x2 <-seq(0.1,max(x),by=0.1)
	y2<-f2(x2);
		
	## blue curve is nonparametric density estimate, 
	## black dotted curve the real incubation density 
	## red dashed curve is the Weibull estimate of the incubation density,
	## computed by the R package coarseDataTools
	
 	plot(c(-1000,-1000),xlim=c(min(x),max(x)), ylim=c(min(y,y1,y2),max(y,y1,y2)), main= "",ylab="",xlab="",bty="n",las=1)
	lines(x,y,lwd=2,col="blue")
  	lines(x1,y1,lwd=2,lty=2,col="red") 
  	lines(x2,y2,lwd=2,lty=3)
 


