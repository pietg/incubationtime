	library(Rcpp)
	sourceCpp("NPMLE_SR_single.cpp")
	
	#95% confidence intervals for 1000 samples of 1000 observations
	output <- CI_NPMLE()
	
	# the MLE for the last sample
	B <- output$MLE
	
	# the 95% confidennce intervals for the last sample
	C <- output$CI_MLE
	
	# coverages of the 1000 confidence intervals
	D <- output$percentages
	
	# the first column of E gives n times the variances over 
	#	the 1000 samples
	# the second column gives the diagonals of the inverse Fisher
	# 	information matrix:
	
	E <- output$Variances
	
	y1<-C[,1]
	y2<-C[,2]
	y3<-C[,3]
	y4<-C[,4]
	y5<-C[,5]
	
	plot(c(-10000,-10000),xlim=c(3,10), ylim=c(0.0,max(y3)), main= "", ylab="",xlab="",bty="n",las=1)
	
	points(y1,y2,col="red",pch = 20, cex = 1.5)
	points(y1,y3,col="black",pch = 20, cex = 1.5)
	segments(y1,y4,y1,y5,lwd=2)
	lines(y1,y2,lwd=2,lty=2,col="red")
	
	u1<-D[,1]
	v1<-D[,2]
	
	plot(c(-10000,-10000),xlim=c(3,10), ylim=c(0.5,1), main= "", ylab="",xlab="",bty="n",las=1)
	lines(u1,v1,lty=1,lwd=2)
	lines(c(min(u1),max(u1)),c(0.95,0.95),col="red",lty=2,lwd=2)
	
	
	
	