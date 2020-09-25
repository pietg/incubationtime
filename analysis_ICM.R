	library(Rcpp)
	sourceCpp("Weibull.cpp")
	sourceCpp("NPMLE_ICM.cpp")

	A<-read.table("inputdata_Wuhan.txt")
	output1 <- NPMLE(A)
	output2 <- Weibull(A)

	B1 <- output1$MLE
	C1 <- output1$SMLE
	D1 <- output1$dens
	
	B2 <- output2$df
	C2 <- output2$dens

	x<-B1[,1]
	y<-B1[,2]
	x1<-B2[,1]	
	y1<-B2[,2]
	
	t<-D1[,1]
	u<-D1[,2]
	t1<-C2[,1]	
	u1<-C2[,2]	
plot(c(-1000,-1000),xlim=c(0,14),ylim=c(0,max(u,u1)), main= "",ylab="",xlab="",bty="n",las=1)

	#lines(x,y,lwd=2,col="red",type='s')
	#lines(x1,y1,lwd=2,col="blue")
	
	lines(t,u,lwd=2, lty=1,col="blue")
	lines(t1,u1,lwd=2, lty=2,col="red")
 


