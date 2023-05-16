	library(Rcpp)
	A<-read.table("inputdata_Wuhan.txt") 
	sourceCpp("bandwidth_df_choice.cpp")
	
	output <- Compute_bandwidth(A)
	
	B <- output$SMLE
	C <- output$MSE
   	x1<-B[,1]
   	y1<-B[,2]
   	
   	x2<-C[,1]
   	y2<-C[,2]
   
    a<-3.035140901
	b<-0.002619475
	f <- function(x) {1-exp(-b*x^a)/(1-exp(-b*20^a))}	
	x0 <-seq(0,20,by=0.01)
	y0<-f(x0)

output$bandwidth

pdf("SMLE.pdf")
	plot(c(-100,-100),xlim=c(0,20), ylim=c(0,max(y0,y1)), main= "",ylab="",xlab="",bty="n",las=1)
	lines(x1,y1,lwd=2,col="blue")
	lines(x0,y0,lwd=2,lty=2,col="red")
dev.off()

pdf("MSE.pdf")
	plot(c(-100,-100),xlim=c(min(x2),max(x2)), ylim=c(min(y2),max(y2)), main= "",ylab="",xlab="",bty="n",las=1)
	lines(x2,y2,lwd=2,col="blue")
dev.off()	
	