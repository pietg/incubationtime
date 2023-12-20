library(Rcpp)
A<-read.table("data.txt") 
sourceCpp("bootstrap_SMLE.cpp")
output <- ComputeIntervals_df(A)
B <- output$CI_df
   x1<-B[,1]
   y1<-B[,2]
   y2<-B[,3]
   y3<-B[,4]
   
    a<-3.035140901
	b<-0.002619475
	f <- function(x) {(1-exp(-b*x^a))/(1-exp(-b*20^a))}	
	x0 <-seq(0,15,by=0.01)
	y0<-f(x0)

pdf("CI_SMLE.pdf")
plot(c(-100,-100),xlim=c(0,15), ylim=c(0,1), main= "",ylab="",xlab="",bty="n",las=1)
lines(x0,y0,lwd=2,col="red")
segments(x1,y2,x1,y3)
dev.off()

