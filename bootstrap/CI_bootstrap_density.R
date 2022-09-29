library(Rcpp)
A<-read.table("data_Wuhan.txt") 
sourceCpp("bootstrap_density.cpp")
output <- ComputeIntervals_density(A)
B <- output$CI_density
   x1<-B[,1]
   y1<-B[,2]
   y2<-B[,3]
   y3<-B[,4]
   
    a<-3.035140901
	b<-0.002619475
	f <- function(x) {a*b*exp(-b*x^a)*x^(-1+a)/(1-exp(-b*20^a))}	
	x0 <-seq(0,15,by=0.01)
	y0<-f(x0)


pdf("CI_density.pdf")
plot(c(-100,-100),xlim=c(0,15), ylim=c(0,max(y3)), main= "",ylab="",xlab="",bty="n",las=1)
lines(x1,y1,lwd=2,col="blue")
lines(x0,y0,lwd=2,lty=2,col="red")
segments(x1,y2,x1,y3)
dev.off()