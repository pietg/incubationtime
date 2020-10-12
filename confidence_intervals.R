library(Rcpp)
sourceCpp("confidence_intervals.cpp")
output <- CI()

C<-output$MLE
D<-output$SMLE
E<-output$dens
CI<-output$CI

a<-3.035140901
b<-0.002619475
f <- function(x) {a*b*exp(-b*x^a)*x^(-1+a)}

x1<-E[,1]
y1<-E[,2]
u1<-CI[,1]
u2<-CI[,2]
u3<-CI[,3]
x0 <-seq(0,max(x1),by=0.01)
y0<-f(x0)

plot(c(-10000,-10000),xlim=c(0,14), ylim=c(0.0,max(u3)), main= "", ylab="",xlab="",bty="n",las=1)

lines(x0,y0,lwd=2,lty=2,col="red")
lines(x1,y1,lty=1,lwd=2,col="blue")
segments(u1,u2,u1,u3)







