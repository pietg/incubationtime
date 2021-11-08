library(Rcpp)
A<-read.table("data_Wuhan.txt") 
sourceCpp("bootstrap_SMLE.cpp")
output <- ComputeIntervals_df(A)
B <- output$CI_df
   x1<-B[,1]
   y1<-B[,2]
   y2<-B[,3]
   y3<-B[,4]

pdf("CI_SMLE.pdf")
plot(c(-100,-100),xlim=c(0,12), ylim=c(0,max(y3)), main= "",ylab="",xlab="",bty="n",las=1)
lines(x1,y1,lwd=2,col="red")
segments(x1,y2,x1,y3)
dev.off()