   B<-read.table("MLE.txt")
   x<-B[,1]
   y<-B[,2]
   a<-3.035140901
   b<-0.002619475
   #a<-3.0360516078
   #b <- 0.0026148586
   f <- function(x) {1-exp(-b*x^a)}/{1-exp(-b*20^a)}
   x0 <-seq(0,max(x),by=0.01)
   y0<-f(x0)
   plot(c(-1000,-1000),xlim=c(min(x),max(x)), ylim=c(min(y,y0),max(y,y0)), main= "",ylab="",xlab="",bty="n",las=1)
   lines(x,y,lwd=2,col="blue",type='s')
   lines(x0,y0,lwd=2,lty=2,col="red")  
	
	
	
	

    