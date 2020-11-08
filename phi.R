### R-code for phi
### Piet Groeneboom


   B<-read.table("phi2.txt")
   x1<-B[,1]
   y1<-B[,2]
   y2<-B[,3]
   a1<-min(x1)
   b1<-max(x1)
   a2<-min(y1,y2)
   b2<-max(y1,y2)
   plot(c(-10000,-10000),xlim=c(0,20), ylim=c(a2,b2), main= "", ylab="",xlab="",bty="n",las=1)
   lines(x1, y1,lwd=2,col="red")
   lines(x1, y2,lwd=2,col="blue",lty=2)


  