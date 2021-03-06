library(lbfgs)
   A<-read.table("inputdata_Wuhan.txt")
   data1<-A[,1]
   data2<-A[,2]
   n<-88
   alpha=vector("numeric",2)
   
criterion=function(alpha)
{
	calc=0;
	for (i in 1:n)
	{
		calc = calc+log(exp(-alpha[2]*data1[i]^alpha[1])-exp(-alpha[2]*data2[i]^alpha[1]));
	}
	return(-calc);
}

grad=function(alpha)
{
	calc1=0;
	calc2=0;
	for (i in 1:n)
	{
		denom = exp(-alpha[2]*data1[i]^alpha[1])-exp(-alpha[2]*data2[i]^alpha[1]);
		if (denom>0)
		{
			if (data1[i]==0)
			{
			  calc1 = calc1+alpha[2]*exp(-alpha[2]*data2[i]^alpha[1])*(data2[i]^alpha[1])*log(data2[i])/denom;
			}
			else
			{
				calc1 = calc1+alpha[2]*exp(-alpha[2]*data2[i]^alpha[1])*(data2[i]^alpha[1])*log(data2[i])/denom-alpha[2]*exp(-alpha[2]*data1[i]^alpha[1])*(data1[i]^alpha[1])*log(data1[i])/denom;
			}
			calc2 = calc2+exp(-alpha[2]*data2[i]^alpha[1])*(data2[i]^alpha[1])/denom-exp(-alpha[2]*data1[i]^alpha[1])*(data1[i]^alpha[1])/denom;
		}
	}
	return(c(-calc1,-calc2));
}

output <- lbfgs(criterion, grad, c(2,0.002))
output$par