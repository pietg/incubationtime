	f1=function(u,v,alpha)
 	{
  		a=0 
 		if (u>0)
 			a <- log(0.5*erf((log(v)-alpha[1])/(alpha[2]*sqrt(2)))-0.5*erf((log(u)-alpha[1])/(alpha[2]*sqrt(2))))
 		else
 			a <- log(0.5+0.5*erf((log(v)-alpha[1])/(alpha[2]*sqrt(2)))) 	
  		return(a)
 	}
 
    
	eval_f1<-function(x,data)
	{
		calc=0;
		for (i in 1:n)
			calc = calc+f1(data[i,1],data[i,2],x)
		return(-calc)
	}

	# Lower and upper bounds
	lb1 <- c(0,0)
	ub1 <- c(5,5)

	eval_g_ineq1<-function(x,data)
	{
		return(c(-x[1],-x[2]))
	}
	