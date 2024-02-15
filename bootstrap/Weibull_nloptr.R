	f0=function(u,v,alpha)
	{
  		a=0 
 		if (u>0)
 			a <- log(exp(-(u/alpha[2])^alpha[1])-exp(-(v/alpha[2])^alpha[1]))
 		else
		 	a <- log(1-exp(-(v/alpha[2])^alpha[1])) 	
  	return(a)
	}
 
    
	eval_f0<-function(x,data)
	{
		calc=0;
		for (i in 1:n)
			calc = calc+f0(data[i,1],data[i,2],x)
		return(-calc)
	}

	# Lower and upper bounds
	lb0 <- c(1,1)
	ub0 <- c(Inf,Inf)

	eval_g_ineq0<-function(x,data)
	{
		return(c(-x[1],-x[2]))
	}
	