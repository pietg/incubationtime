	p0 <- function(s,a,b){0.5+0.5*erf((log(s)-a)/(sqrt(2)*b))}
	p1 <- function(s,a,b){exp(-(log(s)-a)^2/(2*b^2))/sqrt(pi)}
	p2 <- function(s,a,b){(log(s)-a)/(b*sqrt(2))}
	
	# numerator of derivatives w.r.t. a and b of p0
	p2a <- Deriv(p2,"a")
	p2b <- Deriv(p2,"b")
	
	q1a <- function(s,a,b){p2a(s,a,b)*p1(s,a,b)}
	q1b <- function(s,a,b){p2b(s,a,b)*p1(s,a,b)}
	
	q2a <- Deriv(q1a,"a")
	q2b <- Deriv(q1b,"b")
	q2ab <- Deriv(q1a,"b")
	
	# 2nd derivatives w.r.t. a and b of p0
	g11 <- function(s,a,b){-q2a(s,a,b)/p0(s,a,b)+q1a(s,a,b)^2/p0(s,a,b)^2}
	g22 <- function(s,a,b){-q2b(s,a,b)/p0(s,a,b)+q1b(s,a,b)^2/p0(s,a,b)^2}
	g12 <- function(s,a,b){-q2ab(s,a,b)/p0(s,a,b)+q1a(s,a,b)*q1b(s,a,b)/p0(s,a,b)^2}
	
	r0 <- function(e,s,a,b){0.5*(erf((log(s)-a)/(sqrt(2)*b))-erf((log(s-e)-a)/(sqrt(2)*b)))}
	r1 <- function(e,s,a,b){exp(-(log(s-e)-a)^2/(2*b^2))/sqrt(pi)}
	r2 <- function(e,s,a,b){(log(s-e)-a)/(sqrt(2)*b)}
	
	r2a <- Deriv(r2,"a")
	r2b <- Deriv(r2,"b")
	
	s1a <- function(e,s,a,b){p2a(s,a,b)*p1(s,a,b)-r2a(e,s,a,b)*r1(e,s,a,b)}
	s1b <- function(e,s,a,b){p2b(s,a,b)*p1(s,a,b)-r2b(e,s,a,b)*r1(e,s,a,b)}
	
	s2a <- Deriv(s1a,"a")
	s2b <- Deriv(s1b,"b")
	s2ab <- Deriv(s1a,"b")

	# 2nd derivatives w.r.t. a and b of r0
	h11 <- function(e,s,a,b){-s2a(e,s,a,b)/r0(e,s,a,b)+s1a(e,s,a,b)^2/r0(e,s,a,b)^2}
	h22 <- function(e,s,a,b){-s2b(e,s,a,b)/r0(e,s,a,b)+s1b(e,s,a,b)^2/r0(e,s,a,b)^2}
	h12 <- function(e,s,a,b){-s2ab(e,s,a,b)/r0(e,s,a,b)+s1a(e,s,a,b)*s1b(e,s,a,b)/r0(e,s,a,b)^2}

