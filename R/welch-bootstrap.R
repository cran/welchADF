#**** define modules to perform bootstrap**** 
#*** define module to generate bootstrap data**** 
# returns: yb1
.bootdat<-function(y,bobs,ntot,nx){
	
	return.vars = vector("list", 1)
	names(return.vars) = c("yb1")
	return.vars$yb1 = matrix(0, nrow = sum(nx), ncol = ncol(y))

	f=1 
	m=0	
	for(j in 1:bobs){ 	
		l = m + nx[j] 
		temp = y[f:l,] 		
		rows = nrow(temp)
		if(is.vector(temp)){ rows = length(temp) }
		## change the ordering of the rows randomly
		rval = runif(rows)
		rval = ceiling(rval * rows)
		if(is.vector(temp)){
			bval = temp[rval]
		}else{
			bval = temp[rval,]
		}
		
		return.vars$yb1[f:l,] = bval
		m=l 
		f=f+nx[j] 
	}
	
	return(return.vars)
}

## ________________________________________________________

#**** define module to centre the bootstrap data**** 
# returns: yb
.bootcen <- function(yb1,bhat,bobs,ntot,nx){
  return.vars = vector("list", 1)
  names(return.vars) = c("yb")
	m=0 
	f=1 
	return.vars$yb=yb1 

	for(i in 1:bobs){ 
		l=m+nx[i] 
		mval=bhat[i,]		
		ncolumns = 1
		if(!is.vector(mval)){ ncolumns = ncol(mval) }
		else{ 								ncolumns = length(mval) }
		# replicate vector val in all rows
		identical.rows = matrix(mval, nrow = l-f+1, ncol = ncolumns, byrow = TRUE)
		return.vars$yb[f:l,] = return.vars$yb[f:l,] - identical.rows # the same vector is subtracted from every row
		#do q=f to l by 1 
		#	yb[q,]=yb1[q,]-mval 
		#end 
		m=l 
		f=f+nx[i] 
	}
	
	return(return.vars)
}

## ________________________________________________________

#**** define module to compute bootstrap test statistic****
# returns: fstat, dfr1, dfr2
.bootstat <- function(yb,x,opt1,per,r, bobs,wobs,wobs1,ntot,nx){

	result1 = .mnmod(yb, x, opt1, per, bobs, wobs, wobs1, ntot, nx) # returns: bhatb,bhatbw,muhatb,yt,dfr
	result2 = .sigmod(result1$yt,x,result1$bhatw,result1$dfr,bobs,wobs,wobs1,ntot,nx) # returns: sigma,stdizer
	result3 = .testmod(result2$sigma,result1$muhat,r,result1$dfr, bobs,wobs, wobs1, ntot, nx) # returns: fstat,dfr1,dfr2
	
	return(result3)
}

## ________________________________________________________

#**** define module to compute bootstrap effect size**** 
#returns: multp, effsz, stdizer
.bootes <- function(yb,x,standardizer,scale,opt1,per,r,
										bobs,wobs,wobs1,ntot,nx){
										
	result1 = .mnmod(yb, x, opt1, per, bobs, wobs, wobs1, ntot, nx) # returns: bhat,bhatw,muhat,yt,dfr
	result2 = .sigmod(result1$yt,x,result1$bhatw,result1$dfr,bobs,wobs,wobs1,ntot,nx) # returns: sigma,stdizer
  result3 = .wjeffsz(standardizer, scale,opt1,per,r,result1$muhat,result2$stdizer,bobs,wobs,wobs1,ntot,nx)  #returns: multp, effsz
	
	return(result3)
}

## ________________________________________________________

#**** compute measure of effect size and bootstrap confidence interval**** 
#returns: multp, effsz, stdizer
.wjeffsz <- function(standardizer, scale,opt1,per,r,muhat,stdizer, bobs,wobs,wobs1,ntot,nx){ 
	return.vars = vector("list", 2)
	names(return.vars) = c("multp", "effsz")
	
	if(!opt1 || (opt1 && !scale)){ return.vars$multp=1	}  
	else if(opt1 && scale){
		cut=qnorm(per) # quantile function of the standard normal distribution N(0,1)
		# Integrate the probability density function of the standard normal along the intervals [cut, 0] and [0, -1*cut]
		myf = function(z){		(z*z)*(1/(sqrt(2*pi))*(exp(-(1/2)*(z*z))))		}
		# b = integrate(dnorm, lower = cut, upper = 0)$value + integrate(dnorm, lower = 0, upper = -cut)$value		
		b = integrate(myf, lower = cut, upper = 0)$value + integrate(myf, lower = 0, upper = -cut)$value
		winvar = b + per*(2*(cut*cut))
		return.vars$multp = sqrt(winvar)
	}

	num = r %*% muhat 
	stdz = NULL

	if(!standardizer){ stdz = 1 }
	else{
	  # Standardizer is the Square Root of the average of variances
		r2=r*r 	# element wise multiplication (squares every element of r)
		#rvec = shape(r2,bobs#wobs)` 
		rvec = t(as.vector(t(r2))) # column vector
		stdz1 = sqrt(diag(diag(stdizer), nrow(stdizer), ncol(stdizer)))
		stdz = (rvec %*% stdz1 %*% t(rvec))/sum(rvec) # average variance
	}
	return.vars$stdizer = stdz
	return.vars$effsz = return.vars$multp * (num/stdz) 

	return(return.vars)
} 