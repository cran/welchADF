#** define module to check initial specifications****
# returns: ntot, wobs, wobs1, bobs, r, x
.initial <- function(Cmat, Umat = NULL, y, nx, trimming, per=0.2, bootstrap, effect.size, 
										scaling = TRUE, alpha = 0.05) {

	ntot = 0
	bobs = 0
	wobs= 0
	wobs1 = 0
		
	if(is.null(Umat)){ 
		Umat = diag(ncol(y)) 
	}
	else if(ncol(Umat) > nrow(Umat)){ 
		warning("Possible Error: Number Of Columns Of U Exceeds Number Of Rows")
	}

	if(is.null(Cmat)){
	  Cmat = diag(1)
	}
	
	x = rep(NA, sum(nx))	# preallocate vector
  temp = 1
	for(i in 1:length(nx)){
	 if(nx[i] > 0){
		 x[temp:(temp + nx[i] - 1)] = i # fill nx[i] positions with the integer i
		 temp = temp + nx[i]
	 }
	}
	# example: x = c(1,1,2,2,2,3,3,1,2,2,2,2,3,4,4,2,3)

	# The following code is equivalent to the SAS code x=design(x) 
	new.x = matrix(0, nrow = length(x), ncol = max(x))
	for(i in 1:ncol(new.x)){
		new.x[,i] = (x == i)
	}
	x = new.x 	# x is the design matrix

	# simulate global variables with a tagged list returned by the function	
	return.vars = vector("list", 6)
	names(return.vars) = c("x","ntot","wobs","wobs1","bobs","r")	
	
	if(!is.matrix(y)){	
		return.vars$ntot = length(y) 
		return.vars$wobs = 1
	}
	else{								
		return.vars$ntot = nrow(y) 
		return.vars$wobs = ncol(y)
	}		
	return.vars$wobs1 = return.vars$wobs-1 
	return.vars$bobs = ncol(x)
	if(is.matrix(Cmat)){	
		return.vars$r = Cmat %x% t(Umat) 		# Kronecker product	.
	}else{ 
		return.vars$r = t(Cmat) %x% t(Umat) 		# Kronecker product	.
	}
	if(is.vector(return.vars$r)){
    return.vars$r = t(return.vars$r)
	}
	return.vars$x = x

	return(return.vars)
}

# ________________________________________________________________________________

# ****define module to compute least squares or trimmed means**** 
# returns: 	bhat,bhatw, yt, dfr, muhat
.mnmod <- function(y,x,opt1,per, bobs,wobs,wobs1,ntot,nx){

	return.vars = vector("list", 5)
	names(return.vars) = c("bhat", "bhatw", "yt", "dfr", "muhat")
	
	if(!opt1){  				# use ordinary least squares
		return.vars$bhat = solve(t(x) %*% x, tol = 1E-30) %*% t(x) %*% y # solve the normal equations for least squares
		return.vars$bhatw = return.vars$bhat 
		return.vars$yt = y 
		return.vars$dfr = nx-1 
	}
	else{							# use trimmed means
		return.vars$bhat = matrix(0,bobs,wobs) 
		return.vars$bhatw = return.vars$bhat 
		return.vars$yt = matrix(0,ntot,wobs) 
		return.vars$dfr = rep(0, times = bobs) 
		f = 1 
		m = 0
		l = NULL
		for(j in 1:length(nx)){ 
			samp = nx[j] 
			l = m + samp 
			g = as.integer(per*samp) 
			return.vars$dfr[j] = samp - 2*g - 1 
			for(k in 1:ncol(y)){ 
				temp = y[f:l,k] 
				nv = temp 
				temp = sort(temp)				
				trimy = temp[(g+1):(samp-g)] 
				trimmn = sum(trimy)/(return.vars$dfr[j]+1)
				return.vars$bhat[j,k] = trimmn 
				mint = min(trimy) 
				maxt = max(trimy) 
				nv[nv < mint] = mint
				nv[nv > maxt] = maxt
				return.vars$yt[f:l,k] = nv
				return.vars$bhatw[j,k] = sum(nv)/samp
			}
			m = l 
			f = f + nx[j] 
		} 
	}
	return.vars$muhat = as.vector(t(return.vars$bhat)) # turn bhat into a vector taking elements row by row	
	return(return.vars)
}

# ________________________________________________________________________________

# *** define module to compute sigma matrix**** 
# returns: sigma, stdizer
.sigmod <-function (yt,x,bhatw,dfr,bobs,wobs,wobs1,ntot,nx){

  return.vars = vector("list", 2)
  names(return.vars) = c("sigma", "stdizer")

	return.vars$sigma = matrix(0, nrow = wobs*bobs, ncol = wobs*bobs) 
	return.vars$stdizer = return.vars$sigma 	
	for(i in 1:bobs){
	  
    auxmat = (yt*x[,i] - t(t(x[,i])) %*% bhatw[i,])
		sigb= (t(auxmat) %*% auxmat) / ( (dfr[i]+1) * dfr[i] )
		#sigb[sigb==0] = 1E-2 # to avoid cells having 0

		f= i*wobs - wobs1 
		l= i*wobs 
		return.vars$sigma[f:l,f:l]=sigb 
		return.vars$stdizer[f:l,f:l] = sigb * ( (dfr[i]+1) * dfr[i] )/( nx[i]-1 ) 
	}
	return(return.vars)
}

# ________________________________________________________________________________

#**** define module to compute test statistic**** 
# returns: fstat, dfr1, dfr2
.testmod <- function(sigma,muhat,r,dfr,
                    bobs,wobs,wobs1,ntot,nx){

  return.vars = vector("list", 3)
  names(return.vars) = c("fstat", "dfr1", "dfr2")

  if(is.vector(r)){ r = t(r) }

  if(det(r %*% sigma %*% t(r)) == 0){
    stop("R * sigma * R' is a singular matrix so it cannot be inverted")
  }

	temp.t = t(r %*% muhat) %*% solve(r %*% sigma %*% t(r), tol = 1E-30) %*% (r %*% muhat)
	temp.t = as.vector(temp.t)
	a=0 
	imat = diag(wobs) 
	for(i in 1:bobs){
	 qmat = matrix(0, nrow = bobs*wobs, ncol = bobs*wobs) 
	 f = i*wobs-wobs1 
	 l = i*wobs 
	 qmat[f:l,f:l] = imat 
	 prod = (sigma %*% t(r)) %*% solve((r %*% sigma) %*% t(r), tol = 1E-30) %*%r %*% qmat 
	 a = a + ( sum(diag(prod %*% prod)) + sum(diag(prod))*sum(diag(prod)) )/dfr[i] 
	}

	a = a/2
	
	dfr1 = nrow(r) 
	dfr2 = dfr1*(dfr1+2)/(3*a) 
	cval = dfr1 + 2*a-6*a/(dfr1+2) 
	
	return.vars$fstat = temp.t/cval 
	return.vars$dfr1 = dfr1
	return.vars$dfr2 = dfr2
	
	return(return.vars)
}

# ________________________________________________________________________________
  