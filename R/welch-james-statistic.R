wjglm <- function(Cmat, Umat, y, nx, trimming, per, bootstrap, numsim_b, effect.size, numsim_es,
										standardizer, scaling, alpha, seed){

	opt1 = trimming
	opt2 = bootstrap
	opt3 = effect.size

	if(seed != 0){ set.seed(seed) }
	
  if(is.vector(y)){
    y = t(t(y))
  }

	vars.initial = .initial(Cmat,Umat,y,nx,opt1,per,opt2,opt3,scaling,alpha)
	
	ntot = vars.initial$ntot 
	wobs = vars.initial$wobs 
	wobs1 = vars.initial$wobs1 
	bobs = vars.initial$bobs 	
	r = vars.initial$r 		
	x = vars.initial$x
		
	vars.mnmod = .mnmod(y,x,opt1,per, bobs,wobs,wobs1,ntot,nx)
	
	bhat = vars.mnmod$bhat
	bhatw = vars.mnmod$bhatw
	muhat = vars.mnmod$muhat
	yt = vars.mnmod$yt
	dfr = vars.mnmod$dfr

	vars.sigmod = .sigmod(yt,x,bhatw,dfr,bobs,wobs,wobs1,ntot,nx)
	sigma = vars.sigmod$sigma
	stdizer = vars.sigmod$stdizer

	vars.testmod = .testmod(sigma,muhat,r,dfr,bobs,wobs,wobs1,ntot,nx)
	fstat = vars.testmod$fstat
	dfr1 = vars.testmod$dfr1
	dfr2 = vars.testmod$dfr2
	
	multp = NULL
	effsz = NULL
	fmat = NULL
	esmat = NULL
	
	if(opt3){  # compute effect size
		vars.wjeffsz = .wjeffsz (standardizer,scaling,opt1,per,r,muhat,stdizer, bobs,wobs,wobs1,ntot,nx)
		multp = vars.wjeffsz$multp
		effsz = vars.wjeffsz$effsz
		stdizer = vars.wjeffsz$stdizer
	}
	
	if(opt2){ 	# use bootstrap
		fmat = replicate(numsim_b, expr=		
			{	
				yb1 = (.bootdat(y,bobs,ntot,nx))$yb1
				yb = (.bootcen(yb1,bhat,bobs,ntot,nx))$yb
				res = .bootstat(yb,x,opt1,per,r, bobs,wobs,wobs1,ntot,nx)	# value obtained by this block of function calls
				res$fstat
			})
		fmat = sort(fmat)		
	}
	
	if(opt3){   # compute effect size using bootstrap
		esmat = replicate(numsim_es, expr = 
			{	
				yb1 = (.bootdat(y,bobs,ntot,nx))$yb1
				res = .bootes(yb1,x,standardizer,scaling,opt1,per,r, bobs,wobs,wobs1,ntot,nx) # value obtained by this block of function calls
				res$effsz
			})
		esmat = sort(esmat)
	}
	#*** calculate significance level for welch-james statistic****
	result.list = list()

	if(!opt2){ # do not use bootstrap estimators
		# pf is the (cumulative) distribution function of an F distribution
		pval = 1 - pf(fstat,dfr1,dfr2)
		result.list$pval = pval
	} 
	else{		# use bootstrap estimators
		avec = (fmat >= fstat)
		pval = sum(avec)/numsim_b
		result.list$pval = pval
	}
		
	lcl = NULL
	ucl = NULL
	
	if(opt3){		# compute effect size
		ind1 = as.integer(numsim_es * (alpha/2)) + 1
		ind2 = numsim_es - as.integer((numsim_es * (alpha/2)))
		lcl = esmat[ind1]
		ucl = esmat[ind2]
	}
	
	result.list$Cmatrix = Cmat
	result.list$Umatrix = Umat
	result.list$welch.T = fstat
	result.list$numeratorDF = dfr1
	result.list$denominatorDF = dfr2
	result.list$contrast.matrix = r
	result.list$mean.vector = muhat
	result.list$sigma.matrix = sigma
	result.list$trimming = as.logical(opt1)
	result.list$bootstrap = as.logical(opt2)
	result.list$compute.effsz = as.logical(opt3)
	result.list$scaling = as.logical(scaling)
	result.list$standardize.effsz = as.logical(standardizer)
	
	#	****store results****
	if(!opt1){ 
	}
	else{ 
		result.list$trimming.per = per
	}
	if(!opt2 && !opt3){ 
	}
	if(opt2 && !opt3){ 
		result.list$boot.samples = numsim_b
		result.list$seed = seed		
	}
	
	if(!opt3){ 	# don't compute effect size
	}
	else{ 
		cilev = (1 - alpha)*100
		result.list$effsz = effsz
		result.list$stdizer
		
		if(!scaling){ 
		}
		else{ 
			result.list$scaling.factor = multp
		}

		result.list$effsz.samples = numsim_es
		result.list$seed = seed
		result.list$CI.level = cilev
		result.list$CI = c(lcl, ucl)
		
	}
	
	class(result.list) = "WelchTestObj"
	return(result.list)
}

# ______________________________________________________
