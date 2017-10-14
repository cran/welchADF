.start.bootcom <- function(Cmat, Umat, y, nx, trimming, per, bootstrap, numsim_b, effect.size, numsim_es, 
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

  temp = .testmod(sigma, muhat, r, dfr, bobs,wobs,wobs1,ntot,nx)
  test.statistic = unlist(temp)

  fmat = replicate(numsim_b, 
  {
    yb1 = (.bootdat(y,bobs, ntot, nx))$yb1
    yb = (.bootcen(yb1, bhat, bobs, ntot, nx))$yb
    (.bootstat(yb, x, opt1, per, r, bobs,wobs,wobs1,ntot,nx))$fstat
  })

  # -----------------------------
  #   compute effect size in case it is needed
  # -----------------------------
  
  multp = NULL
  effsz = NULL
  esmat = NULL
  lcl = NULL
  ucl = NULL
  
  if(opt3){
    vars.wjeffsz = .wjeffsz (standardizer,scaling,opt1,per,r,muhat,stdizer, bobs,wobs,wobs1,ntot,nx)
    multp = vars.wjeffsz$multp
    effsz = vars.wjeffsz$effsz
    
    esmat = replicate(numsim_es, expr = 
    {	
      yb1 = (.bootdat(y,bobs,ntot,nx))$yb1
      res = .bootes(yb1,x,standardizer,scaling,opt1,per,r, bobs,wobs,wobs1,ntot,nx) # value obtained by this block of function calls
      res$effsz
    })
    esmat = sort(esmat)
    
    ind1 = as.integer(numsim_es * (alpha/2)) + 1
    ind2 = numsim_es - as.integer((numsim_es * (alpha/2)))
    lcl = esmat[ind1]
    ucl = esmat[ind2]
    
  }
  
  result.list = list()
  result.list$Cmatrix = Cmat
  result.list$Umatrix = Umat
  result.list$welch.T = test.statistic["fstat"]
  result.list$numeratorDF = test.statistic["dfr1"]
  result.list$denominatorDF = test.statistic["dfr2"]
  result.list$contrast.matrix = r
  result.list$mean.vector = muhat
  result.list$sigma.matrix = sigma
  result.list$trimming = as.logical(opt1)
  result.list$bootstrap = as.logical(opt2)
  result.list$compute.effsz = as.logical(opt3)
  result.list$scaling = as.logical(scaling)
  result.list$standardize.effsz = as.logical(standardizer)
  result.list$fmat = fmat # bootstrap samples

  if(!opt1){ 
  }
  else{ 
    result.list$trimming.per = per
  }
  if(!opt2 && !opt3){ 
  }
  else if(opt2 && !opt3){ 
    result.list$boot.samples = numsim_b  
    result.list$seed = seed
  }
  else if(opt3){
    cilev = (1 - alpha)*100
    result.list$effsz = effsz
    if(!scaling){ 
    }
    else{ 
      result.list$scaling.factor = multp
    }

    # Standardizer Is always Square Root of Average Variance
    result.list$effsz.samples = numsim_es
    result.list$seed = seed
    result.list$CI.level = cilev
    result.list$CI = c(lcl, ucl)
  }
  
  return(result.list)
}

# ______________________________________________________________________________

.compute.bootcom.muhat.sigma.r <- function(Cmat, Umat, y, nx, trimming, per, bootstrap, numsim_b, effect.size, numsim_es, 
                                          scaling, alpha, seed){
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
  
  return(list(sigma = sigma, muhat = muhat, r = r))
}

# ______________________________________________________________________________

.end.bootcom <- function(fmat, alpha, numsim_b){

  fmax = fmat

  if(is.matrix(fmat)){
    fmat = t(fmat)
    fmax = apply(fmat, 1, max)
  }
  fmax = sort(fmax)
  qcrit = round((1-alpha) * numsim_b)
  critv = fmax[qcrit]
  return(critv)
}

# ______________________________________________________________________________
