
.compute.pairwise <- function(data, response, between.s, within.s, subject, effect, 
				correction, trimming, per, bootstrap, numsim_b, effect.size, numsim_es, 
				standardizer, scaling, alpha, seed, y, nx, levelslist.between.s, levelslist.within.s){

  	frameNames = names(data)	

		ncols = ncol(data)
		factorColumns = (1:ncols)[!(frameNames %in% response)]
		columns.between.s = (1:ncols)[(frameNames %in% between.s)]
		columns.within.s = (1:ncols)[(frameNames %in% within.s)]
		nresponses = length(response)
		
		opt1 = trimming
		opt2 = bootstrap
		opt3 = effect.size

		## ----------------------------------------	
		##	Generate all C and U matrices needed
		## ----------------------------------------	
		## All pairwise comparisons between the levels of the "effect" argument
		# Find out if the effect is a between-subject or a within-subject effect
		effect.between.s = NULL
		effect.within.s = NULL
		all.C.matrices = NULL
		all.U.matrices = NULL
		names.all.C = NULL
		names.all.U = NULL
		
		Cmatrix_bootstrap = NULL
		Umatrix_bootstrap = NULL
		bootstrap_sigma = NULL
		bootstrap_muhat = NULL
		bootstrap_r = NULL
		
		allnamesboot.between.s = NULL
		allnamesboot.within.s = NULL
		
		if(is.null(within.s) && nresponses > 1){ 	# multivariate setting
			all.U.matrices = list(diag(nresponses))
			Umatrix_bootstrap = all.U.matrices[[1]]
		}
		
		if(!is.null(between.s)){
			## -------------------------------------------
			## The model contains between-subjects factors
			## -------------------------------------------
			nlevelsvec = sapply(levelslist.between.s, FUN = length)
			mylist = lapply(nlevelsvec, FUN = rep.int, x = 1)
			mylist = lapply(mylist, FUN = t) # list of vectors, each of the same length as the number of levels of that variable			
			effect.between.s = effect[effect %in% frameNames[columns.between.s]]
			
			if(length(effect.between.s)>0){
			
				## -------------------------------------------------------------------------------
				## One or more between-subjects effects are involved in the requested pairwise comparisons
				## -------------------------------------------------------------------------------
				
				# Put the names in effect.between.s in the same order as the columns of the dataset
				positions = match(frameNames, effect.between.s)
				positions = positions[!is.na(positions)]
				effect.between.s = effect.between.s[positions]
				
				mynames1 = levelslist.between.s[[effect.between.s[1]]]
				
				# Generate contrasts vectors for all pairwise combinations
				comb1 = .generateAllPairwiseCombinations(nlevelsvec[effect.between.s[1]], mynames1)
				
				## Gather the pair of level names involved in each pairwise comparison in a 
				## list of string vectors
				comb1logical = lapply(comb1, FUN = as.logical)
				allnames1 = lapply(comb1logical, FUN = function(v){
					l = as.list(c(mynames1[v], sep = ":"))
					do.call(paste, l)
				})

				# This might be overwritten just inside the following "if"				
				allnamesboot.between.s = as.character(allnames1)

				if(length(effect.between.s) == 2){
					
					mynames2 = levelslist.between.s[[effect.between.s[2]]]
					
					# Generate contrasts vectors for all pairwise combinations
					comb2 = .generateAllPairwiseCombinations(nlevelsvec[effect.between.s[2]], levelslist.between.s[[effect.between.s[2]]])
					
					## Gather the pair of level names involved in each pairwise comparison in a 
					## list of string vectors					
					comb2logical  = lapply(comb2, FUN = as.logical)
					allnames2 = lapply(comb2logical, FUN = function(v){
						l = as.list(c(mynames2[v], sep = ":"))
						do.call(paste, l)
					})
					
					allnamesboot = sapply(allnames1, function(n1){
						paste(n1, allnames2, sep = " x ")
					})
					allnamesboot.between.s = as.character(allnamesboot)
					
					if(bootstrap){	# run bootcom: all pairwise contrast vectors of a factor must be collected in a matrix for that factor
						## For between-subjects factors, every row is a contrast vector
						matcomb1 = matrix(unlist(comb1), ncol = length(comb1[[1]]), byrow = TRUE) 
						matcomb2 = matrix(unlist(comb2), ncol = length(comb2[[1]]), byrow = TRUE)
						mytemplist = mylist
						mytemplist[[effect.between.s[1]]] = matcomb1
						mytemplist[[effect.between.s[2]]] = matcomb2
						Cmatrix_bootstrap = .kronecker.all(mytemplist)
						#names.all.C[[1]] = do.call(paste, c(as.list(effect.between.s), sep = " x "))
					}
					
					# For each pairwise contrasts vector, substitute it in the expression 
					# of the Kronecker product and then perform the product
					j = 1
					for(vec1 in comb1){
						mylist[[effect.between.s[1]]] = vec1
						names1 = as.list(colnames(vec1)[vec1!=0]) 
						names1$sep = ":"
						names1 = do.call(paste, args = names1)
						for(vec2 in comb2){
							mylist[[effect.between.s[2]]] = vec2
							all.C.matrices[[j]] = .kronecker.all(mylist)
							
							names2 = as.list(colnames(vec2)[vec2!=0])	
							names2$sep = ":"
							names2 = do.call(paste, args = names2)
							involved.levels = paste(names1, names2, sep = " x ")							
							names.all.C = c(names.all.C, involved.levels)
							j = j+1
						}
					}
				}
				else{
					if(bootstrap){	# run bootcom: all pairwise contrast vectors of a factor must be collected in a matrix for that factor

						## For between-subjects factors, every row is a contrast vector
						matcomb1 = matrix(unlist(comb1), ncol = length(comb1[[1]]), byrow = TRUE) 
						mytemplist = mylist
						mytemplist[[effect.between.s[1]]] = matcomb1
						Cmatrix_bootstrap = .kronecker.all(mytemplist)
						#names.all.C[[1]] = effect.between.s
					}
				  
					j = 1
					for(vec1 in comb1){
						names1 = as.list(colnames(vec1)[vec1!=0])
						names1$sep = ":"
						names1 = do.call(paste, args = names1)
						names.all.C = c(names.all.C, names1)
						
						mylist[[effect.between.s[1]]] = vec1
						all.C.matrices[[j]] = .kronecker.all(mylist)
						j = j+1
					}
				}
			}
			else{  
			  # No within-subject factors among the effects involved in the pairwise comparisons
				all.C.matrices[[1]] = .kronecker.all(mylist)
				names.all.C = rep(NA, length(all.C.matrices))
				Cmatrix_bootstrap = all.C.matrices[[1]]
			}
		}
		
		if(!is.null(within.s)){
			## -------------------------------------------
			## The model contains within-subjects factors
			## -------------------------------------------
			nlevelsvec = sapply(levelslist.within.s, FUN = length)
			mylist = lapply(nlevelsvec, FUN = rep.int, x = 1)			
			mylist = lapply(mylist, FUN = t) # list of vectors, each of the same length as the number of levels of that variable
			if(nresponses > 1){	# multivariate setting: append a diagonal matrix as the last factor of the Kronecker product
				mylist[[length(mylist)+1]] = diag(nresponses)
			}
			
			effect.within.s = effect[effect %in% frameNames[columns.within.s]]  
			
			if(length(effect.within.s)>0){
				## -------------------------------------------------------------------------------
				## One or more within-subjects effects are involved in the requested pairwise comparisons
				## -------------------------------------------------------------------------------
				# Put the names in effect.within.s in the same order as the columns of the dataset
				positions = match(frameNames, effect.within.s)
				positions = positions[!is.na(positions)]
				effect.within.s = effect.within.s[positions]
				
				mynames1 = levelslist.within.s[[effect.within.s[1]]]		# level names of this within-subjects factor		
				
				comb1 = .generateAllPairwiseCombinations(nlevelsvec[effect.within.s[1]], mynames1)				
				
				comb1logical = lapply(comb1, FUN = as.logical)
				allnames1 = lapply(comb1logical, FUN = function(v){
					l = as.list(c(mynames1[v], sep = ":"))
					do.call(paste, l)
				})
				
				# This might be overwritten just inside the following "if"
				allnamesboot.within.s = allnames1
				
				if(length(effect.within.s) == 2){
					
					mynames2 = levelslist.between.s[[effect.within.s[2]]]
					
					comb2 = .generateAllPairwiseCombinations(nlevelsvec[effect.within.s[2]], mynames2)
					comb2logical  = lapply(comb2, FUN = as.logical)
					## Gather the pair of level names involved in each pairwise comparison in a 
					## list of string vectors
					allnames2 = lapply(comb2logical, FUN = function(v){
						l = as.list(c(mynames2[v], sep = ":"))
						do.call(paste, l)
					})
					
					allnamesboot = sapply(allnames1, function(n1){
						paste(n1, allnames2, sep = " x ")
					})
					allnamesboot.within.s = as.character(allnamesboot)
					
					if(bootstrap){ 	# run bootcom: all pairwise contrast vectors of a factor must be collected by columns in a matrix for that factor
						## For within-subjects factors, every column is a contrast vector
						matcomb1 = matrix(unlist(comb1), nrow = length(comb1[[1]]), byrow = FALSE) # collected by column
						matcomb2 = matrix(unlist(comb2), nrow = length(comb2[[1]]), byrow = FALSE) # collected by column
						mytemplist =  mylist
						mytemplist[[effect.within.s[1]]] = matcomb1
						mytemplist[[effect.within.s[2]]] = matcomb2
						Umatrix_bootstrap = t(.kronecker.all(mytemplist))
						#names.all.U[[1]] = do.call(paste, c(as.list(effect.within.s), sep = " x "))
					}
					
					# For each pairwise contrasts vector, substitute it in the expression 
					# of the Kronecker product and then perform the product
					j = 1					
					for(vec1 in comb1){					
						mylist[[effect.within.s[1]]] = t(t(vec1))
						names1 = as.list(colnames(vec1)[vec1!=0]) 
						names1$sep = ":"
						names1 = do.call(paste, args = names1)
						for(vec2 in comb2){
							mylist[[effect.within.s[2]]] = t(t(vec2))
							all.U.matrices[[j]] = t(.kronecker.all(mylist))
							
							names2 = as.list(colnames(vec2)[vec2!=0])	
							names2$sep = ":"
							names2 = do.call(paste, args = names2)
							involved.levels = paste(names1, names2, sep = " x ")
							names.all.U[[j]] = c(names.all.U, involved.levels)							
							j = j+1
						}
					}
				}
				else{
					if(bootstrap){	# run bootcom: all pairwise contrast vectors of a factor must be collected in a matrix for that factor
						## For within-subjects factors, every column is a contrast vector
						matcomb1 = matrix(unlist(comb1), nrow = length(comb1[[1]]), byrow = FALSE)
						mytemplist = mylist
						mytemplist[[effect.within.s[1]]] = t(matcomb1)
						Umatrix_bootstrap = t(.kronecker.all(mytemplist))					
						#names.all.U[[1]] = effect.within.s
					}
					j = 1
					for(vec1 in comb1){
						names1 = as.list(colnames(vec1)[vec1!=0]) 
						names1$sep = ":"
						names1 = do.call(paste, args = names1)
						names.all.U = c(names.all.U, names1)
						
						mylist[[effect.within.s[1]]] = t(t(vec1))
						all.U.matrices[[j]] = t(.kronecker.all(mylist))
						j = j+1
					}
				}
			}
			else{	 # No within-subjects factors involved in the effect on which pairwise comparisons are being tested
				all.U.matrices[[1]] = t(.kronecker.all(mylist))
				names.all.U = rep(NA, length(all.U.matrices))
				Umatrix_bootstrap = all.U.matrices[[1]]
			}
		}
		## ----------------------------------------	
		##  Run the test with every combination of an
		##  element from all.C.matrices and another from all.U.matrices
		## ----------------------------------------	
		all.results = list()
		if(bootstrap){
		  
		  # Compute the augmented R matrix, muhat and sigma (just for the summary)
		  temp = .compute.bootcom.muhat.sigma.r(Cmatrix_bootstrap, Umatrix_bootstrap, y, nx, trimming, per=per, bootstrap = bootstrap, 
		                                numsim_b = numsim_b, effect.size = effect.size, numsim_es = numsim_es, 
		                                           scaling = scaling, alpha = alpha, seed = seed)
		  bootstrap_r = temp$r
		  bootstrap_muhat = temp$muhat
		  bootstrap_sigma = temp$sigma
		}
		
		i=1
		if(!is.null(between.s)){
			names(all.C.matrices) = names.all.C
			if(!is.null(within.s)){
				## ----------------------------------------	
				## Both between and within-subject effects
				## ----------------------------------------					
				names(all.U.matrices) = names.all.U
				all.results = vector("list", length(all.C.matrices)*length(all.U.matrices))
				for(matC in all.C.matrices){
					for(matU in all.U.matrices){
					  if(bootstrap){
					    # First part of bootcom (compute the test statistic for this matrix combination using bootstrap)
					    all.results[[i]] = .start.bootcom(Cmat = matC, Umat= matU, y = y, nx = nx, trimming = trimming,
                per=per, bootstrap = bootstrap, numsim_b = numsim_b, effect.size = effect.size,
                standardizer = standardizer, numsim_es = numsim_es, scaling = scaling, alpha = alpha, seed = seed)
					  }
					  else{
							all.results[[i]] = wjglm(Cmat = matC, Umat = matU, y = y, nx = nx,
									trimming = trimming, bootstrap = bootstrap, seed = seed, per = per,
									numsim_b = numsim_b, numsim_es = numsim_es,
									standardizer = standardizer, scaling = scaling, alpha = alpha, effect.size = effect.size)
					  }
						i = i+1
					}
				}
				
				mixed.results.names = sapply(names(all.C.matrices), FUN <- function(matname){
					temp = NULL
					if(is.na(matname)){		temp = names(all.U.matrices)					}
					else{									temp = paste(matname, names(all.U.matrices), sep = " x ")	}
					temp[is.na(names(all.U.matrices))] = matname
					temp
				})
				
				mixed.results.names[is.na(mixed.results.names)] = ""
				for(i in 1:length(all.results)){
					all.results[[i]]$effect.name = mixed.results.names[i]
				}
				names(all.results) = mixed.results.names
			}
			else{
				## ----------------------------------------	
				## Only between-subject effects
				## ----------------------------------------	
			  if(bootstrap){
			    # First part of bootcom (compute the test statistic for this matrix combination using bootstrap)
			    all.results = lapply(all.C.matrices, FUN = .start.bootcom, y = y, nx = nx, trimming = trimming,
			                            per=per, bootstrap = bootstrap, numsim_b = numsim_b, effect.size = effect.size,
			                            numsim_es = numsim_es, standardizer = standardizer, scaling = scaling, alpha = alpha, 
			                            seed = seed, Umat= NULL)
			  }
			  else{
					all.results = lapply(all.C.matrices, FUN = wjglm, y = y, nx = nx,
						trimming = trimming, bootstrap = bootstrap, seed = seed, numsim_es = numsim_es,
						per = per, standardizer = standardizer, scaling = scaling, numsim_b = numsim_b,
						alpha = alpha, effect.size = effect.size, Umat = NULL)
			  }
			  names(all.results) = names.all.C
				for(i in 1:length(all.results)){
					all.results[[i]]$effect.name = names.all.C[i]
				}
			}
		}
		else{
			## ----------------------------------------	
			## Only within-subject effects
			## ----------------------------------------	
		  if(bootstrap){
		    # First part of bootcom (compute the test statistic for this matrix combination using bootstrap)
		    all.results = lapply(all.U.matrices, FUN = .start.bootcom, y = y, nx = nx, trimming = trimming,
		                         per=per, bootstrap = bootstrap, numsim_b = numsim_b, effect.size = effect.size,
		                         numsim_es = numsim_es, standardizer = standardizer, scaling = scaling, alpha = alpha, 
		                         seed = seed, Cmat= NULL)
		  }
		  else{
				all.results = lapply(all.U.matrices, FUN = wjglm, y = y, nx = nx,
					trimming = trimming, bootstrap = bootstrap, seed = seed,
					per = per, standardizer = standardizer, scaling = scaling,
					numsim_b = numsim_b, numsim_es = numsim_es,
					alpha = alpha, effect.size = effect.size, Cmat = NULL)
		  }				
			for(i in 1:length(all.results)){
				all.results[[i]]$effect.name = names.all.U[i]
			}
		  names(all.results) = names.all.U
		}
		
		if(bootstrap){
		  ## ----------------------------------------	
	    ##  Compute the single bootstrap critical value for the whole test
		  ## ----------------------------------------	
	    all.fmat = lapply(all.results, `[[`, "fmat")
	    
	    fmat = do.call(rbind, args = all.fmat)
	    critv = .end.bootcom(fmat, alpha = alpha, numsim_b = numsim_b) # critical value
	    attr(all.results, "critv") = critv
	    
	    for(i in 1:length(all.results)){ # enhance the list by setting fmat to NULL
	      all.results$fmat = NULL
	    }
	    
		  mat = sapply(all.results, FUN = function(object, effect.size, critv){ 
		    if(effect.size){
		      c(object$welch.T, object$numeratorDF, object$denominatorDF, as.integer(object$welch.T > critv), object$effsz)
		    }
		    else{
		      c(object$welch.T, object$numeratorDF, object$denominatorDF, as.integer(object$welch.T > critv))
		    }
		    
		  }, effect.size = effect.size, critv = critv)
		  mat = t(mat)
		  
		  # With bootstrap we do not use p-values
		  if(effect.size){
		    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF", "significant?", "eff.size")
		  }
		  else{
		    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF", "significant?")
		  }
		  
		  rownames(mat) = sapply(all.results, FUN = "[[", "effect.name")
		  attr(all.results, "summary") = mat	
		}
		else{
		  ## ----------------------------------------	
		  ##  Adjust the p-values for multiple pariwise comparisons
		  ## ----------------------------------------			
		  all.pvalues = sapply(all.results, "[[", "pval")
		  
		  adj.pvalues = p.adjust(all.pvalues, method = correction)
		  for(i in 1:length(all.results)){
		    all.results[[i]]$adj.pval = adj.pvalues[i]
		  }
		  mat = sapply(all.results, FUN = function(object, effect.size){
		      if(effect.size){
		        c(object$welch.T, object$numeratorDF, object$denominatorDF, object$pval, object$adj.pval, object$effsz)
		      }
	        else{
	          c(object$welch.T, object$numeratorDF, object$denominatorDF, object$pval, object$adj.pval)
	        }
		  }, effect.size = effect.size)
		  
		  mat = t(mat)
		  
		  if(effect.size){
		    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval", "adj.pval", "eff.size")
		  }
		  else{
		    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval", "adj.pval")
		  }
		  
		  rownames(mat) = sapply(all.results, FUN = "[[", "effect.name")
		  attr(all.results, "summary") = mat			  
		}
		class(all.results) = "WelchTestObjList"			
		
		attr(all.results, "type") = "all.pairwise"
		attr(all.results, "correction") = correction
		attr(all.results, "effect") = effect
		attr(all.results, "bootcom") = bootstrap
		
		return(all.results)
}

# ______________________________________________________

# ______________________________________________________
#
# Generates a list of row vectors of length n, each with all components set to 0
# except for two, one of them to 1 and the other to -1
# Example: v = .generateAllPairwiseCombinations(3) yields
# list([0, -1, 1], [1,0,-1], [1, -1, 0])
# ______________________________________________________
.generateAllPairwiseCombinations <- function(n, levelnames){
	vecs = NULL
	if(n == 2){
		vecs = t(c(1, -1))
		colnames(vecs) = levelnames
		temp = list()
		temp[[1]] = vecs
		vecs = temp
	}
	else{
		a = contr.sum(n)
		vecs = NULL
		mycolumns = split(a, col(a))
		for(i in 1:length(mycolumns)){ names(mycolumns[[i]]) = levelnames }
		## transpose to turn them into row vectors
		vecs = c(vecs, lapply(mycolumns, FUN = t)) 
		for(i in (n-1):2){
			a[i,] = a[i+1,]
			a[i+1,] = 0
			a = a[,1:(i-1)]
			if(is.vector(a)){
				names(a) = levelnames
				## transpose to turn them into row vectors
				vecs = c(vecs, list(t(a)))
			}
			else{
				mycolumns = split(a, col(a))
				for(i in 1:length(mycolumns)){  names(mycolumns[[i]]) = levelnames }
				## transpose to turn them into row vectors
				vecs = c(vecs, lapply(mycolumns, FUN = t))
			}
		}	
	}
	return(vecs)
}

# ______________________________________________________

.kronecker.all <- function(alist){
	result = 1
	for(vec in alist){
		result = result %x% vec
	}
	return(result)
}
