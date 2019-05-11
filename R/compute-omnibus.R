
# argument 'only' is a string vector with the effects and combinations specified by the user 
# in the formula interface of welchADF.test. E.g. only = c("a", "b", "c", "a:b", "a:c")
.compute.omnibus <- function(data, response, between.s, within.s, subject, effect,
				trimming, per, bootstrap, numsim_b, effect.size, numsim_es, standardizer, 
				scaling, alpha, seed, y, nx, levelslist.between.s, levelslist.within.s, only){

	frameNames = names(data)	
	ncols = ncol(data)
	factorColumns = (1:ncols)[!(frameNames %in% response)]
	columns.between.s = (1:ncols)[(frameNames %in% between.s)]
	columns.within.s = (1:ncols)[(frameNames %in% within.s)]
	
	all.C.matrices = NULL
	all.U.matrices = NULL
	unique.C.matrix = NULL
	unique.U.matrix = NULL
	
	no.between.s.matrix = 1
	no.within.s.matrix = NULL
	all.between.s.results = NULL
	all.within.s.results = NULL
	all.mixed.results = NULL

	if(length(response) > 1 && is.null(within.s)){
	  no.within.s.matrix = diag(length(response))
	}
				
	if(!is.null(between.s)){ 
		## there exist between-subjects factors
		## ------------------------------------------------
		## Main effects and interactions involving only between-subjects effects
		## ------------------------------------------------
		rhs = as.list(frameNames[columns.between.s])
		rhs$sep = "*"
		formulastring = do.call(paste, rhs) # all effects and interactions
		myformula = as.formula(paste0("y ~ ", formulastring))			
		myterms = terms(myformula)
		terms.matrix.between.s = attr(myterms, "factors")		

		if(is.null(effect)){
			all.C.matrices = .compose.C.matrices(terms.matrix.between.s, levelslist.between.s)
			no.between.s.matrix = all.C.matrices[[1]]
			all.C.matrices = all.C.matrices[-1]	# drop first matrix corresponding to no between-subject effect
			if(!is.null(only)){
			  # Filter out those between-subjects factors not present in only
			  validMatrices = names(all.C.matrices) %in% only
			  all.C.matrices = all.C.matrices[validMatrices]
			}
		}
		else{
			# Select from terms.matrix.between.s the specific between-subject terms involved in the effect argument				
			effect.between.s = effect[effect %in% frameNames[columns.between.s]]
			if(length(effect.between.s) > 0){
				# Now in the user-given effects vector, put the factors in the same order of the dataset
				positions = match(rhs, effect.between.s)
				positions = positions[!is.na(positions)]
				effect.between.s = as.list(effect.between.s[positions])
				effect.between.s$sep = ":"
				between.s.formulastring = do.call(paste, effect.between.s)
				terms.matrix.between.s = t(t(terms.matrix.between.s[,between.s.formulastring])) # select only the column of this effect
				colnames(terms.matrix.between.s) = between.s.formulastring

				all.C.matrices = .compose.C.matrices(terms.matrix.between.s, levelslist.between.s)
				unique.C.matrix = all.C.matrices[[2]]	# drop first matrix corresponding to no between-subject effect 
			}
			else{
				# there exist between-subject factors but they are not involved in the effect queried by the user
				all.C.matrices = .compose.C.matrices(terms.matrix.between.s, levelslist.between.s)
				unique.C.matrix = all.C.matrices[[1]]
			}
		}
	}
	
	if(!is.null(within.s)){
		## ------------------------------------------------
		## Main effects and interactions involving only within-subjects effects
		## ------------------------------------------------
		rhs = as.list(names(data)[columns.within.s])
		rhs$sep = "*"
		formulastring = do.call(paste, rhs) # all effects and interactions
		myformula = as.formula(paste0("y ~ ", formulastring))			
		myterms = terms(myformula)
		terms.matrix.within.s = attr(myterms, "factors")
		if(is.null(effect)){  	
			all.U.matrices = .compose.U.matrices(terms.matrix.within.s, levelslist.within.s, length(response))
			no.within.s.matrix = all.U.matrices[[1]]
			all.U.matrices = all.U.matrices[-1]	# drop first matrix corresponding to no within-subject effect
			if(!is.null(only)){
			  # Filter out those within-subjects factors not present in only
			  validMatrices = names(all.U.matrices) %in% only
			  all.U.matrices = all.U.matrices[validMatrices]
			}
		}
		else{
			# Select from terms.matrix.within.s the specific within-subject terms involved in the effect argument				
			effect.within.s = effect[effect %in% frameNames[columns.within.s]]
			if(length(effect.within.s) > 0){
				# Now in the user-given effects vector, put the factors in the same order of the dataset
				positions = match(rhs, effect.within.s)
				positions = positions[!is.na(positions)]
				effect.within.s = as.list(effect.within.s[positions])
				effect.within.s$sep = ":"
				within.s.formulastring = do.call(paste, effect.within.s)
				terms.matrix.within.s = t(t(terms.matrix.within.s[,within.s.formulastring])) # select only the column of this effect
				colnames(terms.matrix.within.s) = within.s.formulastring

				all.U.matrices = .compose.U.matrices(terms.matrix.within.s, levelslist.within.s, length(response))
				#no.within.s.matrix = all.U.matrices[[1]]		  	
				unique.U.matrix = all.U.matrices[[2]]	# drop first matrix corresponding to no within-subject effect 
				
			}
			else{
				all.U.matrices = .compose.U.matrices(terms.matrix.within.s, levelslist.within.s, length(response))
				unique.U.matrix = all.U.matrices[[1]]
			}
		}
	}
		
	if(is.null(effect)){
		## ------------------------------------------------			
		## Omnibus test of ALL main and interaction effects
	  ## (if argument 'only' is not null, test only those indicated on it)
		## ------------------------------------------------
		if(!is.null(between.s)){
			## ------------------------------------------------
			## Test all between-subjects effects assuming no within-subject effect
			## ------------------------------------------------
			all.between.s.results = lapply(all.C.matrices, FUN = wjglm, y = y, nx = nx,
					trimming = trimming, bootstrap = bootstrap, seed = seed, 
					standardizer = standardizer, scaling = scaling, alpha = alpha,
					numsim_b = numsim_b, numsim_es = numsim_es,
					per = per, effect.size = effect.size, Umat = no.within.s.matrix)
			for(i in 1:length(all.between.s.results)){
				all.between.s.results[[i]]$effect.name = names(all.between.s.results[i])
			}
		}
		if(!is.null(within.s)){
			## ------------------------------------------------
			## Test all within-subjects effects assuming no between-subject effect
			## ------------------------------------------------
			all.within.s.results = lapply(all.U.matrices, FUN = wjglm, y = y, nx = nx,
					trimming = trimming, bootstrap = bootstrap, seed = seed, 
					standardizer = standardizer, scaling = scaling, alpha = alpha,
					numsim_b = numsim_b, numsim_es = numsim_es,
					per = per, effect.size = effect.size, Cmat = no.between.s.matrix)
			for(i in 1:length(all.within.s.results)){
				all.within.s.results[[i]]$effect.name = names(all.within.s.results[i])
			}
		}
		
		if(!is.null(between.s) && !is.null(within.s)){
			## ------------------------------------------------
			## Test all mixed interactions between between- and within-subject effects
			## ------------------------------------------------
			all.mixed.results = vector("list", length(all.C.matrices)*length(all.U.matrices))
			i = 1
			#for( matC in all.C.matrices){
			#	for( matU in all.U.matrices){
			for( x1 in 1:length(all.C.matrices)){
			  for( x2 in 1:length(all.U.matrices)){
			    validCombination = TRUE
			    combinationName = paste(names(all.C.matrices)[x1], names(all.U.matrices)[x2], sep = ":")
			    if(!is.null(only)){
			      # check if this combination of factors is among those specified by the user
			      validCombination = combinationName %in% only
			    }
			    if(validCombination){
  				  matC = all.C.matrices[[x1]]
  				  matU = all.U.matrices[[x2]]
  					all.mixed.results[[i]] = wjglm(Cmat = matC, Umat = matU, y = y, nx = nx,
  						trimming = trimming, bootstrap = bootstrap, seed = seed, 
  						standardizer = standardizer, scaling = scaling, 
  						numsim_b = numsim_b, numsim_es = numsim_es,
  						per = per, alpha = alpha, effect.size = effect.size)
  					
  					# add space to both sides of ":"
  					all.mixed.results[[i]]$effect.name = sub(":", " : ", combinationName)
  					i = i+1
			    }
				}
			}
		}
	
		result = c(all.between.s.results, all.within.s.results, all.mixed.results)
	}
	else{
	  # the result will be a list of one element
		element = wjglm(Cmat = unique.C.matrix, Umat = unique.U.matrix, y = y, nx = nx,
					trimming = trimming, bootstrap = bootstrap, seed = seed, alpha = alpha,
					numsim_b = numsim_b, numsim_es = numsim_es,
					per = per, standardizer = standardizer, scaling = scaling, effect.size = effect.size)
		
		efflist = as.list(effect)
		efflist$sep = ":"						
		element$effect.name = do.call(paste, efflist)
		
		result = list(element)
	}
	
	class(result) = "welchADFt"
	attr(result, "type") = "omnibus"
	attr(result, "bootstrap") = bootstrap
	attr(result, "effect.size") = effect.size
			
	return(result)		
}

# ______________________________________________________

# ======================================================
## Returns a list of C matrices, one for an omnibus test of each between-subject single effect or between-effect interaction appearing in
## the columns of formula.matrix, plus one additional matrix in the first position of the
## list with the C matrix for no-between-subjects effects
# ______________________________________________________
.compose.C.matrices <- function(formula.matrix, levelslist){
		
  n.mainefx = nrow(formula.matrix) - 1 # number of main effects (first row is for the response variable)
  n.efx = ncol(formula.matrix)
	result = vector("list", n.efx + 1) # the first position is for the C matrix where no between-subject effect is involved
	names(result) = c("<none>", colnames(formula.matrix))
	result[1:length(result)] = 1
	
	# drop row 1 which is all 0 corresponding to the non-existent y
	formula.matrix = formula.matrix[-1,] 	
	if(!is.matrix(formula.matrix)){
		formula.matrix = t(t(formula.matrix))
		names(formula.matrix) = names(result)[-1] # drop the first name
	}
			
	nlevelsvec = sapply(levelslist, FUN = length)
	
	## generate matrices for all linearly independent contrasts
	mainC = lapply(nlevelsvec, FUN = function(x) -t(contr.sum(x))) 	
	names(mainC) = names(levelslist)
	for(i in 1:length(mainC)){
		colnames(mainC[[i]]) = levelslist[[i]]
		rownames(mainC[[i]]) = levelslist[[i]][1:(nlevelsvec[i]-1)]
	}	
	
	if(is.matrix(formula.matrix)){
		for(i in 1:n.mainefx){
			for(j in 1:n.efx){
				if(formula.matrix[i,j]){ 		result[[j+1]] = result[[j+1]] %x% mainC[[i]]											}
				else{ 											result[[j+1]] = result[[j+1]] %x% t(rep(1,nlevelsvec[i]))				}
			}
			result[[1]] = result[[1]] %x% t(rep(1,nlevelsvec[i]))
		}		
	}else{
		if(formula.matrix){		result[[1]] = mainC[[1]] 		}
		else{									result[[1]] = t(rep(1,nlevelsvec[1])) }
	}

	return(result)
}

# ______________________________________________________

# ======================================================
## Returns a list of U matrices, one for an omnibus test of each within-subject single effect or within-effect interaction appearing in
## the columns of formula.matrix, plus one additional matrix in the first position of the
## list with the U matrix for no-within-subjects effects
# ______________________________________________________
.compose.U.matrices <- function(formula.matrix, levelslist, number.of.responses){

  n.mainefx = nrow(formula.matrix) - 1 # number of main effects (first row is for the response variable)
  n.efx = ncol(formula.matrix)
	result = vector("list", n.efx + 1) # the first position is for the U matrix where no between-subject effect is involved
	names(result) =  c("<none>", colnames(formula.matrix))
	result[1:length(result)] = 1
		
	# drop row 1 which is all 0 corresponding to the non-existent y
	formula.matrix = formula.matrix[2:nrow(formula.matrix),] 	
	if(!is.matrix(formula.matrix)){
		formula.matrix = t(t(formula.matrix))
		names(formula.matrix) = names(result)[2:length(names(result))]
	}
			
	nlevelsvec = sapply(levelslist, FUN = length)
	
	## generate matrices for all linearly independent contrasts
	mainU = lapply(nlevelsvec, FUN = function(x) -(contr.sum(x))) 	
	names(mainU) = names(levelslist)
	for(i in 1:length(mainU)){
		rownames(mainU[[i]]) = levelslist[[i]]
		colnames(mainU[[i]]) = levelslist[[i]][1:(nlevelsvec[i]-1)]
	}	
	
	if(is.matrix(formula.matrix)){
		for(i in 1:n.mainefx){
			for(j in 1:n.efx){
				if(formula.matrix[i,j]){ 		result[[j+1]] = result[[j+1]] %x% t(mainU[[i]])									}
				else{ 											result[[j+1]] = result[[j+1]] %x% t(rep(1,nlevelsvec[i]))				}
			}
			result[[1]] = result[[1]] %x% t(rep(1,nlevelsvec[i]))
		}
	}else{
		if(formula.matrix){		result[[1]] = t(mainU[[1]]) 		}
		else{									result[[1]] = t(rep(1,nlevelsvec[1])) }
	}
	
	if(number.of.responses > 1){
		result = lapply(result, FUN = `%x%`, Y = diag(number.of.responses))
	}
	
	result = lapply(result, FUN = t)	
	
	return(result)	
}

# ______________________________________________________