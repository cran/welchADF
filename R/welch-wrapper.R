##' General Approximate Degrees of Freedom Solution for Inference and Estimation
##' 
##' Computes the Welch-James statistic based on the SAS code by H.J. Keselman, R.R. Wilcox and L.M. Lix.
##' @param data Data frame with one row for each observation of the response variable(s), and columns to indicate the factor combination to which it corresponds.
##' @param response A string or vector of strings with the name(s) of the column(s) of \code{data} corresponding to the response variable(s). If 
##'		a vector of strings is provided, the responses are taken as a set of repeated measurements or dependent variables.
##' @param between.s Vector of strings with the columns that correspond to between-subjects factors.
##' @param within.s Vector of strings with the columns that correspond to within-subjects factors, if any (defaults to NULL).
##' @param subject Name of the column (if any) containing the subject ID (defaults to NULL).
##' @param contrast One of \code{"omnibus", "all.pairwise"} indicating the type of contrast to be performed with the data.
##' \code{"omnibus"} stands for omnibus contrasts on either all the variables found in the data, including the within-subjects factor if any,
##' or on a single effect specified in the \code{effect} argument.
##' \code{"all.pairwise"} stands for pairwise comparisons of the response(s) variable(s) over all combinations of values of one effect, indicated
##'	 in the \code{effect} parameter.
##'	Defaults to \code{"omnibus"}.
##'	@param effect If \code{contrast = "omnibus"}, then \code{effect} must be a string or vector of strings with the name of the effect(s) for which the \code{"omnibus"} 
##'		contrasts will be applied. 
##'		If \code{contrast = "all.pairwise"}, then \code{effect} must be \emph{a tagged list} in which the names are the factors, and the elements must be vectors of two strings
##'		each, indicating 
##'		If \code{contrast = "omnibus"} but \code{effect} is not specified, then an omnibus contrast will be done separately for each of the effects.
##'		If \code{contrast = "omnibus"} and \code{effect} is a vector of strings, then only \emph{one} contrast will be done, involving the simultaneous interaction of all effects
##'		indicated in the \code{effect} vector (i.e. no other interaction or main effect contrasts will be performed).
##' @param correction The type of p-value correction when applying multiple pariwise comparisons (i.e. when \code{contrast = "all.pairwise"}). Defaults to "hochberg".
##' @param trimming Boolean to control the use of robust estimators. 
##'		\code{FALSE} indicates the ADF test statistic will be computed with the usual Least-Squares Estimators. 
##' \code{TRUE} indicates that trimmed means and Winsorized variances/covariances will be used to compute the ADF test statistic.
##' @param per (Only if trimming is TRUE) The proportion of trimming in EACH tail of the data distribution. Must be between 
##'   0 and .49 (i.e., 49\% trimming in each tail). Recommended per = 0.2 symmetric trimming (i.e., 20\% of the observations in each tail are trimmed).
##' @param bootstrap Boolean; whether bootstrap should be used to compute a critical value for the test statistic produced by WJGLM. \code{FALSE} means a theoretical 
##' critical value for the test statistic TWJ/c (or TWJt/c) will be computed. \code{TRUE} indicates the non-parametric percentile bootstrap method 
##' will be used to compute a critical value. When combined with \code{contrast = "all.pairwise"}, then Family-Wise Error Rate (FWER) control is achieved using
##' bootstrapping to compute an adjusted empirical critical value (refer to Keselman, Wilcox & Lix 2003, page 589 for details).
##' @param numsim_b If \code{bootstrap = TRUE}, this is a positive integer that defines the number of bootstrap samples to generate a critical 
##' value (defaults to 999).  Carpenter and Bithell (2000) recommend between 1000 and 2000 samples for setting confidence intervals around a 
##' parameter estimate, but Efron and Tibshirani (1986) note that for other applications of the bootstrap, such as computing the standard error, as 
##' few as 100 bootstrap samples will be sufficient. Choosing a value of B such that (B + 1)\eqn{\alpha} is an integer value is recommended, because it 
##' avoids the need for interpolation.
##' @param seed Initial seed (positive integer) for the (R default) random number generator. Defaults to NULL (seed taken from the current time). 
##' @param effect.size Boolean; whether effect size estimates should be computed. \code{TRUE} means effect size estimates together with its confidence intervals 
##' will be computed (defaults to FALSE).
##' @param numsim_es Positive integer defining the number of bootstrap samples to generate a CI for the effect size estimate (defaults to 999). Ignored if \code{effect.size = FALSE}.
##' @param alpha Significance threshold for the confidence interval calculation, where 1 - \code{alpha} is the probability coverage. Defaults to 0.05.
##' @param scaling Boolean; whether a scaling factor for the effect size estimator should be used (0.642 for 20\% symmetric trimming) when 
##' robust estimators are adopted. \code{FALSE} means no scaling factor will be used. \code{TRUE} indicates that a scaling factor will be adopted. Defaults to TRUE
##' @param standardize.effsz Boolean: whether the effect size should be standardized by the average of variances or not. Defaults to TRUE.
##' @references Home page of Prof. Lisa M. Lix, author of the original SAS code: \url{http://homepage.usask.ca/~lml321/SAS_Programs.html}
##' @references Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach to hypothesis testing in independent and correlated groups designs. 
##' \emph{Psychophysiology}, 40, 586-596.
##' @references Lix, L. M., & Keselman, H. J. (1995). Approximate degrees of freedom tests: A unified perspective on testing for mean equality. 
##' \emph{Psychological Bulletin}, 117, 547-560. 
##' @references Carpenter, J., & Bithell, J. (2000). Bootstrap confidence intervals: When, which, what? A practical guide for medical statisticians. 
##' \emph{Statistics in Medicine}, 19, 1141-1164. 
##' @references Efron, B., & Tibshirani, R. (1986). Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. 
##' \emph{Statistical Science}, 1, 54-75. 
##' @seealso \code{\link{p.adjust.methods}}
##' @examples
##' # Two-way factorial design. See the vignette
##' omnibus_LSM = welchADF.test(womenStereotypeData, response = "y", between.s = 
##' c("condition", "sex"), contrast = "omnibus")
##' omnibus_trimmed = welchADF.test(womenStereotypeData, response = "y", between.s = 
##' c("condition", "sex"), contrast = "omnibus", trimming = TRUE)
##' pairwise_LSM = welchADF.test(womenStereotypeData, response = "y", between.s = 
##' c("condition", "sex"), contrast = "all.pairwise", effect = c("condition", "sex"))
##' pairwise_trimmed = welchADF.test(womenStereotypeData, response = "y", between.s = 
##' c("condition", "sex"), contrast = "all.pairwise", effect = c("condition", "sex"), 
##'   trimming = TRUE)
##' pairwise_trimmed_boot = welchADF.test(womenStereotypeData, response = "y", between.s = 
##'   c("condition", "sex"), contrast = "all.pairwise", effect = c("condition", "sex"), 
##'   trimming = TRUE, bootstrap = TRUE)
##' summary(omnibus_LSM)
##' summary(pairwise_LSM)
##' summary(pairwise_trimmed_boot)
welchADF.test <- function(data, response, between.s, within.s = NULL, subject = NULL, contrast = c("omnibus", "all.pairwise"), 
			effect = NULL, correction = c("hochberg", "holm"), trimming = FALSE, per=0.2, bootstrap = FALSE, 
			numsim_b = 999, effect.size = FALSE, numsim_es = 999, scaling = TRUE, standardize.effsz = TRUE, alpha = 0.05, seed = 0){
	
	contrast = match.arg(contrast)
	correction = match.arg(correction)
	
	## ----------------------------------------
	##		CHECK ALL PARAMETERS ARE VALID
	## ----------------------------------------
	.check.parameters(data, response, between.s, 
	                  within.s, subject, effect, 
	                  contrast, correction, bootstrap, trimming,
	                  scaling, standardize.effsz, effect.size, per, 
	                  numsim_b, numsim_es, alpha)
	## ----------------------------------------

	if(effect.size && (contrast != "all.pairwise")){
	  effect.size = FALSE
	}
	
	if(!is.null(within.s)){
  	if(within.s == "multivariate"){
  	  # turn several response columns into one, adding an explicit within-suject column
  	  templist = .reshape.implicit.withins(data, response, subject)
  	  data = templist$data
  	  response = templist$response
  	  within.s = templist$within.s
  	  subject = templist$subject
  	}
	}
	
	frameNames = names(data)	
	ncols = ncol(data)
	factorColumns = (1:ncols)[!(frameNames %in% response)]
	columns.between.s = (1:ncols)[(frameNames %in% between.s)]
	columns.within.s = (1:ncols)[(frameNames %in% within.s)]
	
	
	## turn into factor those columns that are not part of the response columns and are not factors yet
	for(i in 1:length(factorColumns)){
		if(!is.factor(data[[factorColumns[i]]])){
			data[[factorColumns[i]]] = as.factor(data[[factorColumns[i]]])
		}
	}
	
	## ------------------------------------------
	if(!is.null(subject)){
		if(!is.factor(subject)){
			data[[subject]] = as.factor(data[[subject]])
		}
	}

	nlevelslist = lapply(data[factorColumns], FUN = nlevels)
	
	formatted.data = .reshape.data(data, response, between.s, within.s, subject, nlevelslist)	
	y = formatted.data[[1]]
	nx = formatted.data[[2]]

  levelslist.between.s = lapply(data[columns.between.s], FUN = levels)  
	levelslist.within.s = lapply(data[columns.within.s], FUN = levels)

	if(contrast == "omnibus"){
		## ----------------------------------------
		##					OMNIBUS CONTRASTS
		## ----------------------------------------		

		omnibus.results = .compute.omnibus(
		  data = data, response = response, between.s = between.s, within.s = within.s, 
		  subject = subject, effect = effect, trimming = trimming, per = per, 
		  bootstrap = bootstrap, numsim_b = numsim_b, effect.size = effect.size, 
		  numsim_es = numsim_es, standardizer = standardize.effsz, scaling = scaling, 
		  alpha = alpha, seed = seed, y = y, nx = nx, 
		  levelslist.between.s = levelslist.between.s, levelslist.within.s = levelslist.within.s)
	
		return(omnibus.results)
	}		
	else if(contrast == "all.pairwise"){
		## ----------------------------------------
		##		ALL PAIRWISE CONTRASTS OF ONE EFFECT LEVELS
		## ----------------------------------------	
		##	NOTE: only pairwise comparisons of marginal means
		##  or 2-way interactions via tetrad contrasts are implemented
		## ----------------------------------------	
		if(bootstrap){
		  numsim_b = 699 # for FWER control via bootstrapping we set numsim_bc to 699 instead of 999
		}
	  
		pairwise.results = .compute.pairwise(
		  data = data, response = response, between.s = between.s, within.s = within.s, 
		  subject = subject, effect = effect, correction = correction, trimming = trimming, 
		  per = per, bootstrap = bootstrap, numsim_b = numsim_b, effect.size = effect.size, 
		  numsim_es = numsim_es, standardizer = standardize.effsz, 
		  scaling = scaling, alpha = alpha, seed = seed, y = y, nx = nx, 
		  levelslist.between.s = levelslist.between.s, levelslist.within.s = levelslist.within.s)
		
		return(pairwise.results)
						
	}	# if contrast == "all.pairwise"
}

# ______________________________________________________

.reshape.data <- function(data, response, between.s, within.s, subject, nlevelslist){

	frameNames = names(data)
	ncols = ncol(data)

	columns.between.s = (1:ncols)[(frameNames %in% between.s)]
	columns.within.s = (1:ncols)[(frameNames %in% within.s)]	
	subject.col = (1:ncols)[frameNames == subject]
	nsubjects = 0
	repetitions.per.subject = NULL
	if(is.null(subject)){ 
		subject = integer(0) 
		nsubjects = 0
	}
	else{
			nsubjects = nlevels(data[[subject]])
			nlevels.within = unlist(nlevelslist[within.s])
			repetitions.per.subject = sum(nlevels.within)
	}
	
	responseColumns = match(response, frameNames)
	
	# IMPORTANT: first, the between-subject factors, then the subject column, 
  # and last the within-subject factors
  
	allData = data[ , c(columns.between.s, subject.col, columns.within.s)]
	
	betweenSubjectData = data[,columns.between.s]
	tablecount = table(betweenSubjectData)
	if(nsubjects>0){
		tablecount = tablecount/repetitions.per.subject
	}
	r = as.data.frame(tablecount)
	
	if(sum(r[["Freq"]] == 0) > 0){
		stop("Data corresponding to some between-subjects factor combination(s) are missing: ")
	}

	orderedfactors = do.call(order, args = as.list(betweenSubjectData))
	if(length(between.s) == 1){
    orderedfactors = do.call(order, args = list(betweenSubjectData))
  }
	orderedcombinations = do.call(order, args = as.list(r[,1:(ncol(r)-1)]))
	if(ncol(r) == 2){	# correct for this special case
		orderedcombinations = order(r[[1]])
	}

	## ----------------------------------------------------------------
	##  COMPUTE VECTOR nx OF NUMBER OF OBSERVATIONS FOR EACH CELL
	##				(i.e. EACH BETWEEN-FACTORS COMBINATION)
	## ----------------------------------------------------------------
	r = r[orderedcombinations,]	
	cnms = colnames(r)	
	colnames(r) = cnms
	nx = r[["Freq"]]
	## ----------------------------------------------------------------

	ordered.all = do.call(order, args = as.list(allData))		
	responses.only = data[ordered.all, response]	# response can be more than one column (multivariate responses)
	y = NULL
	
	if(nsubjects > 0){	# there are repeated measures (within-subjects factors)
		y = matrix(data = as.vector(t(responses.only)), nrow = nsubjects, byrow = TRUE)
	}
	else{	# there are no repeated measures
		y = as.matrix(data[orderedfactors, response])
	}

	return(list(y, nx, r))
}

# ______________________________________________________

.check.parameters <- function(data, response, between.s, within.s, subject, effect, contrast, correction, bootstrap, trimming,
                              scaling, standardize.effsz, effect.size, per, 
                              numsim_b, numsim_es, alpha){

	frameNames = names(data)	
	if(is.null(response)){
		stop("ERROR: at least one non-null response column must always be provided")
	}
	
	if(is.null(between.s) && is.null(within.s)){
		stop("ERROR: there must be at least one between-subject or one within-subject effect. Both cannot be null simultaneously")
	}


  for(i in 1:length(response)){
    if(!(response[i] %in% frameNames)){
      stop("ERROR: response ",response[[i]]," not found in the data")
    }
  }
	
	if(length(response) == 1 && "multivariate" %in% within.s){
	  stop("ERROR: within.s vector contains \"multivariate\" but there is only one response column, hence there is no implicit within-subject factor")
	}

  if(length(effect) > 0){
		for(i in 1:length(effect)){
			if(!(effect[i] %in% frameNames) && effect[i] != "multivariate"){
				stop("ERROR: effect ", effect[i], " not found in the data")
			}
		}
  }

	if(length(between.s) > 0){
		for(i in 1:length(between.s)){
			if(!(between.s[i] %in% frameNames)){
				stop("ERROR: between-subject effect ", between.s[i], " not found in the data")
			}
		}
  }

  if(length(within.s) > 0){
		for(i in 1:length(within.s)){
			if(!(within.s[i] %in% frameNames) && within.s[i] != "multivariate"){
				stop("ERROR: within-subject effect ", within.s[i], " not found in the data")
			}
		}
  }
    
	if(contrast == "all.pairwise" && !(length(effect) == 1 || length(effect)==2)){	
			stop("ERROR: parameter effect must be a vector of 1 or 2 elements to compute all pairwise comparisons between the levels of that effect or interaction")
	}
	
	if(length(subject) > 1){
		stop("ERROR: the subject must be indicated in exactly one column, not more than one")
	}	
	
	if(!is.logical(trimming)){
	  stop("ERROR: argument trimming must be given a boolean value, either TRUE or FALSE")
	}
	if(!is.logical(bootstrap)){
	  stop("ERROR: argument bootstrap must be given a boolean value, either TRUE or FALSE")
	}	
	
	if( bootstrap && !trimming ){
	  warning("attempting to use bootstrapping without trimming. Consider setting trimming = TRUE")
	}
	
	if(!is.logical(scaling)){
	  stop("ERROR: Argument scaling must be given a boolean value, either TRUE or FALSE")
	}
	
	if(!is.logical(standardize.effsz)){
	  stop("ERROR: Argument standardize.effsz must be given a boolean value, either TRUE or FALSE")
	}
	
	if(effect.size && (contrast != "all.pairwise")){
	  warning("attempting to compute effect size in a contrast different than all.pairwise. No effect size will be calculated.")
	}
	
	if(numsim_b < 0 || numsim_es < 0){
	  stop("ERROR: Arguments numsim_b and numsim_es must take a positive integer value")
	}

	if(trimming && (per < 0 || per > 0.49)){
	  stop("ERROR: argument per must be between 0 and 0.49")
	}
	
  valid_condition1 = ( is.null(within.s) && is.null(subject) )
  valid_condition2 = ( !is.null(within.s) && !is.null(subject) )
  valid_condition3 = FALSE
  if( length(within.s) == 1 ){
    valid_condition3 = ( within.s == "multivariate" && is.null(subject) )
  }
  if(! ( valid_condition1 || valid_condition2 || valid_condition3 ) ){
    stop(paste0("ERROR: a non-null subject column must always be accompanied by at least one within-subject column name and viceversa,\n",
                "except when within.s = \"multivariate\". If there are no within-subjects effects, leave the subject argument blank"))
  }
}

# ______________________________________________________

.reshape.implicit.withins <- function(data, response, subject){
  
  mynames = colnames(data)
  newnames = sapply(mynames, FUN = function(x, response){ 
    if(x %in% response){ 
      paste0("response.",x) # so that we later call reshape() and get a single
    }                       # column called "response" and the levels of the 
    else{ x }               # new multivariate factor are each of the .x
  }, response)
  colnames(data) = newnames  
  
  varying = paste0("response.", response)  
  mysubject = subject
  if(is.null(subject)){
    # This means the multivariate is the only within-subjects effects, therefore
    # each row corresponds to a different subject
    data$subject = 1:nrow(data)
    mysubject = "subject"
  }
  
  newdata = reshape(data, varying = varying, direction = "long", timevar = "multivariate")
  rownames(newdata) = NULL
  newdata = subset(newdata, select = -get("id"))
  return(list(data = newdata, response = "response", within.s = "multivariate", subject = mysubject))
}