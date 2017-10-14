## ------------------------------------
##  helper functions to parse a formula
## ------------------------------------
toChar <- function(lang) { 
  if (length(lang) == 3) {
    return(c(toChar(lang[[2]]), as.character(lang[[3]])))
  } else if (length(lang) == 2) {
    return(as.character(lang[[2]]))
  } else {
    return(as.character(lang))
  }
}

## from lme4
RHSForm <- function(form,as.form=FALSE) {
  rhsf <- form[[length(form)]]
  if (as.form) reformulate(deparse(rhsf)) else rhsf
}

## from lme4
`RHSForm<-` <- function(formula,value) {
  formula[[length(formula)]] <- value
  formula
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

## Check parameters of welchADF.test.data.frame
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