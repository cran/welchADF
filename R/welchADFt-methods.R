##' Summarizing Welch Approximate Degrees of Freedom test results
##' 
##' \code{summary} and \code{print} methods for class \code{"welchADF"}
##' @rdname summary.welchADF
##' @param object an object of class "welchADFt" returned by a call to 
##'   \code{\link{welchADF.test}}.
##' @param verbose whether the summary table should be preceded by an extended description of the
##'   welch ADF test performed, including the kind of means (trimmed or not), kind of test
##'   (omnibus or pairwise) and whether bootstrap was used or not. Defaults to FALSE.
##' @param digits the number of significant digits to use when printing. 
##' @return a string which summarizes the info, with class \code{"summary.welchADFt"}
##'   ready to be printed by specific method \code{print.summary.welchADFt}.
##' @export
summary.welchADFt<-function(object, verbose = FALSE, digits = max(4, getOption("digits")), ...){
  
  # Save old digits configuration, set 4 and restore it
  old.digits = getOption("digits")
  options(digits = 4)
  
  type = attr(object, "type")
  bootstrap = attr(object, "bootstrap")
  effect.size = attr(object, "effect.size")
  
  msg = ""
  
  if(verbose){
    # Initial information about the kind of analysis conducted
    msg = format(object, summary = TRUE) # summary=TRUE prevents the pvalues from being included
  }
  else{
    # Only the call but no explanation
    msg = paste0(msg, "Call:\n")
    callstring = do.call(paste, args = c(as.list(deparse(object$call)), sep = "\n"))
    msg = paste0(msg, "   ", callstring, "\n\n")
  }
  
  # filter out properties that are not effects, like $call etc
  mylevels = Filter(function(x) !is.null(x[["welch.T"]]), object) 
  
  if(type == "omnibus"){
    
    mat = sapply(mylevels, FUN = function(object) {
      # in case we have not computed effect size, object$effsz will be NULL and will not be shown
      c(object$welch.T, object$numeratorDF, object$denominatorDF, object$effsz, object$pval, -1)
    })
    mat = data.frame(t(mat))
    effsz.colname = NULL
    if(effect.size){
      effsz.colname = "eff.size"
    }
    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF", effsz.colname, 
                      ifelse(bootstrap, "Pr(>WJ)!", "Pr(>WJ)"), "")
    rownames(mat) = sapply(mylevels, FUN = "[[", "effect.name")
    pval.colname = ifelse(bootstrap, "Pr(>WJ)!", "Pr(>WJ)" )
    mat[[5]] = symnum(mat[[pval.colname]], corr = FALSE, na = FALSE, 
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " "))
    msg = paste0(msg, do.call(paste, args = c(as.list(capture.output(mat)), sep = "\n")))
    msg = paste0(msg, "\n---\n")
    msg = paste0(msg, "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    
    if(bootstrap){
      msg = paste0(msg, "\n(!) Showing bootstrapped p-values\n")
    }
  }
  else if(type == "all.pairwise"){
    
    # Common columns
    mat = sapply(mylevels, FUN = function(object) c(object$welch.T, object$numeratorDF, object$denominatorDF))
    mat = data.frame(t(mat))
    colnames(mat) = c("WJ statistic", "Numerator DF", "Denominator DF")
    
    if(effect.size){
      effsz = sapply(mylevels, FUN = "[[", "effsz")
      mat[["eff.size"]] = effsz
    }    
    
    if(bootstrap){
      critv = attr(object, "critv")
      significant = sapply(mylevels, function(x) c("no", "yes")[1 + as.integer(x$welch.T > critv)])
      mat[["significant?"]] = significant
      msg = paste0(msg, do.call(paste, args = c(as.list(capture.output(mat)), sep = "\n")))
    }
    else{
      adj.pvalues = sapply(mylevels, FUN = "[[", "adj.pval")
      mat[["adj.pval"]] = adj.pvalues
      
      mat[["symbols"]] = symnum(mat[["adj.pval"]], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
      cn = colnames(mat)
      cn[cn == "symbols"] = "" # empty column name for the symbols column
      colnames(mat) = cn
      corr.method = attr(object, "correction")
      corr.method = paste0(toupper(substring(corr.method, 1, 1)), substring(corr.method, 2))
      msg = paste0(msg, do.call(paste, args = c(as.list(capture.output(mat)), sep = "\n")))
      msg = paste0(msg, "\n---\n")
      msg = paste0(msg, "Signif. codes (", corr.method, " p-values):  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
    
    if(bootstrap){
      msg = paste0(msg, "\n\nBootstrap critical value: ", format(attr(object, "critv")), "\n")
    }
  }
  class(msg) = "summary.welchADFt"
  options(digits = old.digits)
  msg
}

## ______________________________________________________________

##' @method print summary.welchADFt
print.summary.welchADFt <- function(x, ...){
  cat(x)
}

## ______________________________________________________________

##' @param x an object of class "welchADFt".
##' @param ... further arguments passed to or from other methods.
##' @rdname summary.welchADF
##' @export
print.welchADFt<-function(x, digits = getOption("digits"), ...){
  # Save old digits configuration, set 4 and restore it
  old.digits = getOption("digits")
  options(digits = digits)
  cat(format(x, ...), "\n")
  options(digits = old.digits)
}

##' @export
format.welchADFt <- function(x, ...){
  arguments = list(...)
  msg = ""
  msg = paste0(msg, "Call:\n")
  callstring = do.call(paste, args = c(as.list(deparse(x$call)), sep = "\n"))
  msg = paste0(msg, "   ", callstring, "\n\n")
  
  msg = paste0(msg, "Welch-James Approximate DF Test (")
  
  if(!x[[1]]$trimming){ msg = paste0(msg, "Least squares means & variances)\n")  }	
  else{ msg = paste0(msg, paste0("Trimmed means [", 100*x[[1]]$trimming.per, "% trimming] & Winsorized variances)\n")) }
  
  type = attr(x, "type")
  bootstrap = attr(x, "bootstrap")
  
  if(type == "omnibus"){
    msg = paste0(msg, "Omnibus test(s) of effect and/or interactions")
    if(bootstrap){
      msg = paste0(msg, " with a bootstraped critical value per effect to compute the p-values")	
    }
    msg = paste0(msg, "\n\n")
    
    # if we are calling summary(), the format function must not show the next lines
    if(is.null(arguments[["summary"]])){
      
      # Print the effect names and their p-values
      effect.names = Filter(function(x) !is.null(x), sapply(x, FUN = "[[", "effect.name"))
      pvalues = Filter(function(x) !is.null(x), sapply(x, FUN = "[[", "pval"))
      minimat = matrix(data = pvalues, nrow = 1, ncol = length(pvalues), byrow = TRUE)
      rownames(minimat) = ifelse(bootstrap, "Bootstrapped Pr(>WJ)", "Pr(>WJ)")
      colnames(minimat) = effect.names
      minimat.string = do.call(paste, args = c(as.list(capture.output(minimat)), sep = "\n"))
      minimat.string = gsub("\"", "", minimat.string)
      msg = paste0(msg, minimat.string)
    }
  }
  else if(type == "all.pairwise"){
    effect = attr(x, "effect")
    if(length(effect) == 1){
      msg = paste0(msg, "Multiple pairwise comparisons of response marginal means with respect to",effect)
    }
    else{
      temp = as.list(effect)				
      temp$sep = " x "
      temp = do.call(paste, temp)
      msg = paste0(msg, "Multiple tetrad interaction contrasts with respect to ", temp, " interaction")
    }
    if(bootstrap){
      msg = paste0(msg, "\nusing a Bootstrap Critical Value for FWER control\n")
    }
    else{
      correction = attr(x,"correction")
      correction= paste(toupper(substring(correction, 1, 1)), substring(correction, 2),
                        sep = "", collapse = " ")
      msg = paste0(msg, "\nusing ",correction," method for adjusting the p-values\n")
    }
    effect.size = attr(x, "effect.size")
    if(effect.size){
      msg = paste0(msg, "Effect size ")
      if(x[[1]]$standardize.effsz){
        msg = paste0(msg, "standardized via square root of the average of cell variances\n")
      }
      else{
        msg = paste0(msg, "with no standardization\n")
      }
    }
    # if we are calling summary(), the format function must not show the next lines
    if(is.null(arguments[["summary"]])){
      
      # Print the effect names and their p-values
      mylevels = Filter(function(x) !is.null(x$welch.T), x)
      significance = NULL
      rows.minimat = ifelse(attr(x, "effect.size"), 2, 1)
      rnames = rep("", rows.minimat)
      if(attr(x, "bootstrap")){ # there are no p-values
        significance = sapply(mylevels, FUN = function(object, critv){ 
          ifelse(object$welch.T > critv, "yes", "no")
        }, critv = attr(x, "critv"))
        rnames[1] = "significant?"
      }
      else{
        significance = sapply(mylevels, FUN = "[[", "adj.pval")
        rnames[1] = "adj.pval"
      }
      
      effect.size.v = NULL
      if(attr(x, "effect.size")){
        effect.size.v = sapply(mylevels, FUN = "[[", "effsz")
        rnames[2] = "eff.size"
      }
        
      minimat = matrix(data = c(significance, effect.size.v), 
                       nrow = rows.minimat, ncol = length(significance), byrow = TRUE)
      rownames(minimat) = rnames
      colnames(minimat) = sapply(mylevels, FUN = "[[", "effect.name")
      minimat.string = do.call(paste, args = c(as.list(capture.output(minimat)), sep = "\n"))
      minimat.string = gsub("\"", "", minimat.string)
      
      msg = paste0(msg, minimat.string)
    }
  }
  msg
}

## ______________________________________________________________

formula.welchADFt <-function(x, ...){
  mycall = getCall(x)
  if(class(eval(mycall$formula)) == "data.frame"){
    # the call was done to welchADF.test.data.frame so there is no formula
    NULL
  }
  else formula(mycall, ...)
}

## ______________________________________________________________

##' Confidence interval for the effect size of model effects
##' 
##' For the effects indicated by the user, extracts confidence intervals 
##' for the effect size from an object of class \code{welchADFt}
##' returned by a call to \code{\link{welchADF.test}}.
##' @param parm a specification of which parameters are to be given confidence intervals, 
##'   either a vector of numbers or a vector of names. If missing, all parameters are considered.
##' @param object a object of (S3 class) welchADF.t
##' @param level confidence level (completely ignored, since the confidence threshold 
##'   had been in the call to \code{\link{welchADF.test}}).
##' @param ... additional argument(s) for methods.
##' @return If  \code{effect.size} was set to \code{TRUE} in the call to \code{welchADF.test}, 
##'   a matrix (or vector) with columns giving lower and upper confidence limits for each parameter. 
##'   These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
##'   If \code{effect.size} was set to \code{FALSE}, \code{NULL} is returned.
##' @seealso \code{\link{welchADF.test}}, \code{\link[stats]{confint.lm}}
confint.welchADFt <- function(object, parm, level = 0.95, ...){
  if(!attr(object, "effect.size")){
    NULL
  }
  else{
    pnames = unlist(sapply(object, "[[", "effect.name"))
    mylevel = unique(unlist(sapply(object, "[[", "CI.level")))/100
    if(!missing(level)){
      warning(paste0("level argument ignored. Setting level = ", mylevel))
    }
    level = mylevel
    ci.pairs = sapply(object, "[[", "CI")
    all.params = missing(parm)
    if(all.params){
      parm = pnames
    }
    else if(is.numeric(parm)){
      parm = pnames[parm]
    }
    
    if(sum(!parm %in% pnames) > 0){
      notfound.indices = which(is.na(match(parm, pnames)))
      notfound = parm[notfound.indices]
      notfound = do.call(paste, args = as.list(c(notfound, sep = ", ")))
      msg = do.call(paste0, args = list(c("effect(s) ", notfound, " not found")))
      stop(msg)
    }
    # ci has dimensions length(parm) x 2 
    ci = matrix(unlist(ci.pairs[parm]), ncol = 2, byrow = TRUE)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    rownames(ci) = parm
    colnames(ci) = pct
    ci
  }
}

## ______________________________________________________________

model.frame.welchADFt <- function(formula, ...){
  if(is.null(formula(formula))){ 
    # the welchADFt object in formula was returned by a call to welchADF.test.data.frame
    mycall = getCall(formula)
    eval(mycall$formula)
  }
  else{
    # the welchADFt object in formula was returned by a call to welchADF.test.formula or similar
    formula$model
  }
}
