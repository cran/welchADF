##' General Approximate Degrees of Freedom Solution for Inference and Estimation
##' 
##' Computes the Welch-James statistic with Approximate Degrees of Freedom, 
##' based on the SAS code by H.J. Keselman, R.R. Wilcox and L.M. Lix.
##'
##' @param formula A data frame or a formula or an lm object returned by \code{\link[stats]{lm}}.
##'   If \code{formula} is a data frame, it must have 
##'   the format described for argument \code{data} of \code{welchADF.test} for class formula.
##'   If \code{formula} is a formula, it must be a two-sided linear formula object describing the 
##'   between-subjects effects of the model, with the response variable(s) (surrounded by \code{cbind()} 
##'   if it is a multivariate response, as \code{cbind(r1, r2, r3)}) on the left of the ~ 
##'   operator and the terms on the right. Within-subject factors must be specified as a parenthesis
##'   term on the rhs with the form (w1*w2*...|subject) with 'subject' being the name of the 
##'   subject ID column, which must be present in the data (cannot be omitted like in the data.frame version)
##'   e.g. \code{cbind(y1,y2) ~ A*B*C + (B*C|subjet)}. In this example, A would be a between-subject
##'   factor because it is the only variable not present again on the left side of the parenthesis term
##'   (where the within-subjects factors are specified). Multi-level or customized structures are not supported, e.g.
##'   those having strata or more than one random factor. The formula does not affect the results of the tests as 
##'   no 'model' is actually fit to the data. It only affects the terms displayed as a result of an omnibus
##'   contrast when no \code{effect} argument is given (which means \emph{run an omnibus test for every effect
##'   or interaction (up to second order)}). In such case, the factors or interactions not specified in the formula
##'   will not be applied an omnibus test. In case of a pairwise contrast, the structure in the formula is ignored
##'   beyond interpreting which variables are between-subject factors and which are within-subjects factors.
##' @param ... Further arguments to be passed to specialized methods. Should be named arguments.
##' @return An object of class "welchADFt" which is a list of lists (one sub-list per effect, even if there is only one).
##'   There are methods \code{\link[welchADF]{print.welchADFt}} and \code{\link[welchADF]{summary.welchADFt}}
##' @seealso \code{\link[welchADF]{print.welchADFt}}, \code{\link[welchADF]{summary.welchADFt}}
##' @export
welchADF.test <- function(formula, ...){
  UseMethod("welchADF.test")
}

##' @describeIn welchADF.test Specialized method that accepts a formula
##' @export
##' @method welchADF.test formula
##' @param data A data.frame with the data, formatted as follows: one row per 
##'   observation of the response variable(s), and as many columns as needed to indicate the 
##'   factor combination to which the observation corresponds. If necessary, an extra column with
##'   the subject ID for designs having within-subjects factors
##'   (can be omitted if there is only one within-subjects factor, see the vignette).
##' @param subset A specification of the rows to be used as in \code{\link[stats]{model.frame}}: defaults to all rows. 
##'   This can be any valid indexing vector (see [.data.frame) for the rows of data 
##'   or if that is not supplied, a data frame made up of the variables used in formula.
##' @examples
##' # Omnibus contrast only of those effects included, namely condition and sex (no interaction)
##' omnibus_LSM_formula <- welchADF.test(y ~ condition + sex, data = womenStereotypeData)
##' # Works well with update.default() method
##' omnibus_interact_formula <- update(omnibus_LSM_formula, . ~ condition*sex)
##' summary(omnibus_LSM_formula)
##' summary(omnibus_interact_formula)
welchADF.test.formula <- function(formula, data, subset, ...){
  
  arguments = list(...)

  mf <- mc <- match.call()
  
  ## should do checks and tests as done in lme4::lFormula

  ## plain copy from lme4::lFormula start
  ## m stores the indices of those arguments
  m <- match(c("data", "subset"), names(mf), 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  fr.form <- subbars(formula) # substitute "|" by "+"
  environment(fr.form) <- environment(formula)

  mf$formula <- fr.form

  fr <- eval(mf, parent.frame())

  ## convert character vectors to factor (defensive)
  fr <- factorize(fr.form, fr, char.only=TRUE)

  ## store full, original formula
  attr(fr,"formula") <- formula
  
  ## plain copy from lme4::lFormula end
  terms <- attr(fr, "terms")
  
  ## find subject and within.s
  within.s <- subject <- NULL
  bars <- findbars(formula)
  if (length(bars) > 1) {
    stop("Please specify only one subject column.")
  } else if (length(bars) == 1) {
    subject <- toChar(bars[[1]][[3]])
    within.s <- toChar(bars[[1]][[2]])
  }
  
  ## response
  response <- toChar(formula[[attr(terms, "response") + 1]])
  
  if (length(response) > 1) { ## multivariate
    ## trick to make fr usuable
    fr <- cbind(fr[[1]], fr[-1])
    response <- response[response %in% colnames(fr)]
  }
  
  ## between.s
  fixedform <- formula
  RHSForm(fixedform) <- nobars(RHSForm(fixedform))
  candidates.between <- attr(terms(fixedform), "term.labels")
  between.s <- candidates.between[!(
    grepl(":", candidates.between, fixed = TRUE) |
    candidates.between %in% within.s
  )]
  
  # compose the complete list of the effects and interactions of the model, to be
  # applied a separate omnibus test if required. Those interactions not present here
  # will not be displayed in the result of the omnibus contrast
  
  effects.within = NULL
  if(length(bars) == 1){
    dummy.expression.within = as.character(as.expression(bars[[1]][[2]]))
    dummy.within.formula = formula(paste0("~", dummy.expression.within))
    effects.within = attr(terms(dummy.within.formula), "term.labels")
  }  
  only = unique(c(candidates.between, effects.within))

  arguments[["formula"]] = fr
  arguments[["response"]] = response
  arguments[["between.s"]] = between.s
  arguments[["within.s"]] = within.s
  arguments[["subject"]] = subject
  arguments[["only"]] = only
  
  wtestobj = do.call(welchADF.test.default, args = arguments)
  # Store the original call so that the user can later call update()
  # To avoid visibility problems, we replace the specialized function name by the generic
  mc[[1]] = quote(welchADF.test)
  wtestobj$call = mc
  wtestobj$model = fr

  wtestobj
}

##' @describeIn welchADF.test Specialized method that accepts a linear model object of class \code{\link[stats]{lm}}
##' @export
##' @method welchADF.test lm
##' @examples
##' 
##' # Fit a linear model using the built-in function stats::lm
##' lm.women <- lm(y ~ condition + sex, womenStereotypeData)
##' 
##' # Fit an Analysis of Variance model using the built-in function stats::aov
##' aov.women <- aov(lm.women)
##' 
##' # Now use the this object to apply a welchADF test to the same formula and data
##' omnibus_no_interact <- welchADF.test(lm.women, contrast = "omnibus")
##' omnibus_no_interactB <- welchADF.test(aov.women) # omnibus as well
##' 
##' # Integrates well with the update.default() method
##' omnibus_interact <- update(omnibus_no_interact, . ~ condition*sex)
##' summary(omnibus_no_interact)
##' summary(omnibus_interact)
welchADF.test.lm <- function(formula, ...) {
  ## get original call
  mc <- getCall(formula)
  
  ## should probably clean up call as to avoid warnings
  ## now add args from ...
  args <- list(...)
  for(arg in names(args)) {
    mc[[arg]] <- args[[arg]]
  }
  mc[[1L]] <- quote(welchADF.test)
  
  wtestobj <- eval(mc, parent.frame())
  
  # store the call so the user can later call update()
  # To avoid visibility problems, we replace the specialized function name by the generic
  wtestobj$call <- mc
  
  wtestobj
}

##' @describeIn welchADF.test Specialized method that accepts an Analysis of Variance Model of class \code{\link[stats]{aov}}
##' @method welchADF.test aov
##' @export
welchADF.test.aov <- welchADF.test.lm

##' @describeIn welchADF.test Specialized method that accepts a Linear Mixed-Effects Model of class \code{\link[lme4]{lmer}}
##' @method welchADF.test lmer
##' @export
welchADF.test.lmer <- welchADF.test.lm

##' @describeIn welchADF.test Default method that accepts a data.frame (see argument \code{formula}) and where 
##'   factors are passed as strings corresponding to column names
##' @export
##' @param response A string or vector of strings with the name(s) of the column(s) of \code{data} corresponding to the response variable(s). If 
##'		a vector of strings is provided, the responses are taken as a set of repeated measurements or dependent variables.
##' @param between.s Vector of strings with the columns that correspond to between-subjects factors, if any (defaults to NULL).
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
##'		If \code{contrast = "omnibus"} but \code{effect} is not specified, then an omnibus contrast will be done separately for each of the effects and also all possible interactions.
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
##' @references Villacorta, P.J. (2017). The welchADF Package for Robust Hypothesis Testing in Unbalanced Multivariate Mixed Models with Heteroscedastic and Non-normal Data. 
##' \emph{The R Journal}, 9:2, 309 - 328.
##' @references Website with the original SAS code and examples: \url{http://supp.apa.org/psycarticles/supplemental/met_13_2_110/met_13_2_110_supp.html}
##' @references Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach to hypothesis testing in independent and correlated groups designs. 
##' \emph{Psychophysiology}, 40, 586-596.
##' @references Lix, L. M., & Keselman, H. J. (1995). Approximate degrees of freedom tests: A unified perspective on testing for mean equality. 
##' \emph{Psychological Bulletin}, 117, 547-560. 
##' @references Carpenter, J., & Bithell, J. (2000). Bootstrap confidence intervals: When, which, what? A practical guide for medical statisticians. 
##' \emph{Statistics in Medicine}, 19, 1141-1164. 
##' @references Efron, B., & Tibshirani, R. (1986). Bootstrap methods for standard errors, confidence intervals, and other measures of statistical accuracy. 
##' \emph{Statistical Science}, 1, 54-75. 
##' @seealso \code{\link{p.adjust.methods}} 
##'   \code{\link{perceptionData}}
##'   \code{\link{adhdData}}
##'   \code{\link{adhdData2}}
##'   \code{\link{womenStereotypeData}}
##'   \code{\link{miceData}}
##' @examples
##' # Two-way factorial design using the default interface. See the vignette.
##' omnibus_LSM <- welchADF.test(womenStereotypeData, response = "y", between.s = 
##'   c("condition", "sex"), contrast = "omnibus")
##' # Method update() also works with the welchADF.test.default interface
##' omnibus_trimmed <- update(omnibus_LSM, effect = "condition", trimming = TRUE)
##' pairwise_LSM <- update(omnibus_LSM, contrast = "all.pairwise", effect = c("condition", "sex"))
##' pairwise_trimmed <- update(pairwise_LSM, trimming = TRUE) # just trimming
##' summary(omnibus_LSM)
##' pairwise_LSM
##' \dontrun{
##'   pairwise_trimmed_boot <- update(pairwise_trimmed, bootstrap = TRUE) # trimming and bootstrapping
##'   summary(pairwise_trimmed_boot)
##' }
welchADF.test.default <- function(formula, response, between.s = NULL, within.s = NULL, subject = NULL, contrast = c("omnibus", "all.pairwise"), 
                                     effect = NULL, correction = c("hochberg", "holm"), trimming = FALSE, per=0.2, bootstrap = FALSE, 
                                     numsim_b = 999, effect.size = FALSE, numsim_es = 999, scaling = TRUE, standardize.effsz = TRUE, alpha = 0.05, seed = 0, ...){
  
  data = formula
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
      levelslist.between.s = levelslist.between.s, levelslist.within.s = levelslist.within.s,
      only = list(...)[["only"]])
    
    omnibus.results$call = match.call() # this might replaced by any other non-default specialized function
                                        # if this function was called from welchADF.formula, welchADF.lm, lme, etc
    
    # replace the specialized function name by the generic to avoid visibility problems
    omnibus.results$call[[1]] = as.symbol("welchADF.test") 
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
    
                                         # save the call
    pairwise.results$call = match.call() # this might be replaced by any other non-default specialized function
                                         # if this function was called from welchADF.formula, welchADF.lm, lme, etc
    
    # replace the specialized function name by the generic to avoid visibility problems
    pairwise.results$call[[1]] = as.symbol("welchADF.test") 
    return(pairwise.results)
    
  }	# if contrast == "all.pairwise"
}
