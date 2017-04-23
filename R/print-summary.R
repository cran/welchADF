print.WelchTestObj<-function(x, ...){
  cat(format(x, ...), "\n")
}

format.WelchTestObj<-function(x,...){
	
  trimming = x$trimming
	bootstrap = x$bootstrap
	effsz = x$compute.effsz	

	#	****print results****
	cat("Welch-James Approximate DF Solution\n")
	if(!trimming){ cat("Least Squares Means & Variances\n") }
	else{ 
		cat("Trimmed Means & Winsorized Variances\n")
		cat("Percentage of Trimming:", x$trimming.per, "\n")
	}
	if(!bootstrap && !effsz){ cat("F Distribution Critical Value\n") }
	if(bootstrap && !effsz){ 
		cat("Bootstrap Critical Value for Single Test Statistic\n")
		cat("Number of Bootstrap Samples:", x$boot.samples, "\n")
	  cat("starting Seed:", x$seed, "\n")
	}
	cat("Contrast Matrix:\n") 	
	print(x$contrast.matrix)
	cat("Mean Vector: ")				
	print(x$mean.vector)
	cat("Sigma Matrix:\n")			
	print(x$sigma.matrix)
		
	if(!effsz){ 	# don't compute effect size
		cat("Significance Test Results:\n")
		if(!is.null(x$pval)){	# the object was created by function wjglm
			
		  results = NULL
		  if(x$compute.effsz){
		    # Effect size
		    results = t(c(x$welch.T, x$numeratorDF, x$denominatorDF, x$pval, x$effsz))
		    colnames(results) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval", "eff.size")
		  }
		  else{
		    # No effect size
		    results = t(c(x$welch.T, x$numeratorDF, x$denominatorDF, x$pval))
		    colnames(results) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval")
		  }
			
			rownames(results) = x$effect.name
			print(results)
		}
		else{	# the object was created by function bootcom
			results = rbind(x$welch.T, x$numeratorDF, x$denominatorDF)
			rownames(results) = c("WJ statistic", "Numerator DF", "Denominator DF")
			print(results)
			cat("Critical value:", x$crit.value, "\n")
		}
	}
	else{
		cat("Effect Size:", x$effsz, "\n")
		if(!x$scaling){ cat("No Scaling Factor\n") }
		else{ 
			cat("Scaling Factor is:", x$scaling.factor, "\n")			
		}
		if(x$standardize.effsz){ 
			cat("Standardizer Is Square Root of Average Variance\n") 
		}
		else{
			cat("No Standardizer\n")
		}
		cat("Number of Effect Size Bootstrap Samples:", x$effsz.samples, "\n")
	  cat("Starting Seed:", x$seed, "\n")
		cat("Confidence Level (%):", x$CI.level, "\n")
		cat("Confidence interval: [", x$CI[1], ",", x$CI[2],"]\n")			
		
	} 
}


summary.WelchTestObj<-function(object,...){  

	cat("Welch-James Approximate DF Test (")
	if(!object$trimming){ cat("Least squares means & variance)\n")  }	
	else{ cat(paste0("Trimmed means [", object$trimming.per, "% trimming] & Winsorized variances)\n")) }
	cat("Significance Test Results")
	if(!is.null(object$pval)){	# the object was created by function wjglm
		cat("\n")
		results = t(c(object$welch.T, object$numeratorDF, object$denominatorDF, object$pval))		
		colnames(results) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval")
		if(!is.null(object$adj.pval)){
			results = t(c(object$welch.T, object$numeratorDF, object$denominatorDF, object$pval, object$adj.pval))
			colnames(results) = c("WJ statistic", "Numerator DF", "Denominator DF", "pval", "adj.pval")
		}
		rownames(results) = object$effect.name
		print(results)
	}
	else{	# the object was created by function bootcom
		cat(" using a bootstrap critical value for FWER control\n")

		significance = object$welch.T > object$crit.value
		results = cbind(object$welch.T, object$numeratorDF, object$denominatorDF, significance)
		colnames(results) = c("WJ statistic", "Numerator DF", "Denominator DF", "significant?")

		if(!is.null(object$effect.name)){		
			rownames(results) = object$effect.name
		}
		print(results)
		cat("Bootstrap critical value:", object$crit.value, "\n")
	}
}

summary.WelchTestObjList<-function(object, ...){

	cat("Welch-James Approximate DF Test (")
	if(!object[[1]]$trimming){ cat("Least squares means & variances)\n")  }	
	else{ cat(paste0("Trimmed means [", 100*object[[1]]$trimming.per, "% trimming] & Winsorized variances)\n")) }
	type = attr(object, "type")
	if(type == "omnibus"){
		cat("Omnibus test(s) of effect and/or interactions")
		if(object[[1]]$bootstrap){
			cat(" with Bootstrap Critical Value for each Test Statistic")	
		}
		cat("\n")
	}
	else if(type == "all.pairwise"){
		effect = attr(object, "effect")
		if(length(effect) == 1){
			cat("Multiple pairwise comparisons of response marginal means with respect to",effect,"main effect")
		}
		else{
			temp = as.list(effect)				
			temp$sep = " x "
			temp = do.call(paste, temp)
			cat("Multiple tetrad interaction contrasts with respect to", temp, "interaction")
		}
		if(attr(object, "bootcom")){
			cat("\nusing Bootstrap Critical Value for FWER control")
		}
		else{
		  correction = attr(object,"correction")
		  correction= paste(toupper(substring(correction, 1, 1)), substring(correction, 2),
		                    sep = "", collapse = " ")
		  cat("\nusing",correction,"method for adjusting the p-values\n")
		}
		if(object[[1]]$compute.effsz){
		  cat("Effect size ")
		  if(object[[1]]$standardize.effsz){
		    cat("standardized via square root of the average of cell variances\n")
		  }
		  else{
		    cat("with no standardization\n")
		  }
		}
		cat("\n")
	}
	mat = attr(object, "summary")
	print(mat)
	bootcom = attr(object, "bootcom")
	if(!is.null(bootcom)){
	  if(bootcom == TRUE){
		  critv = attr(object, "critv")
		  cat("Bootstrap critical value:", critv, "\n")
	  }
	}
}
