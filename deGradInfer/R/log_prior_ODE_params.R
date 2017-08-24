# Calculate Log Prior for ODE Parameters

calculateLogParamPrior <- function(params,auxVars) {

	if(is.character(auxVars$logPrior)) {

	  if(auxVars$logPrior=="Mixed") {
	    ### Mixed prior from Campbell and Steele
	    return(sum(dnorm(params[1:2],0,0.4,log=TRUE),dchisq(params[3],2,log=TRUE)))
	  } else if(auxVars$logPrior=="Gamma") {
	  		# Gamma Prior
  			return(sum(dgamma(params, 4, 2, log=TRUE)))
 		} else if(auxVars$logPrior=='Uniform') {
 		    # Uniform prior
 		    return(0)
 		}

	} else if(is.function(auxVars$logPrior)) {
		return(sum(auxVars$logPrior(params)))
	}

   # Failure case if not a function or matching string
	 stop('logPrior must be one of "Uniform", "Mixed" or "Gamma", or a
	       user-specified function.')

}

