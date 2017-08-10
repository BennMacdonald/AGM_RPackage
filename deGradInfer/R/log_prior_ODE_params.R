# Calculate Log Prior for ODE Parameters

calculateLogParamPrior <- function(params,auxVars) {
  
	if(!is.null(auxVars$defaultLogParamPrior)){	
		if(auxVars$defaultLogParamPrior=="Uniform") {
			# Uniform prior
			return(0)
		}


		if(auxVars$defaultLogParamPrior=="Gamma") {
	  		# Gamma Prior
  			return(sum(dgamma(params, 4, 2, log=TRUE)))
			#sum(dgamma(params, 1.2, 0.75, log=TRUE))
 		}


		if(auxVars$defaultLogParamPrior=="Mixed") {
			### Mixed prior from Campbell and Steele
			return(sum(dnorm(params[1:2],0,0.4,log=TRUE),dchisq(params[3],2,log=TRUE)))
		}
	} else{

		return(sum(userLogPrior(params)))

	} # end else statement for if(!is.null(auxVars$defaultLogParamPrior)

}

