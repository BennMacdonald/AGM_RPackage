# Sample from proposal distribution for MCMC (currently uniform)
proposeParamsMCMC <- function(oldParams, inferredParams, width,
                              covEstimate, explicit=F) {
  
  changed = matrix(0, length(oldParams), 1)
  newParams = oldParams

  if(F && explicit) {
    newParams[inferredParams] = newParams[inferredParams] + 
      width[inferredParams]*mvrnorm(1, rep(0, length(inferredParams)),
                            covEstimate)
    changed[inferredParams] = 1
  } else {
    for(i in 1:length(inferredParams)) {
      choice = resample(inferredParams, 1) 
     
      changed[choice] = 1  
      param = choice 

      # Gaussian random walk proposal (like in Calderhead)
      newParams[param] = newParams[param] + width[param]*rnorm(1)
    
      # Log Gaussian random walk proposal
      # newParams[param] = exp(log(newParams[param]) + width[param]*rnorm(1))
    }
  }

  return(list(params = newParams, changed=changed, oldProb=1, newProb=1))
}

# Draw a new parameter uniformly based on the previous parameter value
uniformProposal <- function(oldParam, range, limit) {
  if(runif(1) < 0.5) {
    newParam = oldParam + runif(1)*range
  } else {
    newParam = oldParam - runif(1)*range
  } 

  while(newParam > limit) {
    newParam = newParam - limit
  } 

  while(newParam < 0) {
    newParam = limit + newParam
  }
  return(newParam)
}

# Log likelihood of one species
calculateLogLikelihoodMCMC <- function(params, gpFit, X, lambda, timePoints, 
  auxVars, species, chain) {
  # Calculate summary statistics
  gpSummary = likelihoodUtil(params, X, lambda, timePoints, auxVars, gpFit, species, chain)
  
  noiseA = gpSummary$noiseA
  gradDiff = gpSummary$gradDiff
  invK = gpSummary$invK
  m = gpSummary$m

  if(gpSummary$error) {
    logLikelihood = -1e6
    gpXPrior = 1e6
  } else { 
    # Prior for latent variables
    gpXPrior = t.default(X[,species,drop=F]) %*% invK %*% X[,species,drop=F]
  
    # Main log likelihood term
    
    if(auxVars$Kchanged == species ||
       auxVars$lambdaChanged == species) {
      invNoiseA = NULL
      try(invNoiseA <- solve.default(noiseA), silent=T)
      
      if(!is.null(invNoiseA)) 
        invNoiseA.temp <<- c(invNoiseA)
    } else {
      invNoiseA = matrix(invNoiseA.rec[[chain]][, species], 
                         length(timePoints), 
                         length(timePoints))
    }
    
    if(!is.null(invNoiseA)) {
      #browser()
      prodXdot = t.default(gradDiff) %*% invNoiseA %*% gradDiff
      logLikelihood = - 0.5 * (prodXdot)
    } else {
      logLikelihood = -1e6
      gpXPrior = 1e6
      gpSummary$error = T
    }
  }

  return(list(gpXPrior = -0.5 * gpXPrior, logLikelihood = logLikelihood, m = m, 
              gradDiff=gradDiff, invK=invK, A=gpSummary$A, K=gpSummary$K, 
              deriv.m = gpSummary$deriv.m, 
              noiseA=noiseA, error=gpSummary$error))
}
