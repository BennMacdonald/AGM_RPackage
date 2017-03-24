# Sample ODE Mismatch parameter

sampleLambda <- function(lambda, gpFit, x, parameters, 
                    timePoints, auxVars, species, chain, chainTemp) {
  # Propose from prior
  proposal = proposeLambda(lambda[species], auxVars)
  new.lambda = lambda
  new.lambda[species] = proposal$lambda

  #oldLL = calculateLogLikelihoodSpecies(parameters, gpFit[[species]], x, lambda[species],
  #          timePoints, auxVars, 
  # 		     'MCMC', species, chain)
  oldLL = calculateLogLikelihood(parameters, gpFit, x, lambda, timePoints, auxVars, 
                                   'MCMC', chain, includeDet=T) 

  #det.old = determinant.matrix(oldLL$noiseA, logarithm=T)
  prior.old = dexp(lambda[species], 1, log=T)
  
  auxVars$lambdaChanged = species
  #newLL = calculateLogLikelihoodSpecies(parameters, gpFit[[species]], x, proposal$lambda, 
  #         timePoints, auxVars, 
  #  	     'MCMC', species, chain)
  newLL = calculateLogLikelihood(parameters, gpFit, x, new.lambda, timePoints, auxVars, 
                                   'MCMC', chain, includeDet=T)  

  #det.new = determinant.matrix(newLL$noiseA, logarithm=T)
  prior.new = dexp(proposal$lambda, 1, log=T)

  ratio = exp(chainTemp * (newLL$LL - 0.5*newLL$log.det - 
    (oldLL$LL - 0.5*oldLL$log.det)) + prior.new - prior.old)

  #ratio = exp(chainTemp*(newLL$logLikelihood - 0.5*det.new$modulus - 
  #          (oldLL$logLikelihood - 0.5*det.old$modulus)) + prior.new - prior.old)

  #lL.old = oldLL$logLikelihood - 0.5*det.old$modulus
  lL.old = oldLL$LL - 0.5 * oldLL$log.det + oldLL$gpXPrior

  if(min(1, ratio) > runif(1)) {
    sampled.lambda = proposal$lambda
    accept = 1 
    invNoiseA.rec[[chain]][,species] <<- invNoiseA.temp
    #lL.new = newLL$logLikelihood - 0.5*det.new$modulus 
    lL.new = newLL$LL - 0.5 * newLL$log.det + newLL$gpXPrior
    auxVars$lambdaChanged = 0
    newLL.temp = calculateLogLikelihood(parameters, gpFit, x, new.lambda, timePoints, auxVars, 
                                   'MCMC', chain, includeDet=T)    
    lL = newLL.temp$LL - 0.5*newLL.temp$log.det + newLL.temp$gpXPrior
  } else {
    sampled.lambda = lambda[species]
    accept = 0
    lL.new = lL.old
    lL = lL.old
  }


  return(list(lambda=sampled.lambda, accept=accept, lL.old = lL.old, lL.new = lL.new, lL=lL))
}

proposeLambda <- function(lambda, auxVars) {
  #lambda.new = rgamma(1, 0.5, 3)
  #lambda.new = uniformProposal(lambda, 0.005, 0.1) 
  lambda.new = rexp(1,10)
  #rexp(1, auxVars$lambdaPrior)

  return(list(lambda=lambda.new))
}
