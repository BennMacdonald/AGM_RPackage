
# Main MCMC function: Runs the MCMC for the specified number of iterations and returns the 
# sampled parameter values 
doMCMC <- function(timePoints, data, auxVars, options) {

  # Initialisations
  iterationNum = options$iterations 
  chainNum = options$chainNum 
  paramNum = length(options$paramsInit)
  y = data
  swappedChains <- 0

  init = setupChains(timePoints, data, auxVars, options)
  auxVars = init$auxVars; tuning = init$tuning;
  x = init$x; gpFit = init$gpFit; parameters = init$parameters; 
  xRec = init$xRec; gpRec=init$gpRec; lambda = init$lambda; sigma=init$sigma
  speciesNum = init$speciesNum; paramsTempRec = init$paramsTempRec
  options$proposalGPTuning = init$proposalGPTuning; gradientMismatchParameterRec=init$lambda

  # When to monitor acceptance rates
  monitorRate = ifelse('monitorRate' %in% names(options), options$monitorRate, 
                       1000)
  # When to record current sample
  recordRate = ifelse('recordRate' %in% names(options), options$recordRate,
                      100)
  # When to save all samples
  saveRate = ifelse('saveRate' %in% names(options), options$saveRate, 1e4)

  #temp.exponent = 5
  #temperatures = seq(0, 1, length.out=chainNum)^temp.exponent

	powers <- 1:chainNum
	temp.exponent <- powers[options$temps]
	temperatures <- c()
	for(CHAIN in 1:chainNum)
	{
		temperatures[CHAIN] <- (CHAIN/chainNum)^temp.exponent
	}

  lLchains = matrix(-1e6, chainNum, 1); lLRec = matrix(0, 0, 1)
  dataLLchains = matrix(0, chainNum, speciesNum)
  dataLLchains[,auxVars$observedSpeciesList] = -1e6
  paramsRec = matrix(0, 0, paramNum)
  noiseRec = matrix(0, 0, speciesNum)
  lLAllChains <- matrix(0,0,options$chainNum)

  lastMove = ''
  pick = 1
  for(i in 1:iterationNum) {

    if(options$showProgress && i %% recordRate == 0) {
      cat((i/iterationNum)*100, '% Top Chain:', parameters[chainNum,options$inferredParams], "\n")
      cat('Noise', sigma[chainNum,], '\n') 
      cat('Lambda', lambda[chainNum,], '\n')
    }
    
    
   ### Specify lambda values if tempering the mismatch parameter
	if (auxVars$Mismatch$Tempering) { # For users who want to specify their own values
		if (!is.null(auxVars$Mismatch$lambdaValues)){
			lambda <- auxVars$Mismatch$lambdaValues
		}
		if (is.null(auxVars$Mismatch$lambdaValues)){ # Using the default LB2 or LB10 values
			lambda <- options$lambda
		}
	}


	# Sample Lambda
	if (!auxVars$Mismatch$Tempering){
       if(!options$explicit && runif(1) > 0.99) {
         pick = resample(auxVars$speciesList, 1)
         lambda.sampling = sampleLambda(lambda[chain,], gpFit[[chain]], x[[chain]],
                                        parameters[chain,], timePoints, auxVars, pick,
                                        chain, chainTemp)
         lambda[chain,pick] = lambda.sampling$lambda
         lLchains[chain] = lambda.sampling$lL.new
        
         cat(lambda.sampling$lL.old, lambda.sampling$accept, lambda.sampling$lL.new, '\n')
        
         if(lambda.sampling$accept) lastMove = 'LambdaAccept'
         else lastMove = 'LambdaReject'
        
         tuning$proposeLambdaTemp[chain,pick] = tuning$proposeLambdaTemp[chain,pick] + 1
         tuning$acceptLambdaTemp[chain,pick] = tuning$acceptLambdaTemp[chain,pick] + lambda.sampling$accept
        
         if(lLchains[chain] > 1e4) browser()
       }
	}


    # Chain-specific Inference
    for(j in 1:chainNum) {

      chain = resample(1:chainNum, 1)
      chainParams = parameters[chain,]
      chainTemp = temperatures[chain]
      u = runif(1,0,1)

      # GP Parameter Inference, coupled with latent variable inference
      if(!options$explicit && i > options$burnin && u > 0.98) {
        pick = resample(auxVars$speciesList, 1)
 
        
        gp.sampling = sampleGPX(gpFit[[chain]], sigma[chain, pick], x[[chain]], y[,pick],
                              lambda[chain,], parameters[chain,], timePoints, auxVars, pick,
                              options$proposalGPTuning[[chain]][,pick], chain, chainTemp, options)

        gpFit[[chain]] = gp.sampling$gp
      
        x[[chain]] = gp.sampling$x
     
        lLchains[chain] = gp.sampling$lL
        
        if(lLchains[chain] > 1e4) browser()
        
        dataLLchains[chain, pick] = gp.sampling$data.LL
        auxVars = gp.sampling$auxVars
        if(chain == chainNum) {
          #cat(lLchains[chain], sum(gp.sampling$accept), '\n')
          if(any(gp.sampling$accept==T)) lastMove = 'GPAccept'
          else lastMove = 'GPReject'
        }

        tuning$acceptGPTemp[[chain]][, pick] = 
          tuning$acceptGPTemp[[chain]][, pick] + gp.sampling$accept 
        tuning$proposeGPTemp[[chain]][, pick] = 
          tuning$proposeGPTemp[[chain]][, pick] + gp.sampling$changed
      } else if(i > options$noiseBurnin && u > 0.9) {
        # Noise Inference
        for(species in auxVars$observedSpeciesList) {
          pick = resample(auxVars$observedSpeciesList, 1) 
          
          noise.sampling = sampleNoise(sigma[chain, pick], x[[chain]][,pick], y[,pick],
                               options, chain, pick, temperatures,auxVars)
          sigma[chain, pick] = noise.sampling$sigma
          dataLLchains[chain, pick] = noise.sampling$data.LL
          tuning$proposeNoiseTemp[chain, species] = 
            tuning$proposeNoiseTemp[chain, species] + 1
          tuning$acceptNoiseTemp[chain, species]  = 
            tuning$acceptNoiseTemp[chain, species] + 
                                      noise.sampling$accept  
          if(chain == chainNum) lastMove = 'Noise'
        }
      # X Inference 
      } else if(!options$explicit && i > options$burnin && u > 0.75) { 
          pick = resample(auxVars$speciesList, 1)
          
          tuning$proposeXTemp[chain] = tuning$proposeXTemp[chain] + 1
          X.sampling = sampleX(chainParams, chain, temperatures,
                             gpFit[[chain]], x[[chain]], y[,pick], lambda[chain,], 
                             sigma[chain, pick], timePoints, 
                             auxVars, options, pick)

          x[[chain]][,pick] = X.sampling$x
          if(chain == chainNum) lastMove = 'X'
          lLchains[chain] = X.sampling$lL 
          dataLLchains[chain, pick] = X.sampling$data.LL

          tuning$acceptXTemp[chain, pick] = tuning$acceptXTemp[chain,pick] + 
            X.sampling$accept
          auxVars = X.sampling$auxVars
          
          if(lLchains[chain] > 1e4) browser()
      # Parameter Inference
      } else {
        sampling = sampleParams(chainParams, gpFit[[chain]], x[[chain]], y, lambda[chain,], sigma[chain,], 
                                  timePoints, 
                                  chainTemp,
                                  auxVars, options, chain)

        parameters[chain,] = sampling$parameters
        #if(chain == chainNum && (lLchains[chain] - sampling$lL) > 20) browser()
        lLchains[chain] = sampling$lL
        auxVars = sampling$auxVars
          
          if(chain == chainNum) lastMove = 'params'
        tuning$acceptTemp[chain,] = tuning$acceptTemp[chain,] + sampling$accepted
        tuning$proposeTemp[chain,] = tuning$proposeTemp[chain,] + sampling$proposed
        paramsTempRec[[chain]] = rbind(paramsTempRec[[chain]], 
          parameters[chain, options$inferredParams])
 
        chainParams = parameters[chain,]
         
        if(options$explicit) x[[chain]] = sampling$x
        
        if(lLchains[chain] > 1e4) browser()
      }        

    }
	if (!auxVars$Mismatch$Tempering){
        gradientMismatchParameterRec <- rbind(gradientMismatchParameterRec,lambda) # Record all lambda values
	}

    # Exchange Chains
    if(T) {
    for(j in 1:chainNum) {
      chain1 = resample(1:(chainNum-1), 1)
      chain2 = chain1 + 1
      
      if(chain2 == chainNum) lastMove = 'exchange'       
      exchange = exchangeChains(chain1, chain2, parameters, gpFit, x, lambda, 
                                timePoints, auxVars, temperatures, 
                                lLchains)
                                #lLchains+apply(dataLLchains, 1, sum))

      if(exchange$accepted) {
        params.temp = parameters[chain1,]
        x.temp = x[[chain1]]
        sigma.temp = sigma[chain1,]
        lL.temp = lLchains[chain1]
        data.ll.temp = dataLLchains[chain1,]
        gpFit.temp = gpFit[[chain1]]
        lambda.temp = lambda[chain1,]

        parameters[chain1,] = parameters[chain2,]    
        x[[chain1]] = x[[chain2]]
        lLchains[chain1] = lLchains[chain2]
        sigma[chain1,] = sigma[chain2,]
        dataLLchains[chain1,] = dataLLchains[chain2,]
        gpFit[[chain1]] = gpFit[[chain2]]
        lambda[chain1,] = lambda[chain2,]

        parameters[chain2,] = params.temp
        x[[chain2]] = x.temp
        lLchains[chain2] = lL.temp 
        sigma[chain2,] = sigma.temp
        dataLLchains[chain2,] = data.ll.temp
        gpFit[[chain2]] = gpFit.temp
        lambda[chain2,] = lambda.temp
        
        A.rec.temp = A.rec[[chain1]]
        deriv.m.rec.temp = deriv.m.rec[[chain1]]
        invK.rec.temp = invK.rec[[chain1]]
        K.rec.temp = K.rec[[chain1]]
        invNoiseA.rec.temp = invNoiseA.rec[[chain1]]

        A.rec[[chain1]] <<- A.rec[[chain2]]
        deriv.m.rec[[chain1]] <<- deriv.m.rec[[chain2]]
        invK.rec[[chain1]] <<- invK.rec[[chain2]]
        K.rec[[chain1]] <<- K.rec[[chain2]]
        invNoiseA.rec[[chain1]] <<- invNoiseA.rec[[chain2]]

        A.rec[[chain2]] <<- A.rec.temp
        deriv.m.rec[[chain2]] <<- deriv.m.rec.temp
        invK.rec[[chain2]] <<- invK.rec.temp
        K.rec[[chain2]] <<- K.rec.temp
        invNoiseA.rec[[chain2]] <<- invNoiseA.rec.temp

	  swappedChains <- swappedChains + 1

      } 
 
      tuning$proposeExchangeTemp[chain2] = tuning$proposeExchangeTemp[chain2] + 1
      tuning$acceptExchangeTemp[chain2] = tuning$acceptExchangeTemp[chain2] + 
                                            exchange$accepted

    }}
    
    if(i %% saveRate == 0) {
		if (!auxVars$Mismatch$Tempering){
      params = list(parameters=parameters, tuning=tuning, 
              paramsRec=paramsRec, lLRec=lLRec, xRec=xRec, gpRec=gpRec, timePoints=timePoints,
              noiseRec=noiseRec, swappedChains=swappedChains, 
		  chainNums=options$chainNum, maxIterations=options$iterations,
		  gradientMismatchParameterRec=gradientMismatchParameterRec,
		  lLAllChains=lLAllChains)
		}

		if (auxVars$Mismatch$Tempering){
      params = list(parameters=parameters, tuning=tuning, 
              paramsRec=paramsRec, lLRec=lLRec, xRec=xRec, gpRec=gpRec, timePoints=timePoints,
              noiseRec=noiseRec, swappedChains=swappedChains, 
		  chainNums=options$chainNum, maxIterations=options$iterations,
		  lLAllChains=lLAllChains)
		}

      save(params, file=options$saveFile)
    }
    
    if(i %% recordRate == 0) {
      paramsRec = rbind(paramsRec, parameters[chainNum,])
      noiseRec = rbind(noiseRec, sigma[chainNum,]) 

      for(species in auxVars$speciesList) {
        xRec[[species]] = rbind(xRec[[species]], t(x[[chainNum]][,species]))
        gpRec[[species]] = rbind(gpRec[[species]], 
                                 gpFit[[chainNum]][1:dim(gpRec[[species]])[2],species])
      }
      
	### Record Log Likelihood and Log Likelihood for all temperature chains
  lLRec = rbind(lLRec, lLchains[chainNum] + sum(dataLLchains[chainNum,]))
  
	
  # Vectorized version of what was previously a loop
  # Also fixed bug where variable i was used to index both chains and 
  # iterations
  lLAllStore = c(lLchains + rowSums(dataLLchains))
  
  lLAllChains <- rbind(lLAllChains,lLAllStore)
  
  if(options$showPlot) { 
    if(length(auxVars$speciesList) <= 6)
      par(mfrow=c(1+ceiling(length(auxVars$speciesList)/2),2))
    else par(mfrow=c(1,2))
    
    plot(lLRec, main=paste('Iterations:', i)) 
    boxplot(paramsRec[,options$inferredParams])    
    
    if(length(auxVars$speciesList) <= 6) {
      for(species in auxVars$speciesList) {
        y.max = max(c(x[[chainNum]][,species], y[,species]))
        y.min = min(c(x[[chainNum]][,species], y[,species]))
        plot(timePoints, x[[chainNum]][,species], ylim=c(y.min,y.max))
        points(timePoints, auxVars$y.true[,species], type='l')
        points(timePoints, y[,species], pch=2)
        text(timePoints[1], 0.05, 
             paste(round(gpFit[[chainNum]][,species], digits=2), collapse=' '), adj=c(0,0))
      }
    }
    
    Sys.sleep(0.005)
  } 
  
    }

    # Adjust proposal to get better acceptance rates
    if(i %% monitorRate == 0) {
      acceptRatio = tuning$acceptTemp[chainNum,]/tuning$proposeTemp[chainNum,]
      acceptNoiseRatio = tuning$acceptNoiseTemp[chainNum,]/tuning$proposeNoiseTemp[chainNum,]
      #acceptGPRatio[tuning$proposeGPTemp[[chainNum]] == 1] = 0.25
      #print(acceptRatio) 
      #print(options$proposalTuning[chainNum,])
      
      #print(acceptNoiseRatio) 
      #print(options$proposalNoiseTuning[chainNum,])

      #print(lLchains + apply(dataLLchains, 1, sum))

      new.tuning = adjustProposal(tuning, options)
      options$proposalTuning = new.tuning$tuningParams

      if(i > options$burnin) {
        options$proposalGPTuning = new.tuning$tuningGP
        options$proposalXTuning = new.tuning$tuningX
        options$proposalNoiseTuning = new.tuning$tuningNoise
      } 
      
      for(chain in 1:chainNum) {
        auxVars$paramsCovEstimate[[chain]] = cov(paramsTempRec[[chain]])
        paramsTempRec[[chain]] = matrix(0, 0, length(options$inferredParams))
      }


      temperatures = adjustExchangeProposal(tuning,
                                            temperatures, temp.exponent)
      tuning = resetAndUpdateRates(tuning, chainNum, paramNum, speciesNum) 

    }
 
  }
  
		if (!auxVars$Mismatch$Tempering){
  return(list(parameters=parameters, tuning=tuning, 
              paramsRec=paramsRec, lLRec=lLRec, xRec=xRec, gpRec=gpRec, timePoints=timePoints,
              noiseRec=noiseRec, swappedChains=swappedChains, 
		  chainNums=options$chainNum, maxIterations=options$iterations,
		  gradientMismatchParameterRec=gradientMismatchParameterRec,
		  lLAllChains=lLAllChains))
		}

		if (auxVars$Mismatch$Tempering){
  return(list(parameters=parameters, tuning=tuning, 
              paramsRec=paramsRec, lLRec=lLRec, xRec=xRec, gpRec=gpRec, timePoints=timePoints,
              noiseRec=noiseRec, swappedChains=swappedChains, 
		  chainNums=options$chainNum, maxIterations=options$iterations,
		  lLAllChains=lLAllChains))
		}

}


# Sample X from Y
sampleX <- function(parameters, chain, temperatures, gpFit, x, y, lambda, sigma, 
                    timePoints, auxVars, options, species, flag=F) {

  samplingStrategy = options$samplingStrategy
  
  if(samplingStrategy == 'mixed') {
    samplingStrategy = ifelse(runif(1) < 0.3, 'MCMC', 'HMC')
  }
   
  # Hamiltonian Monte Carlo 
  if(samplingStrategy == 'HMC') {
    proposal = proposeXHMC(x, y, timePoints, parameters, 
                 lambda[species], sigma, auxVars, gpFit[[species]], 
                 options$proposalXTuning, chain, species, temperatures[chain])
    proposal$invM = 1 
  # Hamiltonian Monte Carlo with mass matrix
  } else if(samplingStrategy == 'HMCM') {
    proposal = proposeXHMCM(x, y, timePoints, parameters, lambda[species], sigma, auxVars, gpFit[[species]], 
                           options$proposalXTuning, chain, species, temperatures[chain])
  # Markov Chain Monte Carlo
  } else {
    proposal = proposeX(x, species, options$proposalXTuning, chain)
    proposal$invM = 1
    proposal$p = 0
  }
  
  proposal.X = proposal$x
  
  oldLL = calculateLogLikelihood(parameters, gpFit, x, lambda,
            timePoints, auxVars, 
  		      samplingStrategy, chain, proposal$p, proposal$invM, includeDet=T)
  
  newLL = calculateLogLikelihood(parameters, gpFit, proposal.X, lambda,
           timePoints, auxVars, 
    	     samplingStrategy, chain, proposal$p, proposal$invM, includeDet=T)

  auxVars = newLL$auxVars  
  
  oldPrior = sum(dnorm(y, x[,species], sigma, log=T))
  newPrior = sum(dnorm(y, proposal.X[,species], sigma, log=T))

  ratio = exp((newPrior - oldPrior) + temperatures[chain] * 
      ((newLL$gpXPrior - oldLL$gpXPrior) + 
      (newLL$LL - oldLL$LL)))

  if(is.na(ratio)) {
    browser() 
  } else if((options$allowNeg || all(proposal.X[,species] > 0)) && 
             min(1, ratio) > runif(1)) {
    sampled.X = proposal.X[,species]
    accept = 1 
    lL = newLL$LL + newLL$gpXPrior - 0.5 * newLL$log.det
    data.LL = newPrior 
  } else {
    sampled.X = x[,species]
    accept = 0
    lL = oldLL$LL + oldLL$gpXPrior - 0.5 * oldLL$log.det
    data.LL = oldPrior
  }

  return(list(x=sampled.X, accept=accept, lL=lL, auxVars=auxVars, data.LL=data.LL))
}


# Sample GP and X
sampleGPX <- function(gpFit, sigma, x, y, lambda, parameters, 
                    timePoints, auxVars, species, proposal.width, chain, chainTemp, options) {
  
  # Propose new GP parameters
  proposal = proposeGP(gpFit, proposal.width, species)
  accept = matrix(0, length(proposal$changed), 1)

  # Calculate covariance matrix and derivatives
  gpCovs = getGPCovs(proposal$gp[,species], auxVars)

  # Estimate nu 
  lower = t.default(chol(matrix(K.rec[[chain]][,species], length(timePoints), length(timePoints))))
  nu = solve(lower) %*% x[,species] 

  # Check if Cholesky decomposition possible for the new parameters
  upper.new = NULL 
  try(upper.new <- chol(gpCovs$K), silent=T)
  
  oldLL = calculateLogLikelihood(parameters, gpFit, x, lambda,
           timePoints, auxVars, 
  		     'MCMC', chain, includeDet=T)
  
  if(species %in% auxVars$observedSpeciesList) 
    x.ll = sum(dnorm(x[,species], y, sigma, log=T)) 
  else
    x.ll = 0

  if(!is.null(upper.new)) {

    gp.prior.old = calculateLogGPPrior(gpFit[,species], auxVars$covtype)
    gp.prior.new = calculateLogGPPrior(proposal$gp[,species], auxVars$covtype)

    # Update latent variables
    lower.new = t(upper.new)
    x.new = x
    x.new[,species] = lower.new %*% nu

    # Check fit to data
    if(species %in% auxVars$observedSpeciesList)
      x.ll.new = sum(dnorm(x.new[,species], y, sigma, log=T))
    else
      x.ll.new = 0

    auxVars$Kchanged = species
    newLL = calculateLogLikelihood(parameters, proposal$gp, x.new, lambda,
           timePoints, auxVars, 
    	     'MCMC', chain, includeDet=T)

    ratio = exp(gp.prior.new - gp.prior.old + x.ll.new - x.ll + 
                chainTemp*(newLL$LL - oldLL$LL +
                newLL$gpXPrior - oldLL$gpXPrior - 0.5*(newLL$log.det - oldLL$log.det)))
  } else {
    # New parameters do not produce positive definite K
    ratio = 0
  }
  
  if(!is.nan(ratio) &&  min(1, ratio) > runif(1) && 
    (options$allowNeg || all(x.new[,species] > 0))) {
    sampled.gp = proposal$gp; accept = proposal$changed 
    lL = newLL$LL - 0.5 * newLL$log.det + newLL$gpXPrior
    data.LL = x.ll.new
    auxVars = newLL$auxVars; x = x.new
    A.rec[[chain]][,species] <<- A.temp
    deriv.m.rec[[chain]][,species] <<- deriv.m.temp
    invK.rec[[chain]][,species] <<- invK.temp
    K.rec[[chain]][,species] <<- K.temp
    invNoiseA.rec[[chain]][,species] <<- invNoiseA.temp
    auxVars$Kchanged = 0 
  } else {
    sampled.gp = gpFit; 
    lL = oldLL$LL - 0.5 * oldLL$log.det + oldLL$gpXPrior
    auxVars = oldLL$auxVars; x = x
    data.LL = x.ll
  }

  return(list(gp=sampled.gp, accept=accept, changed=proposal$changed, lL=lL, auxVars=auxVars, 
              x=x, data.LL=data.LL))

}

# Calculate the prior probabilities for the the hyperparameters of the GP
calculateLogGPPrior <- function(params, covtype) {
  log.gp.prior = 0 

  if(covtype == 'sigmoidVar') {
    #log.gp.prior = log.gp.prior + dgamma(params[1], 1, 1, log=T)
    log.gp.prior = log.gp.prior + dgamma(params[2], 5, 4, log=T)
  } 

  if(covtype == 'sigmoidVar') {
    log.gp.prior = log.gp.prior + dgamma(params[3], 5, 4, log=T)
  }
  
  return(log.gp.prior)
}

# Sample X from Y
sampleGP <- function(gpFit, sigma, x, y, lambda, parameters, 
                    timePoints, auxVars, species) {
  proposal = proposeGP(gpFit[[species]])

  oldLL = calculateLogLikelihoodSpecies(parameters, gpFit[[species]], x, lambda,
           timePoints, auxVars, 
  		     'MCMC', species)
   
  det.old = -0.5*(determinant(oldLL$noiseA, logarithm=T)$modulus + 
                  determinant(oldLL$K, logarithm=T)$modulus)
 
  auxVars$Kchanged[species] = T
  newLL = calculateLogLikelihoodSpecies(parameters, proposal$gp, x, lambda,
           timePoints, auxVars, 
    	     'MCMC', species)
  det.new = -0.5*(determinant(newLL$noiseA, logarithm=T)$modulus + 
                  determinant(newLL$K, logarithm=T)$modulus) 

  ratio = exp((newLL$gpXPrior - oldLL$gpXPrior) + 
      (newLL$logLikelihood - oldLL$logLikelihood) + (det.new - det.old))

  if(min(1, ratio) > runif(1)) {
    sampled.gp = proposal$gp
    accept = 1 
    lL = newLL$LL + newLL$gpXPrior
    auxVars = newLL$auxVars  
  } else {
    sampled.gp = gpFit[[species]]
    accept = 0
    lL = oldLL$LL + oldLL$gpXPrior
    auxVars = oldLL$auxVars  
  }
  
  return(list(gp=sampled.gp, accept=accept, lL=lL, auxVars=auxVars))
}

# Normal proposal changing one of the GP parameters
proposeGP <- function(gp.orig, proposal.width, species) {

  gp.new = gp.orig
  changed = matrix(0, length(proposal.width), 1)

  #for(i in 1:length(gp.new[[species]]$params)) {
    choice = resample(1:length(gp.new[,species]), 1) 
    gp.new[choice,species] = gp.new[choice, species] + 
      rnorm(1, 0, proposal.width[choice])
  
    gp.new[choice, species] = abs(gp.new[choice, species])

    changed[choice] = 1
  #}

  return(list(gp=gp.new, changed=changed))
}

# Uniform proposal changing one or more of the X values
proposeX <- function(x.orig, species, proposal.width, chain) {
  x.new = x.orig
  
  changed = matrix(0, length(x.orig[,species]), 1)
  
  for(i in 1:length(x.orig)) {
    choice = sample(1:length(x.orig[,species]), 1) 

    x.new[choice, species] = x.orig[choice, species] + 
                               proposal.width[chain, species]*rnorm(1)

    if(abs(x.new[choice, species]) > 100) {
       x.new[choice,species] = x.orig[choice, species]
    } 

    changed[choice] = 1
  }
  
  return(list(x.new=x.new, changed=changed))
}

# Exchange move for two chains      
exchangeChains <- function(chain1, chain2, parameters, gpFit, x, lambda, timePoints, 
                           auxVars, temperatures, lL) {
        
  chain1Params = parameters[chain1,]
  chain2Params = parameters[chain2,]
  chain1X = x[[chain1]]
  chain2X = x[[chain2]]
  chain1Lambda = lambda[chain1,]
  chain2Lambda = lambda[chain2,]
  chain1Temp = temperatures[chain1]
  chain2Temp = temperatures[chain2]
  
  #chain2LL = calculateLogLikelihood(chain2Params, gpFit, chain2X, chain2Lambda, 
  #                                  timePoints, auxVars, 'MCMC', includeDet=T)
  
 
  #oldTempLL = chain1Temp*(chain1LL$gpXPrior + chain1LL$LL) + 
  #            chain2Temp*(chain2LL$gpXPrior + chain2LL$LL)

  #newTempLL = chain2Temp*(chain1LL$gpXPrior + chain1LL$LL) + 
  #            chain1Temp*(chain2LL$gpXPrior + chain2LL$LL)
 
  oldTempLL = chain1Temp*lL[chain1] + chain2Temp*lL[chain2]
  newTempLL = chain2Temp*lL[chain1] + chain1Temp*lL[chain2]
 
  ratio = exp(newTempLL - oldTempLL)
  
  exchange = list()
  
  if(min(1, ratio) > runif(1)) {  
    exchange$chain1 = chain2
    exchange$chain2 = chain1
    #exchange$lLNewChain1 = chain2Temp*(chain1LL$gpXPrior + chain1LL$LL - 0.5*chain1LL$log.det)
    #exchange$lLNewChain2 = chain1Temp*(chain2LL$gpXPrior + chain2LL$LL - 0.5*chain2LL$log.det)
    exchange$accepted = 1
  } else {
    exchange$chain1 = chain1
    exchange$chain2 = chain2    
    #exchange$lLNewChain1 = chain1Temp*(chain1LL$gpXPrior + chain1LL$LL - 0.5*chain1LL$log.det)
    #exchange$lLNewChain2 = chain2Temp*(chain2LL$gpXPrior + chain2LL$LL - 0.5*chain2LL$log.det)
    exchange$accepted = 0
  }
  
  return(exchange)
}

# Sample ODE parameters
sampleParams <- function(oldParams, gpFit, data, y, lambda, sigma, timePoints, temperature, 
                         auxVars, options, chain) {
  # Decide sampling Strategy
  samplingStrategy = options$samplingStrategy

  if(samplingStrategy == 'mixed') {
    samplingStrategy = ifelse(runif(1) < 0.95, 'HMC', 'MCMC')   
  }
 
  if(samplingStrategy == 'HMCM') samplingStrategy = 'MCMC'

  # Propose new set of parameters 
  proposal = proposeParams(data, timePoints, oldParams, lambda, auxVars, gpFit, 
                           samplingStrategy, options, chain, temperature)

  params = oldParams
  accept = matrix(0, length(oldParams), 1)
  proposed = matrix(0, length(oldParams), 1)

  # Log likelihood of old parameters 
  oldPrior = calculateLogParamPrior(oldParams[options$inferredParams],auxVars)
  
  if(options$explicit) {
    oldLL = calculateLogLikelihoodExplicit(oldParams, y, timePoints, sigma, auxVars)
  } else {
    oldLL = calculateLogLikelihood(oldParams, gpFit, data, lambda, timePoints, auxVars, 
			     samplingStrategy, chain,
			     proposal$old.p, includeDet=T)
  }

  x = oldLL$x 
  
  
  gpXPrior = oldLL$gpXPrior; auxVars = oldLL$auxVars

  lL = (gpXPrior + oldLL$LL - 0.5 * oldLL$log.det)
  proposed = proposal$changed
  
  # Check if proposal makes sense
  if(all(proposal$params >= 0)) {
  #if(all(proposal$params > 0.0001)){ #&& all(proposal$params < 10)) {
    # Log Likelihood of new parameters 
    newPrior = calculateLogParamPrior(proposal$params[options$inferredParams],auxVars)
    error = F
    if(options$explicit) {
      newLL = calculateLogLikelihoodExplicit(proposal$params, y, timePoints, sigma, auxVars)
      error = newLL$error
    } else {
      newLL = calculateLogLikelihood(proposal$params, gpFit, data, lambda, timePoints, 
                                 auxVars, samplingStrategy, chain, proposal$p)
    }
    
    ratio = (newPrior - oldPrior) + temperature * (newLL$LL - oldLL$LL) + 
            (proposal$oldProb - proposal$newProb)

    if(!error && !is.nan(ratio) && min(ratio, 0) > log(runif(1))) {  
      params = proposal$params
      accept = proposal$changed
      lL = (gpXPrior + newLL$LL - 0.5 * oldLL$log.det)
      x = newLL$x
    }
  } 
   
  return(list(parameters = params, accepted = accept, proposed=proposed, lL=lL, 
              auxVars = auxVars, x=x, oldLL=oldLL))  
}

calculateLogLikelihoodExplicit <- function(params, y, time, sigma, auxVars) {
  solution = solveODE(dim(y)[2], time, auxVars$ode.system, params) 
  
  param.solution = solution$x 

  explicit.ll = 0
 
  for(species in auxVars$observedSpeciesList) {
    explicit.ll = explicit.ll + sum(dnorm(param.solution[,species], y[,species], 
                                      sigma[species], log=T))     
  }  
        
  return(list(LL=explicit.ll, auxVars=auxVars, gpXPrior=0, x=param.solution, log.det=0,
              error=solution$error))
} 
# Calculate Log Likelihood of the model 
calculateLogLikelihood <- function(params, gpFit, X, lambda, timePoints, auxVars, 
                                   samplingStrategy, chain, p = 0, invM = 1, includeDet=F) {
  LL = 0; gpXPrior = 0 
  gradSum = 0; log.det = 0
 
  # Calculate Log Likelihood for each species
  for(species in auxVars$speciesList) {
    LL_temp = calculateLogLikelihoodSpecies(params, gpFit[,species], 
           X, lambda[species], timePoints, auxVars, samplingStrategy, species, chain, p)
    LL = LL + LL_temp$logLikelihood
    
    if(includeDet) {
      log.det = log.det + determinant.matrix(LL_temp$noiseA, logarithm=T)$modulus +
                determinant.matrix(LL_temp$K, logarithm=T)$modulus
    }
    
    gpXPrior = gpXPrior + LL_temp$gpXPrior
    gradSum = gradSum + sum(abs(LL_temp$gradDiff))
  }
  
  p.LL = 0

  # Add likelihood of momentum variable (if applicable)
  if(samplingStrategy == 'HMC') {
    p.LL = (t(p) %*% p) / 2
    LL = LL - p.LL
  } else if(samplingStrategy == 'HMCM') {
    p.LL = (t(p) %*% invM %*% p) / 2
    LL = LL - p.LL
  } 

  return(list(LL=LL, auxVars=auxVars, gpXPrior = gpXPrior, p.LL=p.LL,
              gradSum=gradSum, log.det=log.det))

}

# Calculate Log Likelihood of one species
calculateLogLikelihoodSpecies <- function(params, gpFit, X, lambda, timePoints, auxVars, 
                                   samplingStrategy, species, chain, p = NULL) {
  
  LL = calculateLogLikelihoodMCMC(params, gpFit, X, lambda, timePoints, auxVars, species, chain)  
  
  return(LL)

}

# Sample Observation Noise Variance
sampleNoise <- function(sigma, x, y, options, chain, species, temperatures,auxVars) {
	if(auxVars$sigmaInfer==TRUE){
  sigma.max = sqrt(var(y))

  interval = options$proposalNoiseTuning[chain, species]
  
  if(interval > sigma.max)
    interval = sigma.max - 0.001
  
  sigma.proposal = uniformProposal(sigma, 
                     interval, sigma.max) 
   
  old.prior = dgamma(sigma, shape=0.5, scale=1, log=T)
  new.prior = dgamma(sigma.proposal, shape=0.5, scale=1, log=T)
  
  old.ll = sum(dnorm(y, x, sigma, log=T))
  new.ll = sum(dnorm(y, x, sigma.proposal, log=T))
  
  if(is.nan(old.ll) || is.nan(new.ll))
    browser()

  #ratio = exp(temperatures[chain] * (new.ll - old.ll) + new.prior - old.prior)
  ratio = exp((new.ll - old.ll) + new.prior - old.prior)

  if(min(1,ratio) > runif(1)) {
    sampled.sigma = sigma.proposal
    accept=1
    data.LL = new.ll
  } else {
    sampled.sigma = sigma
    accept = 0
    data.LL = old.ll
  }
  
  return(list(sigma=sampled.sigma, accept=accept, data.LL=data.LL))}

	if(auxVars$sigmaInfer==FALSE){
	sigmaT <- auxVars$sigmaTrue
	accept <- 1
	data.LL <- sum(dnorm(y, x, sigmaT, log=T))
	return(list(sigma=sigmaT, accept=accept, data.LL=data.LL))
	}
} 



# Propose new parameter set
proposeParams <- function(X, timePoints, oldParams, lambda, auxVars, 
                          gpFit, samplingStrategy,
                          options, chain, temperature) {
  if(samplingStrategy == 'MCMC' || samplingStrategy == 'HMCM') {
    res = proposeParamsMCMC(oldParams, options$inferredParams, options$proposalTuning[chain,],
                            auxVars$paramsCovEstimate[[chain]], options$explicit)
  } else if(samplingStrategy == 'HMC') { 
    res = proposeParamsHMC(X, timePoints, oldParams, lambda, auxVars, gpFit, 
                           options$inferredParams, options$proposalTuning[chain,],
                           temperature)

  }

  return(res)
}




