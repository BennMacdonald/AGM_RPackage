# Initialisation and Utility Functions

# Initialise chains and other things
setupChains <- function(timePoints, data, auxVars, options) {
  chainNum = options$chainNum
  speciesNum = dim(data)[2]
  paramNum = length(options$paramsInit)

	if(auxVars$sigmaInfer==TRUE){
  sigma = matrix(rgamma(speciesNum*chainNum, shape=1, scale=0.1), 
                 chainNum, speciesNum)}

	if(auxVars$sigmaInfer==FALSE){
  sigma = matrix(auxVars$sigmaTrue, 
                 chainNum, speciesNum)}

  if(!('speciesList' %in% names(auxVars)) || 
    is.null(auxVars$speciesList)) {
    auxVars$speciesList = 1:speciesNum
  }
  
  if(!('observedSpeciesList' %in% names(auxVars))) {
    auxVars$observedSpeciesList = auxVars$speciesList
  }

  # Use ML to get GP parameters
  if(options$gpInit) {
    gpFit.temp = fitGP(timePoints, data, auxVars$observedSpeciesList, options, auxVars$covtype)
    if(options$showProgress) print(gpFit.temp) 
  }

  # Setup latent values
  x = list()
  
  auxVars$paramCovEstimate = list()
  paramsTempRec = list()
  proposalGPTuning = list()
  gpFit = list()
  
  for(i in 1:chainNum) {
    x[[i]] = as.matrix(data)
    auxVars$paramsCovEstimate[[i]] = diag(length(options$inferredParams))
    paramsTempRec[[i]] = matrix(0, 0, length(options$inferredParams))
    proposalGPTuning[[i]] = options$proposalGPTuning 
    gpFit[[i]] = matrix(0, 3, speciesNum)
         
    for(species in auxVars$speciesList) {
      if(species %in% auxVars$observedSpeciesList) {
        if(options$gpInit) {
          x[[i]][,species] = gpFit.temp[[species]]$x
          gpFit[[i]][1:length(gpFit.temp[[species]]$params),species] = 
            gpFit.temp[[species]]$params
          sigma[i, species] = gpFit.temp[[species]]$noise
        } else {
          x[[i]][,species] = data[,species]
          gpFit[[i]][,species] = 1
          sigma[i, species] = 0.0001
        }
      } else {
        gpFit[[i]][1:2,species] = 1
        x[[i]][,species] = data[,species]
        sigma[i, species] = 0.0001
      }
    }
  }
 
 
  # Setup precalculated values
  precalc.temp = covUtilPreCalc(timePoints)
  auxVars$diffmatrix = precalc.temp$diffmatrix
  auxVars$sqmatrix = precalc.temp$sqmatrix
  auxVars$tmatrix = precalc.temp$tmatrix
  auxVars$prodmatrix = precalc.temp$prodmatrix

  auxVars$I = diag(dim(data)[1]); 
  
  A.rec <<- list(); invA.rec <<- list();
  K.rec <<- list(); invK.rec <<- list()
  invNoiseA.rec <<- list(); deriv.m.rec <<- list()
  
  # Initialise Parameters
  parameters = t(matrix(options$paramsInit, paramNum, chainNum))
  
  if(!options$keepInit) {
    parameters[,options$inferredParams] = 
      rgamma(chainNum*length(options$inferredParams), 0.5, 1)
    
    if(options$explicit) {
      parameters[,(length(options$paramsInit) - dim(data)[2] + 1):length(options$paramsInit)] = 
        t(matrix(data[1,], speciesNum, chainNum))
    }
  }  

	if (auxVars$Mismatch$Tempering) { # For users who want to specify their own values
		if (!is.null(auxVars$Mismatch$lambdaValues)){
			lambda <- auxVars$Mismatch$lambdaValues
		}
		if (is.null(auxVars$Mismatch$lambdaValues)){ # Using the default LB2 or LB10 values
			lambda <- options$lambda
		}
	}

	if (!auxVars$Mismatch$Tempering){
		lambda = matrix(0.001, chainNum, speciesNum) 
	}


  # Acceptance rate tuning   
  tuning = list()
  tuning$accepted = matrix(0, chainNum, paramNum)
  tuning$acceptTemp = matrix(1, chainNum, paramNum)
  tuning$proposed = matrix(0, chainNum, paramNum)
  tuning$proposeTemp = matrix(1, chainNum, paramNum)
  tuning$acceptedExchange = matrix(0, chainNum, 1)
  tuning$acceptExchangeTemp = matrix(0, chainNum, 1)
  tuning$proposedExchange = matrix(0, chainNum, 1)
  tuning$proposeExchangeTemp = matrix(1, chainNum, 1)
  tuning$acceptXTemp = matrix(0, chainNum, speciesNum)
  tuning$proposeXTemp = matrix(1, chainNum, speciesNum)
  tuning$acceptedX = matrix(0, chainNum, speciesNum)
  tuning$proposedX = matrix(0, chainNum, speciesNum)
  tuning$acceptNoiseTemp = matrix(0, chainNum, speciesNum) 
  tuning$proposeNoiseTemp = matrix(1, chainNum, speciesNum) 
  tuning$acceptedNoise = matrix(0, chainNum, speciesNum) 
  tuning$proposedNoise = matrix(0, chainNum, speciesNum) 
  tuning$acceptGPTemp = list()
  tuning$proposeGPTemp = list()
  tuning$acceptedGP = list()  
  tuning$proposedGP = list()
  
  print('Setup')
  for(chain in 1:chainNum) {
    A.rec[[chain]] <<- matrix(0, dim(data)[1]^2, dim(data)[2])
    invNoiseA.rec[[chain]] <<- matrix(0, dim(data)[1]^2, dim(data)[2])
    K.rec[[chain]] <<- matrix(0, dim(data)[1]^2, dim(data)[2]) 
    invK.rec[[chain]] <<- matrix(0, dim(data)[1]^2, dim(data)[2])
    deriv.m.rec[[chain]] <<- matrix(0, dim(data)[1]^2, dim(data)[2])
    
    for(species in auxVars$speciesList) {
      auxVars$Kchanged = species
      calculateLogLikelihoodMCMC(parameters[chain,], gpFit[[chain]], x[[chain]], lambda[chain,species], timePoints, 
                                             auxVars, species, chain) 
      A.rec[[chain]][,species] <<- A.temp
      invNoiseA.rec[[chain]][,species]  <<- invNoiseA.temp
      K.rec[[chain]][,species]  <<- K.temp
      invK.rec[[chain]][,species]  <<- invK.temp
      deriv.m.rec[[chain]][,species]  <<- deriv.m.temp
    }
    
    auxVars$Kchanged = 0 
    auxVars$lambdaChanged = 0
        
    tuning$acceptGPTemp[[chain]] = matrix(0, 3, speciesNum) 
    tuning$proposeGPTemp[[chain]] = matrix(1, 3, speciesNum)
    tuning$acceptedGP[[chain]] = matrix(0, 3, speciesNum)
    tuning$proposedGP[[chain]] = matrix(1, 3, speciesNum)
  }
  print('End Setup')
  tuning$proposeLambdaTemp = matrix(1, chainNum, speciesNum)
  tuning$acceptLambdaTemp = matrix(0, chainNum, speciesNum)
  tuning$acceptedLambda = matrix(0, chainNum, speciesNum)
  tuning$proposedLambda = matrix(1, chainNum, speciesNum)


  # Set up record of samples for X and GP parameters
  xRec  = list(); gpRec = list()
  
  for(species in auxVars$speciesList) {
    xRec[[species]] =  matrix(0, 0, dim(data)[1])
    gpRec[[species]] =  matrix(0, 0, 3)
  }

  return(list(auxVars=auxVars, tuning=tuning, x=x, gpFit=gpFit, parameters=parameters,
              xRec=xRec, lambda=lambda, sigma=sigma, 
              speciesNum=speciesNum, gpRec=gpRec, paramsTempRec=paramsTempRec,
              proposalGPTuning=proposalGPTuning)) 
  
}

# Reset acceptance and proposal rates for tuning
resetAndUpdateRates <- function(tuning, chainNum, paramNum, speciesNum) {

   tuning$acceptedExchange =  tuning$acceptedExchange + tuning$acceptExchangeTemp
   tuning$acceptExchangeTemp = matrix(0, chainNum, 1)
   tuning$proposedExchange =  tuning$proposedExchange + tuning$proposeExchangeTemp
   tuning$proposeExchangeTemp = matrix(1, chainNum, 1)
      
   tuning$accepted = tuning$accepted + tuning$acceptTemp
   tuning$acceptTemp = matrix(0, chainNum, paramNum)
   tuning$proposed = tuning$proposed + tuning$proposeTemp
   tuning$proposeTemp = matrix(1, chainNum, paramNum)
   tuning$acceptedX = tuning$acceptedX + tuning$acceptXTemp
   tuning$acceptXTemp = matrix(0, chainNum, speciesNum)
   tuning$proposedX = tuning$proposedX + tuning$proposeXTemp
   tuning$proposeXTemp = matrix(1, chainNum, speciesNum)
   tuning$acceptedNoise = tuning$acceptedNoise + tuning$acceptNoiseTemp
   tuning$acceptNoiseTemp = matrix(0, chainNum, speciesNum)
   tuning$proposedNoise = tuning$proposedNoise + tuning$proposeNoiseTemp
   tuning$proposeNoiseTemp = matrix(1, chainNum, speciesNum)

   for(chain in 1:chainNum) {
     tuning$acceptedGP[[chain]] = tuning$acceptedGP[[chain]] + tuning$acceptGPTemp[[chain]]
     tuning$proposedGP[[chain]] = tuning$proposedGP[[chain]] + tuning$proposeGPTemp[[chain]]
     tuning$acceptGPTemp[[chain]] = matrix(0, 3, speciesNum) 
     tuning$proposeGPTemp[[chain]] = matrix(1, 3, speciesNum)
   }

   return(tuning)
}

# Adjust proposals
adjustProposal <- function(tuning, options) {
  acceptParams = tuning$acceptTemp
  proposeParams = tuning$proposeTemp
  acceptX = tuning$acceptXTemp
  proposeX = tuning$proposeXTemp
  acceptNoise = tuning$acceptNoiseTemp
  proposeNoise = tuning$proposeNoiseTemp
  acceptGP = tuning$acceptGPTemp
  proposeGP = tuning$proposeGPTemp

  decrement = 1 - runif(1, 0, 0.5)
  increment = 1 + runif(1, 0, 0.5)
 
  if(options$samplingStrategy == 'HMC' || T) {
    acceptRatio = acceptParams / proposeParams
    acceptXRatio = acceptX / proposeX
    acceptNoiseRatio = acceptNoise / proposeNoise

    tuning  = options$proposalTuning
    tuningX = options$proposalXTuning
    tuningNoise = options$proposalNoiseTuning
    tuningGP = options$proposalGPTuning
    
    for(chain in 1:dim(acceptParams)[1]) {
      
      chainNoiseAccept = acceptNoiseRatio[chain,]
      chainXAccept = acceptXRatio[chain,] 

      for(species in 1:dim(acceptNoise)[2]) {
        
        speciesXAccept = chainXAccept[species]
        if(speciesXAccept < 0.2) {
          tuningX[chain,species] = tuningX[chain,species]*decrement
        } else if(speciesXAccept > 0.3) {
          tuningX[chain,species] = tuningX[chain,species]*increment
        }

        speciesNoiseAccept = chainNoiseAccept[species]
        if(speciesNoiseAccept < 0.2) {
          tuningNoise[chain, species] = tuningNoise[chain, species]*decrement
        } else if(speciesNoiseAccept > 0.3) {
          tuningNoise[chain, species] = tuningNoise[chain, species]*increment
        }

      }
 
      chainAccept = acceptRatio[chain,]

      for(param in 1:length(chainAccept)) { 
        if(chainAccept[param] < 0.2 && tuning[chain, param] > 0.001) { 
          tuning[chain, param] = tuning[chain, param]*decrement
        } else if(chainAccept[param] > 0.3) {
          tuning[chain, param] = tuning[chain, param]*increment 
        }
      }
      

      acceptGPRatio = acceptGP[[chain]]/proposeGP[[chain]]
      acceptGPRatio[proposeGP[[chain]]==1] = 0.25

      for(species in 1:dim(acceptGP[[chain]])[2]) { 
        for(param in 1:dim(acceptGP[[chain]])[1]) { 
          if(acceptGPRatio[param, species] < 0.2) {
            tuningGP[[chain]][param, species] = tuningGP[[chain]][param, species]*decrement
          } else if(acceptGPRatio[param, species] > 0.3) {
            tuningGP[[chain]][param, species] = tuningGP[[chain]][param, species]*increment
          } 
        }
      }
  
    }
    
    
  } else {
    tuning = options$proposalTuning
    tuningX = options$proposalXTuning
    tuningNoise = options$proposalNoiseTuning
  } 
  
  return(list(tuningParams=tuning, tuningX=tuningX, tuningNoise=tuningNoise, tuningGP=tuningGP))
}

# Adjust proposal acceptance rate for exchange moves (by changing temperatures)
adjustExchangeProposal <- function(tuning, temperatures, temp.exponent) {
  accept  = tuning$acceptExchangeTemp
  propose = tuning$proposeExchangeTemp
  
  differences = temperatures[2:length(temperatures)] - 
    temperatures[1:(length(temperatures)-1)]
  
  acceptRatio = accept / propose
  
  for(dist.i in 1:length(differences)) {
    decrement = 1 - runif(1, 0, 0.1)
    increment = 1 + runif(1, 0, 0.1)
     
    dist.temp = differences[dist.i]
    
    if(acceptRatio[dist.i + 1] < 0.45) {
      dist.temp = dist.temp*decrement  
    } else if(acceptRatio[dist.i + 1] > 0.55) {
      dist.temp = dist.temp*increment
    } 
    
    differences[dist.i] = dist.temp
  }
 
  differences = differences / sum(differences)
  temperatures = c(0, cumsum(differences))

  return(temperatures)
}

# Calculate summary statistics for likelihood
likelihoodUtil <- function(params, X, lambda, timePoints, auxVars, gpFit, species, chain) {
  odeNoiseParam = lambda
  error = FALSE
  # Gradient from the ODE system
  f = getODEGradient(X, timePoints, params, auxVars, species)
  f = as.matrix(f) 
  
  # If GP parameters have changed, recalculate K and derivatives, and A
  if(auxVars$Kchanged == species) {
    gpCovs = getGPCovs(gpFit, auxVars)
    invK = NULL
    
    # Try inverse of K
    try(invK <- solve.default(gpCovs$K), silent=T)
    
    if(is.null(invK)) {
      invK = auxVars$I
      error = T
    }
 
    tempA = gpCovs$Kstar %*% invK %*% gpCovs$starK 
    tempA = (tempA + t(tempA)) / 2
  
    A = gpCovs$starKstar - tempA
    deriv.m = gpCovs$Kstar %*% invK 
    K = gpCovs$K
   
    if(!error) {
      A.temp <<- c(A)
      deriv.m.temp <<- c(deriv.m)
      invK.temp <<- c(invK)
      K.temp <<- c(K)
    }
  } else { 
    invK = matrix(invK.rec[[chain]][, species], length(f), length(f))
    deriv.m = matrix(deriv.m.rec[[chain]][, species], length(f), length(f))
    A = matrix(A.rec[[chain]][, species], length(f), length(f))
    K = matrix(K.rec[[chain]][, species], length(f), length(f))
    
  }

  m = deriv.m %*% X[,species,drop=F]

  I = auxVars$I 

  noiseA = A + (odeNoiseParam+1e-3)*I
  return(list(m=m, noiseA=noiseA, gradDiff = f-m, 
              error=error, invK=invK, deriv.m=deriv.m, A=A, K=K))
}

# Precalculate some matrices for GP covariance calculations
covUtilPreCalc <- function(input) {
  inputm = matrix(rep(input, length(input)), length(input), length(input))
  
  diffmatrix = inputm - t(inputm)
  sqmatrix = diffmatrix^2
  prodmatrix = inputm*t(inputm) 

  return (list(diffmatrix=diffmatrix, sqmatrix = sqmatrix, 
               tmatrix=inputm, prodmatrix=prodmatrix))
}

# Calculate covariance matrices for GP with RBF kernel
getGPCovs <- function(gpParams, auxVars) {
   if(auxVars$covtype == 'rbf') {
     res = getGPCovsRBF(gpParams, auxVars$diffmatrix, auxVars$sqmatrix)
   } else if(auxVars$covtype == 'sigmoid') {
     res = getGPCovsSigmoid(gpParams, auxVars$tmatrix, auxVars$prodmatrix)
   } else if(auxVars$covtype == 'sigmoidVar') {
     res = getGPCovsSigmoidVar(gpParams, auxVars$tmatrix, auxVars$prodmatrix)
   }

   return(res)
}

# Calculate covariance matrices for GP with RBF kernel
getGPCovsSigmoidVar <- function(gpParams, tmatrix, prodmatrix) {
  a = gpParams[1]
  b = gpParams[2]
  c = gpParams[3] 

  t1 = tmatrix
  t2 = t(tmatrix)

  num = a + b*prodmatrix
  denom1 = a + b*t1^2 + 1
  denom2 = a + b*t2^2 + 1 
  denom = sqrt(denom1 * denom2)
  Z = num/denom

  K = c*asin(Z)

  # For numerical stability
  K = K + diag(1e-8, dim(tmatrix)[1])
 
  asin.deriv = 1 / sqrt(1 - Z^2)
  Z.deriv1 = b*(t2/denom - Z*t1/denom1)
  Z.deriv2 = b*(t1/denom - Z*t2/denom2)

  Kstar = c * asin.deriv * Z.deriv1 # w.r.t t1
  
  starK = c * asin.deriv * Z.deriv2 # w.r.t. t2
 
  asin.dderiv = Z / (1 - Z^2) 
  
  Z.dderiv = b/denom - b^2*t2^2*denom1/(denom1*denom2*denom) - 
             (b*t1/denom1) * Z.deriv2
 
  starKstar = c * asin.deriv * ( asin.dderiv * Z.deriv1 * Z.deriv2 + Z.dderiv)

  return(list(K=K, starK=starK, Kstar=Kstar, starKstar=starKstar))
}

# Calculate covariance matrices for GP with RBF kernel
getGPCovsSigmoid <- function(gpParams, tmatrix, prodmatrix) {
  a = gpParams[1]
  b = gpParams[2]

  t1 = tmatrix
  t2 = t(tmatrix)

  num = a + b*prodmatrix
  denom1 = a + b*t1^2 + 1
  denom2 = a + b*t2^2 + 1 
  denom = sqrt(denom1 * denom2)
  Z = num/denom

  K = asin(Z)

  # For numerical stability
  K = K + diag(1e-10, dim(tmatrix)[1])
 
  asin.deriv = 1 / sqrt(1 - Z^2)
  Z.deriv1 = b*(t2/denom - Z*t1/denom1)
  Z.deriv2 = b*(t1/denom - Z*t2/denom2)

  Kstar = asin.deriv * Z.deriv1 # w.r.t t1
  
  starK = asin.deriv * Z.deriv2 # w.r.t. t2
 
  asin.dderiv = Z / (1 - Z^2) 
  
  Z.dderiv = b/denom - b^2*t2^2*denom1/(denom1*denom2*denom) - 
             (b*t1/denom1) * Z.deriv2
 
  starKstar = asin.deriv * ( asin.dderiv * Z.deriv1 * Z.deriv2 + Z.dderiv)

  return(list(K=K, starK=starK, Kstar=Kstar, starKstar=starKstar))
}

# Calculate covariance matrices for GP with RBF kernel
getGPCovsRBF <- function(gpParams, diffmatrix, sqmatrix) {
  invcharlength = gpParams[1]
  sigma = gpParams[2]

  K = sigma * exp( - sqmatrix * 0.5 * invcharlength )
 
  # For numerical stability
  K = K + diag(1e-10, dim(diffmatrix)[1])
  
  Kstar = - invcharlength * (diffmatrix) * K
  
  starK = invcharlength * (diffmatrix) * K
  
  starKstar = (invcharlength - invcharlength^2 * sqmatrix) * K
  return(list(K=K, starK=starK, Kstar=Kstar, starKstar=starKstar))
}

