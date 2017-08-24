#' Main function for adaptive gradient matching
#'
#' @description Function agm uses adaptive gradient matching to infer the parameters of a user-defined ODE system from data. For details on AGM,
#' see e.g. Dondelinger et al. (2013), MacDonald (2017).
#'
#'
#' @param data A matrix of observations of the ODE system over time. The number of rows is equal to the number time points and the number of columns is equal to the number of variables in the system.
#' @param time A vector containing the time points at which the observations were made.
#' @param ode.system A function describing the ODE system. See Details for more information.
#' @param numberOfParameters A scalar specifying the number of parameters in the ODE system. If explicitly solving the ODE system, the number of parameters will (usually) be equal to the number of ODE parameters plus the number of initial conditions of the system.
#' @param noise.sd A scalar specifiying the value at which to fix the standard deviation of the observational noise. Default \code{noise.sd=1e-3}.
#' @param observedVariables A vector specifying which variables are observed in the system. Default is \code{observedVariables=1:ncol(data)} (fully observed system).
#' @param temperMismatchParameter Logical: whether tempering of the gradient mismatch parameter be carried out? Default is \code{temperMismatchParameter=FALSE}.
#' @param initialisedParameters A vector containing ODE parameters at which to intialise the MCMC. Can be set as \code{NULL} to initialise with a random draw from the prior distribution. Default is \code{initialisedParameters=NULL}.
#' @param chainNum A scalar specifying the number of parallel temperature chains. Default is \code{chainNum=20}.
#' @param gpCovType A string specifying the choice of kernel for the Gaussian process. Currently, there are two: \code{gpCovType="rbf"} and \code{gpCovType="sigmoidVar"}.
#' @param saveFile A string specifying the path and name of the file containing the result output. Can be set as \code{NULL} to save as "AGM Results.R" to the current working directory. Default is \code{saveFile=NULL}.
#' @param defaultTemperingScheme A string indicating which of the two default gradient mismatch parameter value ladders to use. Choices are \code{defaultTemperingScheme="LB2"} or \code{defaultTemperingScheme="LB10"}. Should only be used when \code{temperMismatchParameter=TRUE}. Default is \code{defaultTemperingScheme=NULL}.
#' @param maxIterations A scalar specifying the number of total MCMC iterations. Default is \code{maxIterations=300000}.
#' @param showPlot Logical: whether plots of the MCMC progress should be displayed. Default is \code{showPlot=FALSE}.
#' @param showProgress Logical: whether \% completion and various parameter values should be printed to the workspace. Default is \code{showProgress=FALSE}.
#' @param mismatchParameterValues A matrix containing user specified values for the gradient mismatch parameter. The number of rows should be equal to \code{chainNum} and the number of columns should be equal to the number of variables in the system. A typical ladder should have the largest value in the first row and the smallest value in the last row. Should only be used when \code{defaultTemperingScheme=NULL}. Default is \code{mismatchParameterValues=NULL}.
#' @param originalSignalOnlyPositive Logical: whether all signals observed should be non-negative. When \code{originalSignalOnlyPositive=TRUE}, any negative values of the sampled interpolant will be set to zero. Default is \code{originalSignalOnlyPositive=FALSE}.
#' @param logPrior A string specifying whether one of the default log priors for the ODE parameters should be used, or a user-specified function. Current choices for the default prior are "Uniform", "Gamma" (shape=4, rate=2) and "Mixed" (3 ODE parameters; N(mean=0, sd=0.4), N(mean=0, sd=0.4). Alternatively the user may specify a function for calculating the prior, see Details below. Default is \code{logPrior='Uniform'}.
#' @param explicit Logical: whether the ODE system should be explicitly solved, rather than doing gradient matching. This means that the Gaussian process model is ignored, and the ODE system is directly fitted to the observed data. Default is \code{explicit=FALSE}.
#' @param explicitNoiseInfer Logical: whether the standard deviation of the observational noise should be inferred when using the method that explicitly solves the ODEs. Only considered when \code{explicit=TRUE}. Default is \code{explicitNoiseInfer=TRUE}.
#'
#' @details
#' The parameters \code{ode.system} should be a function of the form \code{f(t, X, params)} where t is the time point vector for which the derivatives should be calculated, X is a T by p matrix containing the values of the variables in the system at time \code{t}, and params is a vector with the current estimated parameter values. The function should return a matrix with the derivatives of x with respect to time (in the same order as in x). Note that in order to be consistent with the \code{ode} in package \code{deSolve}, we require that the function also works for input at a single time point.
#'
#' For specifying a custom prior on the parameters, the user should write their function to take as input a vector of parameters, and return a vector of log densities for a given parameter set. For example, \code{logPrior = function(params) c(dgamma(params,1,1,log=TRUE)} defines a Gamma parameter prior with shape and scale 1.
#' @return Function returns a list with elements \code{posterior.mean}, the mean of the parameter samples from the posterior (after burning of 1/4 of the number of samples taken), \code{posterior.sd}, the standard deviation, \code{posterior.samples}, the parameter samples, \code{ll} the log likelihood, \code{x.samples}, the samples of the latent variables, \code{gp.samples}, the samples of the GP hyperparameters, \code{noise.samples} the samples of sigma, \code{swappedChains} the number of times the chains have been swapped, \code{ll.all.chain}, the log likelihood for all chains and \code{tuning}, the inferred tuning parameters for acceptance of the MCMC moves.
#' @export
#' @importFrom deSolve ode
#' @import gptk stats graphics
#'
#' @examples
#'
#' dataTest <- LV_example_dataset$data
#' timeTest <- LV_example_dataset$time
#' noiseTest <- LV_example_dataset$noise
#'
#' LV_func = function(t, X, params) {
#' 	dxdt = cbind(
#' 	  X[,1]*(params[1] - params[2]*X[,2]),
#'  	- X[,2]*(params[3] - params[4]*X[,1])
#' 	)
#' 	return(dxdt)
#' }
#'
#' param.result = agm(data=dataTest,time=timeTest,noise.sd=0.31,ode.system=LV_func,
#'     numberOfParameters=4,temperMismatchParameter=TRUE,
#'     chainNum=5, maxIterations=200,originalSignalOnlyPositive=TRUE,
#'     logPrior="Gamma",defaultTemperingScheme="LB10")
#'
#' print(param.result$posterior.mean)
#'
agm <- function(data,time,ode.system,numberOfParameters,noise.sd=1e-3, observedVariables=1:ncol(data),
                temperMismatchParameter=FALSE,
                initialisedParameters=NULL,
                chainNum=20,gpCovType="rbf",saveFile=NULL,
                defaultTemperingScheme=NULL,maxIterations=300000,showPlot=FALSE,
                showProgress=FALSE,mismatchParameterValues=NULL,
                originalSignalOnlyPositive=FALSE,logPrior='Uniform',
		            explicit=FALSE, explicitNoiseInfer=TRUE)
{ # Start function agm

  ### For the time being, users will not be able to use the option to infer the
  ### standard deviation/variance of the observational noise. This intention is to
  ### avoid issues of flattening (interpolant). Since this is not an issue for the
  ### explicit solution, we want users to be able to have the option to infer the noise
  ### or hold it fixed when using the explicit solution (new arguement explicitNoiseInfer
  ### added)

  noiseInfer <- FALSE

  if(explicit==TRUE && explicitNoiseInfer==TRUE)
  {
    noiseInfer <- TRUE
  }

  ### Make sure that lambda values will only come from one place (either the
  ### default LB2 and LB10 schemes or through the user's specifications

  if (!is.null(mismatchParameterValues)){
    defaultTemperingScheme <- NULL
  }

  temperatureExponentChains <- 5 # Friel and Pettit

  time <- as.matrix(time)



  dataConstants <- matrix(NA,nrow=nrow(data),ncol=ncol(data))
  for (i in 1:ncol(dataConstants))
  {
    dataConstants[,i] <- mean(data[,i])
  }

  dataset <- data - dataConstants


  covtypes <- c("sigmoidVar","rbf","rbf","sigmoidVar","sigmoidVar","sigmoidVar",
                "sigmoidVar","sigmoidVar","rbf","rbf","rbf","sigmoidVar",
                "rbf")



  ### Specify the parameters to be inferred.

  inferredParams <- 1:numberOfParameters

  ### Specifiy how the initial parameters are chosen. Sometimes it is desired to
  ### set values of the initial parameters, so that the MCMC algorithm has
  ### a good place to start and sometimes it is desired to randomly draw
  ### the starting values.

  if (is.null(initialisedParameters)){
    paramsInit <- c()

    for (i in 1:numberOfParameters)
    {
      paramsInit[i] <- rgamma(1,4,2)
    }
    # The above is for the initial paramaters to be randomly initialised,
    # drawn from a gamma; shape = 4, rate = 0.5
    # Note that users can use the initialisedParameters arguement to set
    # whatever values they want. Hence, if they wish to draw from a
    # different distriubtion, they can do so before using the package
    # and then pass those values to the arguement "initialisedParameters"

  } else{
    paramsInit <- initialisedParameters
  }


  ### Specify the tuning that will take place in the MCMC algorithm

  proposalTuning <- matrix(0.01,nrow=chainNum,ncol=numberOfParameters)
  # A matrix of tuning values that have worked well before, with the
  # matrix having the same number of columns as the number of
  # parameters

  proposalXTuning <- matrix(0.01,nrow=chainNum,ncol=dim(dataset)[2])
  # A matrix of X tuning values that have worked well before, with the
  # matrix having the same number of columns as the number of
  # species

  proposalNoiseTuning <- matrix(0.1,nrow=chainNum,ncol=dim(dataset)[2])
  # A matrix of noise tuning values that have worked well before, with the
  # matrix having the same number of columns as the number of
  # species

  proposalGPTuning <- matrix(0.05,nrow=3,ncol=dim(dataset)[2])
  # A matrix of GP tuning values that have worked well before, with the
  # matrix having the same number of columns as the number of
  # species. Here nrow = 3 since currently the only GP
  # covariance kernels that are supported have < 4 hyperparameters

  ### Set-up lambda if tempering using the LB2 or LB10 schemes

  if (temperMismatchParameter==TRUE && !is.null(defaultTemperingScheme)){
    if (defaultTemperingScheme=="LB10"){
      Base <- 10
      # The base that will be used to create the increments for the log scale
    }
    if (defaultTemperingScheme=="LB2"){
      Base <- 2
      # The base that will be used to create the increments for the log scale
    }

    startLambda <- 1
    # The starting value for lambda (gamma), corresponding to the lowest
    # chain/temperature

    LambdaTemperatures <- c()
    LambdaTemperatures[1] <- startLambda

    for (i in 1:(chainNum-1))
    {
      LambdaTemperatures[i+1] <- LambdaTemperatures[1]*(Base^i)
    }

    Lambda <- matrix(LambdaTemperatures,nrow=chainNum,ncol=dim(dataset)[2])
    # Lambda is now a matrix of the different Lambda temperatures, with
    # number of rows corresponding to number of chains/temperatures
    # and number of columns corresponding to the number of species

    Lambda <- 1/Lambda
  }


  ### Set up auxiliary variables

  auxVars <- list(speciesList=1:ncol(dataset),originalPositive=originalSignalOnlyPositive,
                  Mismatch=list(Tempering=temperMismatchParameter,
                                lambdaValues=mismatchParameterValues),
                  covtype=gpCovType,observedSpeciesList=observedVariables,
                  constant=dataConstants,sigmaInfer=noiseInfer,
			           logPrior=logPrior)

  auxVars$ode.system = ode.system

  if (!is.null(noise.sd)){
    auxVars$sigmaTrue <- noise.sd
  }


  ### Specify the options for the population MCMC algorithm to operate. The
  ### below are the defaults

  options <- list(inferredParams=inferredParams,samplingStrategy="MCMC",
                  samplingStrategyX="mixed",iterations=maxIterations,
                  proposalTuning=proposalTuning,inferX=TRUE,
                  proposalXTuning=proposalXTuning,proposalNoiseTuning=proposalNoiseTuning,
                  proposalGPTuning=proposalGPTuning,burnin=1,
                  chainNum=chainNum,paramsInit=paramsInit,
                  noiseBurnin=1,showPlot=showPlot,
                  saveFile=saveFile,lambda.max=0.1,
                  showProgress=showProgress,keepInit=TRUE,
                  explicit=explicit,allowNeg=T,
                  gpInit=T,temps=temperatureExponentChains
  )

  if (temperMismatchParameter==TRUE && !is.null(defaultTemperingScheme)){
    options$lambda <- Lambda
  }

  ### Run population MCMC algorithm, Gaussian Process approach with Gradient
  ### Matching. At present, only the top chain is saved (the posterior) and
  ### the algorithm itself runs iteratively. Every 10,000
  ### iterations, the function saves the results of the MCMC to the
  ### location specified in "saveFile"

  paramsMCMC <- doMCMC(time,dataset,auxVars,options)

  # Clean up output
  after.burnin = paramsMCMC$paramsRec[
    round(dim(paramsMCMC$paramsRec)[1]/4):dim(paramsMCMC$paramsRec)[1],]

  result = list(posterior.mean=colMeans(after.burnin),
                posterior.sd=apply(after.burnin,2,sd),
                posterior.samples=paramsMCMC$paramsRec,
                ll=paramsMCMC$lLRec,
                x.samples=paramsMCMC$xRec,
                gp.samples=paramsMCMC$gpRec,
                noise.samples=paramsMCMC$noiseRec,
                swappedChains=paramsMCMC$swappedChains,
                ll.all.chains=paramsMCMC$lLAllChains,
                tuning=paramsMCMC$tuning)

  ### Save the final results of the MCMC to the location specified in "saveFile"
  if (!is.null(saveFile)) save(result,file=saveFile)

  return(result)

} # End function agm

