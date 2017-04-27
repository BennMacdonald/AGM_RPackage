# Main function for adaptive gradient matching

#' Title
#'
#' @param data 
#' @param time 
#' @param ode.system A function describing the ODE system. See Details for more information.
#' @param numberOfParameters 
#' @param temperMismatchParameter 
#' @param initialisedParameters 
#' @param noiseInfer 
#' @param noiseFixed 
#' @param chainNum 
#' @param gpCovType 
#' @param saveFile 
#' @param defaultTemperingScheme 
#' @param maxIterations 
#' @param showPlot 
#' @param showProgress 
#' @param mismatchParameterValues 
#' @param originalSignalOnlyPositive 
#' @param defaultPrior 
#' @param explicit Explicitly solve the ODE system, rather than doing gradient matching. This means that the Gaussian process model is ignored, and the ODE system is directly fitted to the observed data. Default \code{explicit=FALSE}.
#'
#' @details 
#' The parameters \code{ode.system} should be a function of the form \code{f(t, X, params)} where t is the time points for which the derivatives should be calculated, X is T by p matrix containing the values of the variables in the system at time \code{t}, and params is a vector with the current estimated parameter values. The function should return a matrix with the derivatives of x with respect to time (in the same order as in x). Note that in order to be consistent with the \code{ode} in package \code{deSolve}, we require that the function also works for input at a single time point. 
#' @return
#' @export
#'
#' @examples
#' 
#' load("LV SD Noise 0.31 Average SNR 10 1.RData")
#' dataTest <- dataset$data
#' timeTest <- dataset$time
#' noiseTest <- dataset$noise
#' 
#' LV_func = function(t, X, params) {
#' 	dxdt = cbind(
#' 	  X[,1]*(params[1] - params[2]*X[,2]),
#'  	- X[,2]*(params[3] - params[4]*X[,1])
#' 	)
#' 	return(dxdt)
#' }
#' 
#' agm(data=dataTest,time=timeTest,ode.system=LV_func,numberOfParameters=4,
#'     temperMismatchParameter=TRUE,
#'     maxIterations=1000,originalSignalOnlyPositive=TRUE,
#'     defaultPrior="Gamma",defaultTemperingScheme="LB10")
#' 
agm <- function(data,time,ode.system,numberOfParameters,temperMismatchParameter=FALSE,
                initialisedParameters=NULL,noiseInfer=TRUE,
                noiseFixed=NULL,chainNum=20,gpCovType="rbf",saveFile=NULL,
                defaultTemperingScheme=NULL,maxIterations=300000,showPlot=TRUE,
                showProgress=FALSE,mismatchParameterValues=NULL,
                originalSignalOnlyPositive=FALSE,defaultPrior=NULL, explicit=FALSE)
{ # Start function agm
  
  ### Make sure that lambda values will only come from one place (either the
  ### default LB2 and LB10 schemes or through the user's specifications
  
  if (!is.null(mismatchParameterValues)){
    defualtTemperingScheme <- NULL
  }
  
  temperatureExponentChains <- 5 # Friel and Pettit
  
  time <- as.matrix(time)
  
  
  
  dataConstants <- matrix(,nrow=nrow(data),ncol=ncol(data))
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
  
  ### Specify the save file location if none given, for the results of the
  ### population MCMC. Defaults to working directory
  
  if (is.null(saveFile)){
    saveFile <- paste("",getwd(),"/AGM Results.R",sep="")
  }
  
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
                  covtype=gpCovType,observedSpeciesList=1:ncol(dataset),
                  constant=dataConstants,sigmaInfer=noiseInfer)
  
  auxVars$ode.system = ode.system
  
  if (!is.null(noiseFixed)){
    auxVars$sigmaTrue <- noiseFixed
  }
  
  if(!is.null(defaultPrior)){
    auxVars$defaultLogParamPrior <- defaultPrior
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
                  showProgress=showProgress,keepInit=F,
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
  
  ### Save the final results of the MCMC to the location specified in "saveFile"
  
  save(paramsMCMC,file=saveFile)
  
  
} # End function agm

