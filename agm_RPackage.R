agm <- function(data,time,numberOfParameters,temperMismatchParameter=FALSE,
			initialisedParameters=NULL,defaultODE=NULL,noiseInfer=TRUE,
			noiseFixed=NULL,chainNum=20,gpCovType="rbf",saveFile=NULL,
			defaultTemperingScheme=NULL,maxIterations=300000,showPlot=TRUE,
			showProgress=FALSE,mismatchParameterValues=NULL,
			originalSignalOnlyPositive=FALSE,defaultPrior=NULL)
{ # Start function agm




### Make sure ALL packages are installed and loaded
source("required_packages.R")
required_packages(c("gdata","gptk","deSolve"))


source("population_mcmc_RPackage.R")

if (!temperMismatchParameter){
	source("samplingMismatchParameter_RPackage.R")
} 


if (!is.null(defaultPrior)) {
	source("log_prior_ODE_params_RPackage.R")
}

if (!is.null(defaultODE)){
	source("gradients_ODEs_RPackage.R")
}


### Make sure that lambda values will only come from one place (either the
	### default LB2 and LB10 schemes or through the user's specifications

if (!is.null(mismatchParameterValues)){
	defualtTemperingScheme <- NULL
}


if(defaultODE=="LotkaVolterra"){
	source("Lotka_Volterra_Model.R")
	modelName <- "LotkaVolterraModel"
}
if(defaultODE=="FitzHughNagumo"){
	source("Fitz_Hugh_Nagumo_Model.R")
	modelName <- "FitzHughNagumoModel"

}
if(defaultODE=="VG1Full"){
	source("V&G_Model_1_Full.R")
	modelName <- "VGModel1Full"
}
if(defaultODE=="VG2Full"){
	source("V&G_Model_2_Full.R")
	modelName <- "VGModel2Full"
}
if(defaultODE=="VG3Full"){
	source("V&G_Model_3_Full.R")
	modelName <- "VGModel3Full"
}
if(defaultODE=="VG4Full"){
	source("V&G_Model_4_Full.R")
	modelName <- "VGModel4Full"
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


if (!is.null(defaultODE)){
	auxVars$modelName <- modelName
}

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
			explicit=F,allowNeg=T,
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
#browser()
paramsMCMC <- doMCMC(time,dataset,auxVars,options)

### Save the final results of the MCMC to the location specified in "saveFile"

save(paramsMCMC,file=saveFile)


} # End function agm

