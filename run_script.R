library('devtools')


### Load test dataset
load("LV SD Noise 0.31 Average SNR 10 1.RData")

dataTest <- dataset$data
timeTest <- dataset$time
noiseTest <- dataset$noise

# Simulates installation and loading of package
load_all('deGradInfer')

# Lotka volterra function as input
LV_func = function(t, X, params) {
	dxdt = cbind(
	  X[,1]*(params[1] - params[2]*X[,2]), # Sheep
 	- X[,2]*(params[3] - params[4]*X[,1])  # Wolves
	)
	return(dxdt)
}

#timeTest = seq(0,2,0.1)
#dataTest = ode(c(5,3), timeTest, function(t,y,params) list(LV_func(t,matrix(y,1,length(y)),params)), c(2,1,4,1))
#dataTest = dataTest[,2:3] + rnorm(dim(dataTest)[1]*2,0,0.1)


# AGM
agm(data=dataTest,time=timeTest,ode.system=LV_func,numberOfParameters=4,temperMismatchParameter=TRUE,
   maxIterations=100000,originalSignalOnlyPositive=TRUE, showProgress=TRUE, returnResults=FALSE, noise.sd=0.1,
   defaultPrior="Gamma",defaultTemperingScheme="LB10")

agmResults <- agm(data=dataTest,time=timeTest,ode.system=LV_func,numberOfParameters=4,temperMismatchParameter=TRUE,
   maxIterations=100000,originalSignalOnlyPositive=TRUE, showProgress=TRUE, noise.sd=0.1,
   defaultPrior="Gamma",defaultTemperingScheme="LB10")


# Explicit ODE solution
agm(data=dataTest,time=timeTest,ode.system=LV_func,observedSpecies=1,numberOfParameters=6,temperMismatchParameter=TRUE,
    maxIterations=100000,originalSignalOnlyPositive=TRUE, 
    showProgress=TRUE, returnResults=FALSE, noiseInfer=TRUE,
    defaultPrior="Gamma",defaultTemperingScheme="LB10",explicit=TRUE)


