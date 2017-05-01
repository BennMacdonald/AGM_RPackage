library('devtools')


### Load test dataset
load("LV SD Noise 0.31 Average SNR 10 1.RData")

dataTest <- dataset$data
timeTest <- dataset$time
noiseTest <- dataset$noise

# Simulates installation and loading of package
load_all('deGradInfer/')

# Lotka volterra function as input
LV_func = function(t, X, params) {
	dxdt = cbind(
	  X[,1]*(params[1] - params[2]*X[,2]), # Sheep
 	- X[,2]*(params[3] - params[4]*X[,1])  # Wolves
	)
	return(dxdt)
}


# AGM
agm(data=dataTest,time=timeTest,ode.system=LV_func,numberOfParameters=4,temperMismatchParameter=TRUE,
   maxIterations=1000,originalSignalOnlyPositive=TRUE,
   defaultPrior="Gamma",defaultTemperingScheme="LB10")

# Explicit ODE solution
agm(data=dataTest,time=timeTest,ode.system=LV_func,numberOfParameters=6,temperMismatchParameter=TRUE,
    maxIterations=100000,originalSignalOnlyPositive=TRUE, 
    showProgress=TRUE, noiseInfer=TRUE,
    defaultPrior="Gamma",defaultTemperingScheme="LB10",explicit=TRUE)


