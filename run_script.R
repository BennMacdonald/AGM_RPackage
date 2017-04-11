library('devtools')


### Load test dataset
# Frank: It's good practice to give your R data files the *.RData extension, 
# so I changed the extension

load("LV SD Noise 0.31 Average SNR 10 1.RData")


dataTest <- dataset$data
timeTest <- dataset$time
noiseTest <- dataset$noise

rm(dataset)


# Simulates installation and loading of package
load_all('deGradInfer/')

agm(data=dataTest,time=timeTest,numberOfParameters=4,temperMismatchParameter=TRUE,
    defaultODE="LotkaVolterra",maxIterations=1000,originalSignalOnlyPositive=TRUE,
    defaultPrior="Gamma",defaultTemperingScheme="LB10")


