# Tests that for a reasonable number of iterations, the log likelihood improves
library(deGradInfer)

set.seed(10)

context("Likelihood Improvement, 500 Iterations")

# Run deGradInfer on test dataset

dataTest <- LV_example_dataset$data
timeTest <- LV_example_dataset$time
noiseTest <- LV_example_dataset$noise

LV_func = function(t, X, params) {
  dxdt = cbind(
    X[,1]*(params[1] - params[2]*X[,2]),
    - X[,2]*(params[3] - params[4]*X[,1])
  )
  return(dxdt)
}

result.file = paste("",getwd(),"/tests_temp.Rdata",sep="")

agm(data=dataTest,time=timeTest,noise.sd=0.31,ode.system=LV_func,
     numberOfParameters=4,temperMismatchParameter=TRUE, saveFile=result.file,
     showPlot=FALSE,
     chainNum=5, maxIterations=500,originalSignalOnlyPositive=TRUE,
     defaultPrior="Gamma",defaultTemperingScheme="LB10")

load(result.file)
latest = length(paramsMCMC$lLRec)


test_that("Likelihood is improved", {
  expect_gt(paramsMCMC$lLRec[latest], paramsMCMC$lLRec[2])
})
