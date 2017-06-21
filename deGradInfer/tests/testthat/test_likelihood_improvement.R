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

agm(data=dataTest,time=timeTest,noiseFixed=0.31,ode.system=LV_func,
     numberOfParameters=4,temperMismatchParameter=TRUE, saveFile=result.file,
     showPlot=FALSE,
     chainNum=5, maxIterations=500,originalSignalOnlyPositive=TRUE,
     defaultPrior="Gamma",defaultTemperingScheme="LB10")

load(result.file)
latest = length(paramsMCMC$lLRec)


test_that("Likelihood is improved", {
  expect_gt(paramsMCMC$lLRec[latest], paramsMCMC$lLRec[2])
})

# Temporary test to check everything remains the same
test_that("Something changed", {
  expect_lt(abs(paramsMCMC$lLRec[latest] - -6.381104), 1e-5)
  expect_lt(abs(paramsMCMC$parameters[5,1]-2.290987), 1e-5)
  expect_lt(abs(paramsMCMC$parameters[5,2]-1.155992), 1e-5)
  expect_lt(abs(paramsMCMC$parameters[5,3]-1.722465), 1e-5)
  expect_lt(abs(paramsMCMC$parameters[5,4]-0.6981029), 1e-5)
})
