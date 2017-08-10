# Tests that the code does/does not produce an error for different arguments
# to the agm function
library(deGradInfer)

set.seed(10)

context("Arguments Tests, 50 Iterations")

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


# Running without specifying default priors, should work.
test_that("Default prior works",{
expect_silent(agm(data=dataTest,time=timeTest,noise.sd=0.31,ode.system=LV_func,
     numberOfParameters=4, saveFile=result.file,maxIterations=50))
})

# Running without specifying ODE system, should not work.
test_that("ODE system missing throws error",{
  expect_error(agm(data=dataTest,time=timeTest,noise.sd=0.31,
                   numberOfParameters=4, saveFile=result.file,
                   maxIterations=50))
})

