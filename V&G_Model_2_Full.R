### This script creates the Vyshemirsky And Girolami model 2. It sets up the 
	### odes and begins preparation for the MCMC to be set-up/run. Note:
	### Correction has been made from original publication and this model
	### follows the CORRECTED description


### Load the "deSolve" package

library(deSolve)


### Creating the system of ordinary differential equations for the V&G
	### model 2

### Species: There are 4 species and 3 species from models 1 and 4

### Species in order are:
###  1. S protein
###  2. DS degraded S protein
###  3. R inactive protein
###  4. RS binded protein (R and S proteins)
###  5. Rpp active state protein
###  6. PhA phosphotase
###  7. RppPhA protein


### Parameters: There are 5 kinetic parameters

### Parameters in order are:
### 1. k1
### 2. k2
### 3. k3
### 4. V1
### 5. V2



### Calculating gradients for the VG model. This consists of setting up the
	### system of odes. dydt contains the values of dy by dt
	### (differentiating y with respect to t, where y is the response i.e.
	### the species, and t is time). params is a vector of the parameter
	### values

gradientsVGModel2Full <- function(params,X,timepoints=NULL,auxVars,species
					  )
{

	rat1 <- (params[4]*X[,3]*X[,1])/(params[2]+X[,3])

	rat2 <- (params[5]*X[,5])/(params[3]+X[,5])

	dydt <- switch(species,

		- params[1]*X[,1] 
	,
		params[1]*X[,1]
	,
		- rat1 + rat2
	,
		rep(0,dim(X)[1])
	,
		rat1 - rat2
	,
		rep(0,dim(X)[1])
	,
		rep(0,dim(X)[1])
			  )

	return(dydt)
}


### Creating function odeVG. This will be called in a later stage, during the 
	### population MCMC stage and is similar to the gradient
	### specification above. If using the function "cbind", then
	### each element will be one equation.

odeVGModel2Full <- function(t,X,params)
{

	rat1 <- (params[4]*X[3]*X[1])/(params[2]+X[3])

	rat2 <- (params[5]*X[5])/(params[3]+X[5])


	dydt <- c(

		- params[1]*X[1] 
	,
		params[1]*X[1]
	,
		- rat1 + rat2
	,
		0
	,
		rat1 - rat2
	,
		0
	,
		0

			 )

	grads = list()
	grads[[1]] = dydt

  	return(grads)
}


### Creation of simulation VG model 2. It uses the function "ode" in the cran
	### package "deSolve" to solve the ode V&G model function above. xInit
	### is the variable containing the initial values for the
	### species concentrations, t is the timepoints at which each species is
	### simulated and params is the vector of parameter values

simulateVGModel2Full <- function(xInit,t,params,...)
{
	solution <- ode(xInit,t,odeVGModel2Full,params,...)

	return(solution)
}




