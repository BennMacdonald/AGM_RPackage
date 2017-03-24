### This script creates the Lotka-Volterra model. It sets up the odes and
	### begins preparation for the MCMC to be set-up/run


### Load the "deSolve" package

#library(deSolve)


### Creating the system of ordinary differential equations for the
	### Lotka-Volterra model

### Species: There are 2 species

### Species in order are:
### 1. Sheep (Prey)
### 2. Wolves (Predators)


### Parameters: There are 4 parameters
### 1. Alpha
### 2. Beta
### 3. Gamma
### 4. Delta


### Calculating gradients for the Lotka-Volterra model. This consists of 
	### setting up the system of odes. dydt contains the values of dy by dt
	### (differentiating y with respect to t, where y is the response i.e.
	### the species, and t is time). params is a vector of the parameter 
	### values

gradientsLotkaVolterraModel <- function(params,X,timepoints=NULL,
							auxVars,species)
{

	dydt <- switch(species,

		X[,1]*(params[1] - params[2]*X[,2])
	,
		- X[,2]*(params[3] - params[4]*X[,1])

			  		)

	return(dydt)
}


### Creating function odeLotkaVolterra. This will be called in a later stage,
	### during the population MCMC stage and is similar to the gradient
	### specification above. If using the function "c", then each element
	### will be one equation.

odeLotkaVolterraModel <- function(t,X,params)
{
	dydt <- c(

		X[1]*(params[1] - params[2]*X[2])
	,
		- X[2]*(params[3] - params[4]*X[1])

			      )

	grads = list()
	grads[[1]] = dydt

  	return(grads)
}


### Creation of simulation Lotka-Volterra. It uses the function "ode" in the
	### cran package "deSolve" to solve the odeLotkaVolterra model function
	### above. xInit is the variable containing the initial values for the
	### species concentrations, t is the timepoints at which each species is
	### simulated and params is the vector of parameter values

simulateLotkaVolterraModel <- function(xInit,t,params,...)
{
	solution <- ode(xInit,t,odeLotkaVolterraModel,params,...)

	return(solution)
}







