### This script creates the Fitz-Hugh Nagumo model from the Liang & Wu paper,
	### 2008, "Parameter estimation for differential equation models, using a
	### a framework of measurement error in regression models". It sets up the 
	### odes and begins preparation for the MCMC to be set-up/run


### Load the "deSolve" package

library(deSolve)


### Creating the system of ordinary differential equations for the
	### Fitz-Hugh Nagumo model

### Species: There are 2 species

### Species in order are:
### 1. V (Voltage)
### 2. R (Recovery Variable)


### Parameters: There are 3 parameters
### 1. Alpha
### 2. Beta
### 3. Gamma


### Calculating gradients for the Fitz-Hugh Nagumo model. This consists of 
	### setting up the system of odes. dydt contains the values of dy by dt
	### (differentiating y with respect to t, where y is the response i.e.
	### the species, and t is time). params is a vector of the parameter 
	### values

gradientsFitzHughNagumoModel <- function(params,X,timepoints=NULL,
							auxVars,species)
{

	dydt <- switch(species,

		params[3]*(X[,1] + X[,2] - ((X[,1]^3)/3))
	,
		- (X[,1] - params[1] + params[2]*X[,2])/params[3]

			  		)

	return(dydt)
}


### Creating function odeFitzHughNagumo. This will be called in a later stage,
	### during the population MCMC stage and is similar to the gradient
	### specification above. If using the function "c", then each element
	### will be one equation.

odeFitzHughNagumoModel <- function(t,X,params)
{
	dydt <- c(

		params[3]*(X[1] + X[2] - ((X[1]^3)/3))
	,
		- (X[1] - params[1] + params[2]*X[2])/params[3]

			      )

	grads = list()
	grads[[1]] = dydt

  	return(grads)
}


### Creation of simulation Fitz-Hugh Nagumo. It uses the function "ode" in the
	### cran package "deSolve" to solve the odeFitzHughNagumo model function
	### above. xInit is the variable containing the initial values for the
	### species concentrations, t is the timepoints at which each species is
	### simulated and params is the vector of parameter values

simulateFitzHughNagumoModel <- function(xInit,t,params,...)
{
	solution <- ode(xInit,t,odeFitzHughNagumoModel,params,...)

	return(solution)
}










