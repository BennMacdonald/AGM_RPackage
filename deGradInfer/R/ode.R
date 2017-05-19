# Utility functions for explicit solving of ODEs

#' Solve ODE system explicitly.
#'
#' @param num.species Number of variables (species) in the system.
#' @param timePoints Time points at which to evaluate the ODE system.
#' @param ode.system Function for calculating the derivatives of the ODE system.
#' @param params Current values for the ODE parameter estimates.
#'
#' @return A list with two elements: \code{x} contains the results of integrating the ODE at the given time points, and \code{error} flags if there has been an error while invoking \code{deSolve}.
#'
solveODE <- function(num.species, timePoints, ode.system, params) {
  initial.conditions = (length(params) - num.species+1):length(params)
  
  solution = ode(params[initial.conditions], timePoints, 
                 function(t, X, params) {
                   # ensure X is a matrix
                   if(length(t)==1) X = matrix(X,1,length(X)) 
                   list(ode.system(t,X,params))
                 }, 
                 params)
  
  # Check for error in solving ODE system
  error = FALSE 
  result = matrix(0, length(timePoints), num.species)
  result.dims = dim(result)
  result.dims[2] = result.dims[2] + 1
  
  if(any(dim(solution) != result.dims)) error = TRUE
  
  result[1:nrow(solution), 2:ncol(solution)-1] = solution[,2:ncol(solution)]
  
  return(list(x=result,error=error))
}

#' Calculate gradients from ODE system
#'
#' @param X Latent values for the species
#' @param timePoints Times at which to calculate the ODE gradients
#' @param params Current parameter estimates
#' @param auxVars Auxiliary variables (including function for ODE gradients)
#' @param species Which species to return (default=all)
#'
#' @return
#' A T by \code{length(species)} matrix with the gradients calculated at each time point for the specified species.
getODEGradient <- function(X, timePoints, params, auxVars, species=1:dim(X)[2]) {

  # Why the offset?
  XODE <- X + auxVars$constant
  
  if (auxVars$originalPositive==TRUE){
    XODE[XODE < 0] <- 0
  }
  
  # Calculate gradients
  # Currently calculates all gradients, then retains only the ones for one species;
  # this will be inefficient for very large ODE systems
  gradient = auxVars$ode.system(timePoints, XODE, params)[,species,drop=FALSE]
  
  return(gradient)
}