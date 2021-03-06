% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_mcmc.R
\name{doMCMC}
\alias{doMCMC}
\title{Main MCMC function
Runs the MCMC for the specified number of iterations and returns the
sampled parameter values}
\usage{
doMCMC(timePoints, data, auxVars, options)
}
\arguments{
\item{timePoints}{Measured time points for the ODE system.}

\item{data}{Observed data.}

\item{auxVars}{Auxiliary variables.}

\item{options}{Options for MCMC run.}
}
\value{
Function returns a list with elements \code{parameters}, the sampled ODE parameters for the current MCMC iteration, \code{tuning}, the inferred tuning parameters for acceptance of the MCMC moves, \code{paramsRec}, the sampled ODE parameters recorded over all MCMC iterations, \code{lLRec}, the log likelihood recorded over all MCMC iterations, \code{xRec}, the samples from the Gaussian process of the latent variables recorded over all MCMC iterations, \code{gpRec}, the samples of the hyperparameters for the Gaussian process recorded over all MCMC iterations, \code{timePoints}, the time points, \code{noiseRec}, the standard deviation of the observational noise recorded over all MCMC iterations (currently fixed), \code{swappedChains}, the number of times the chains have been swapped, \code{chainNums}, the number of chains, \code{maxIterations}, the total number of MCMC iterations and \code{lLAllChains}, the log likelihood for all chains recorded over all MCMC iterations. If the user specifies \code{temperMismatchParameter=FALSE}, the function additionally returns \code{gradientMismatchParameterRec}, the sampled gradient mismatch parameters recorded over all MCMC iterations.
}
\description{
Main MCMC function
Runs the MCMC for the specified number of iterations and returns the
sampled parameter values
}
