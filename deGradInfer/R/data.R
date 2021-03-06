# Documentation for data files

#' Data from a Lotka-Volterra ODE system with 2 species and 4 parameters.
#' Species in order are:
#' 1. Sheep (Prey)
#' 2. Wolves (Predators)
#'
#' @format A list with 5 components
#' \describe{
#'   \item{data}{Observed species with observation noise.}
#'   \item{data.true}{Observed species without observation noise.}
#'   \item{time}{Observed time points.}
#'   \item{params}{Observed parameters.}
#'   \item{noise}{Variance for the Gaussian observation noise.}
#' }
#' @source Generated by function simulateLotkaVolterraModel.
"LV_example_dataset"