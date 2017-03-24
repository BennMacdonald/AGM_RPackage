#library('gptk')

# Fit a GP to the one-dimensional data, using the gptk package
# time and y need to be of dimensionality Nx1 where N is the number of datapoints
#
# returns parameters variance and inverse width for rbfkernel = var*exp(-(x1-x2)^2*0.5*invwidth)
# as well as the standard deviation of the associated white noise parameter
fitGP <- function(time, y, speciesList, optionsMCMC, covtype) {
  options = gpOptions()
  if(covtype == 'mixture') {
    options$kern$comp = list("sigmoidVar", "white", "rbf")
  } else {
    options$kern$comp = list(covtype, "white")
  }
  #options$kern = 'sigmoid'
  options$optimiser = "SCG"
  gpFit = list()
  
  if(optionsMCMC$showPlot && length(speciesList) <=6) 
    par(mfrow=c(ceiling(length(speciesList)/2),2))
  else if(optionsMCMC$showPlot)
    par(mfrow=c(1,2))
 
  for(species in speciesList) { 
    # Initialise
    model = gpCreate(1, 1, time, y[,species,drop=F], options)
    
    # Find best parameters with ML
    model = fitModel(model, display=0, covtype)

    gpFit[[species]] = list()
    gpFit[[species]]$params = 
        kernExtractParam(model$kern$comp[[1]], untransformed.values=T)
    
    # DEBUG
    timeC = seq(time[1], time[length(time)], 0.1)
    time.temp = matrix(timeC, length(timeC), 1)
    res = gpPosteriorMeanVar(model, time.temp, varsigma.return = TRUE)

    if(optionsMCMC$showPlot) gpPlot(model, time.temp, res$mu, res$varsigma)

    res = gpPosteriorMeanVar(model, time, varsigma.return = TRUE)
    
    gpFit[[species]]$x = res$mu
    gpFit[[species]]$noise = sqrt(model$kern$comp[[2]]$variance)
    
    Sys.sleep(0.005)
  }

  return(gpFit)
}

fitModel <- function(model, display, covtype) {
    starts = 3
    ll.rec = matrix(0, starts, 1)
    models = list()

    for(i in 1:starts) {
      model$kern$comp[[2]]$variance = runif(1, 0, 0.001)

      if(covtype == 'rbf') {
        model$kern$comp[[1]]$inverseWidth = runif(1, 0, 0.1) 
        model$kern$comp[[1]]$variance = runif(1, 0, 1) 
      } else if(covtype == 'sigmoid') {
        model$kern$comp[[1]]$a = runif(1, 0, 1)
        model$kern$comp[[1]]$b = runif(1, 0, 1)
      } else if(covtype == 'sigmoidVar') {
        model$kern$comp[[1]]$a = runif(1, 0, 10)
        model$kern$comp[[1]]$b = runif(1, 0, 10)
        model$kern$comp[[1]]$var = runif(1, 0, 2*var(model$y))
      } else if(covtype == 'mixture') {
        #model$kern$comp[[3]]$inverseWidth = runif(1, 0, 0.01) 
        #model$kern$comp[[3]]$variance = runif(1, 0, 0.01) 
        model$kern$comp[[1]]$a = runif(1, 0, 1)
        model$kern$comp[[1]]$b = runif(1, 0, 1)
        model$kern$comp[[1]]$var = runif(1, 0, 0.01)
        model$kern$comp[[2]]$a = runif(1, 0, 1)
        model$kern$comp[[2]]$b = runif(1, 0, 1)
        model$kern$comp[[2]]$var = runif(1, 0, 0.01)

      }
 
      
      models[[i]] = gpOptimise(model, display=display)
      ll.rec[i] = gpLogLikelihood(models[[i]]) 
    }
    print(ll.rec)
    best = which.max(ll.rec)
    
    return(models[[best]]) 
}

# Fit a GP to the one-dimensional data, using the gptk package
# time and y need to be of dimensionality Nx1 where N is the number of datapoints
#
# returns parameters variance and inverse width for rbfkernel = var*exp(-(x1-x2)^2*0.5*invwidth)
fitGP_experimental <- function(time, y, save.plots=F, filename='', 
                               species.names=NULL, covtype='sigmoidVar', plot.title='') {
  options = gpOptions()
  options$kern$comp = list(covtype, "white")
  #options$kern = 'sigmoid'
  options$optimiser = "SCG"
  gpFit = list()

  if(is.null(species.names)) {
    species.names = paste('Species', 1:dim(y)[2])
  } 

  for(species in 1:dim(y)[2]) { 
    model = gpCreate(1, 1, time, y[,species,drop=F], options)

    model = fitModel(model, display=0, covtype)
    
    # DEBUG
    timeC = seq(time[1], time[length(time)], 0.1)
    time.temp = matrix(timeC, length(timeC), 1)
    res = gpPosteriorMeanVar(model, time.temp, varsigma.return = TRUE)
    
    if(save.plots) {
      postscript(file=
        paste(filename, '_species_', species, '.eps', sep=''), pointsize=25)
    }

    gpPlot(model, time.temp, res$mu, res$varsigma, ylab=species.names[species],
           xlab='Time')
    title(plot.title)
    
    
    if(save.plots) {
      dev.off()
    } else {
      browser()
    }
  }
  
  return(gpFit)
}

sigmoidKernParamInit <- function(kern) {
  kern$a <- 0.1 
  kern$b <- 0.1 
  kern$nParams <- 2
  kern$paramNames <- c("a", "b")
  kern$isStationary <- F
  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && 
    kern$options$isNormalised) 
    kern$isNormalised <- TRUE
  else kern$isNormalised <- FALSE
  
  kern$transforms <- list(list(index = c(1, 2), type = "positive"))
  return(kern)
}

sigmoidKernExtractParam <- function (kern, only.values = TRUE, 
                                      untransformed.values = TRUE) {
  params <- c(kern$a, kern$b)
  
  if (!only.values) 
    names(params) <- c("a", "b")
  return(params)
}

sigmoidKernExpandParam <- function (kern, params)  {
  if (is.list(params)) 
    params <- params$values
  kern$a <- params[1]
  kern$b <- params[2]
  return(kern)
}

sigmoidKernCompute <- function (kern, x, x2 = NULL) {

  if (nargs() < 3) {
    n2 <- .prod2(x, x) 
    n11 = matrix(x^2, length(x), length(x))
    n22 = t(n11)
  }
  else {
    n2 <- .prod2(x, x2)
    n11 = matrix(x^2, length(x), length(x2))
    n22 = t(matrix(x2^2, length(x2), length(x)))
  }
  
  k <- asin((kern$a + kern$b*n2)/sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))) 

  return(k)
}

sigmoidKernGradient <- function (kern, x, x2, covGrad) 
{
  
  n11 = matrix(x^2, length(x), length(x))

  if (nargs() == 3) {
    n2 <- .prod2(x, x) 
    n22 = t(n11)
    covGrad <- x2
  }
  else if (nargs() == 4) {
    n2 <- .prod2(x, x2) 
    n22 = t(matrix(x2^2, length(x2), length(x2)))
  }
  
  num = kern$a + kern$b*n2
  denom = sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))
  ratio = num/denom
  
  #ratio[ratio==1] = 0.99999999
  
  deriv.asin = 1 / sqrt(1 - (ratio)^2)
  deriv.denom = denom^3
  
  g <- array()
  
  g[1] <- sum(covGrad * deriv.asin * 
    (1/denom - 0.5*num*(2*kern$a + 2 + kern$b*n11 + kern$b*n22)/deriv.denom))
  
  g[2] <- sum(covGrad * deriv.asin * 
    (n2/denom - 0.5*num*((kern$a + kern$b*n11 + 1)*n22 + 
                         (kern$a + kern$b*n22 + 1)*n11)/deriv.denom))
 
  if (any(is.nan(g))) {
    warning("g is NaN.")
    browser()
  }
  
  return(g)
}

sigmoidKernDiagCompute <- function (kern, x) 
{
  x2 = x^2
  k <- matrix(asin((kern$a + kern$b*x2)/(kern$a + 1 + kern$b*x2)), 
              dim(as.array(x))[1], 1)

  return(k)
}

sigmoidVarKernParamInit <- function(kern) {
  kern$a <- 1
  kern$b <- 1  
  kern$var <- 1
  kern$nParams <- 3
  kern$paramNames <- c("a", "b", "var")
  kern$isStationary <- F
  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && 
    kern$options$isNormalised) 
    kern$isNormalised <- TRUE
  else kern$isNormalised <- FALSE
  
  kern$transforms <- list(list(index = c(1, 2, 3), type = "positive"))
  return(kern)
}

sigmoidVarKernExtractParam <- function (kern, only.values = TRUE, 
                                      untransformed.values = TRUE) {
  params <- c(kern$a, kern$b, kern$var)
  
  if (!only.values) 
    names(params) <- c("a", "b", "var")
  return(params)
}

sigmoidVarKernExpandParam <- function (kern, params)  {
  if (is.list(params)) 
    params <- params$values
  kern$a <- params[1]
  kern$b <- params[2]
  kern$var <- params[3]
  return(kern)
}

sigmoidVarKernCompute <- function (kern, x, x2 = NULL) {

  if (nargs() < 3) {
    n2 <- .prod2(x, x) 
    n11 = matrix(x^2, length(x), length(x))
    n22 = t(n11)
  }
  else {
    n2 <- .prod2(x, x2)
    n11 = matrix(x^2, length(x), length(x2))
    n22 = t(matrix(x2^2, length(x2), length(x)))
  }
  
  k <- (2/pi) * kern$var * asin((kern$a + kern$b*n2)/sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))) 

  return(k)
}

sigmoidVarKernGradient <- function (kern, x, x2, covGrad) 
{
  
  n11 = matrix(x^2, length(x), length(x))

  if (nargs() == 3) {
    n2 <- .prod2(x, x) 
    n22 = t(n11)
    covGrad <- x2
    k = sigmoidVarKernCompute(kern, x)
  }
  else if (nargs() == 4) {
    n2 <- .prod2(x, x2) 
    n22 = t(matrix(x2^2, length(x2), length(x2)))
    k = sigmoidVarKernCompute(kern, x, x2)
  }
  
  num = kern$a + kern$b*n2
  denom = sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))
  ratio = num/denom
  
  deriv.asin = 1 / sqrt(1 - (ratio)^2)
  deriv.denom = denom^3
  
  g <- array()
  
  g[1] <- kern$var * sum(covGrad * deriv.asin * 
    (1/denom - 0.5*num*(2*kern$a + 2 + kern$b*n11 + kern$b*n22)/deriv.denom))
  
  g[2] <- kern$var * sum(covGrad * deriv.asin * 
    (n2/denom - 0.5*num*((kern$a + kern$b*n11 + 1)*n22 + 
                         (kern$a + kern$b*n22 + 1)*n11)/deriv.denom))
 
  g[3] <- sum(covGrad * k) / kern$var
  
  if (any(is.nan(g))) {
    warning("g is NaN.")
    browser()
  }
  
  return(g)
}

sigmoidVarKernDiagCompute <- function (kern, x) 
{
  x2 = x^2
  k <- matrix((2/pi) * kern$var * asin((kern$a + kern$b*x2)/(kern$a + 1 + kern$b*x2)), 
              dim(as.array(x))[1], 1)

  return(k)
}

sigmoidWhiteKernParamInit <- function(kern) {
  kern$a <- 1
  kern$b <- 0.1
  kern$c <- 0.01 
  kern$nParams <- 3
  kern$paramNames <- c("a", "b", "c")
  kern$isStationary <- F
  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && 
    kern$options$isNormalised) 
    kern$isNormalised <- TRUE
  else kern$isNormalised <- FALSE
  
  kern$transforms <- list(list(index = c(1, 2), type = "positive"))
  return(kern)
}

sigmoidWhiteKernExtractParam <- function (kern, only.values = TRUE, 
                                      untransformed.values = TRUE) {
  params <- c(kern$a, kern$b, kern$c)
  
  if (!only.values) 
    names(params) <- c("a", "b", "c")
  return(params)
}

sigmoidWhiteKernExpandParam <- function (kern, params)  {
  if (is.list(params)) 
    params <- params$values
  kern$a <- params[1]
  kern$b <- params[2]
  kern$c <- params[3]
  return(kern)
}

sigmoidWhiteKernCompute <- function (kern, x, x2 = NULL) {

  if (nargs() < 3) {
    n2 <- .prod2(x, x) 
    n11 = matrix(x^2, length(x), length(x))
    n22 = t(n11)
    white = kern$c*diag(length(x))
  }
  else {
    n2 <- .prod2(x, x2)
    n11 = matrix(x^2, length(x), length(x2))
    n22 = t(matrix(x2^2, length(x2), length(x)))
    white = 0
  }
  
  k <- asin((kern$a + kern$b*n2)/sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))) 
         + white

  return(k)
}

sigmoidWhiteKernGradient <- function (kern, x, x2, covGrad) 
{
  
  n11 = matrix(x^2, length(x), length(x))

  if (nargs() == 3) {
    n2 <- .prod2(x, x) 
    n22 = t(n11)
    covGrad <- x2
    white.grad = sum(diag(as.matrix(covGrad)))    
  }
  else if (nargs() == 4) {
    n2 <- .prod2(x, x2) 
    n22 = t(matrix(x2^2, length(x2), length(x2)))
    white.grad = 0
  }
  
  num = kern$a + kern$b*n2
  denom = sqrt((kern$a + 1 + kern$b*n11)*(kern$a + 1 + kern$b*n22))
  ratio = num/denom
  
  #ratio[ratio==1] = 0.99999999
  
  deriv.asin = 1 / sqrt(1 - (ratio)^2)
  deriv.denom = (denom)^3
  
  g <- array()
  
  g[1] <- sum(covGrad * deriv.asin * 
    (1/denom - 0.5*num*(2*kern$a + 2 + kern$b*n11 + kern$b*n22)/deriv.denom))
  
  g[2] <- sum(covGrad * deriv.asin * 
    (n2/denom - 0.5*num*((kern$a + kern$b*n11 + 1)*n22 + 
                         (kern$a + kern$b*n22 + 1)*n11)/deriv.denom))
 
  g[3] <- white.grad 

  if (any(is.nan(g))) {
    warning("g is NaN.")
    browser()
  }
  
  return(g)
}

sigmoidWhiteKernDiagCompute <- function (kern, x) 
{
  x2 = x^2
  k <- matrix(asin((kern$a + kern$b*x2)/(kern$a + 1 + kern$b*x2)) 
         + kern$c, dim(as.array(x))[1], 1)

  return(k)
}

.prod2 <- function (x, x2) 
{
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))
  xMat <- array(apply(as.matrix(x), 1, sum), c(xdim[1], 
                                                   x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2), 1, sum), c(x2dim[1], 
                                                        xdim[1])))
  if (xdim[2] != x2dim[2]) 
    stop("Data dimensions are not matched.")
  n2 <- xMat * x2Mat
  return(n2)
}




periodicKernParamInit <- function(kern) {
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$invperiod <- 3#1/(2*pi)
  kern$nParams <- 3
  kern$paramNames <- c("inverseWidth", "variance", "invperiod")
  kern$isStationary <- TRUE
  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && 
    kern$options$isNormalised) 
    kern$isNormalised <- TRUE
  else kern$isNormalised <- FALSE
  if ("options" %in% names(kern) && "inverseWidthBounds" %in% 
    names(kern$options)) {
    kern$transforms <- list(list(index = 1, type = "bounded"), 
                            list(index = 2, type = "positive"),
                            list(index = 3, type = "bounded"))
    kern$transformArgs <- list()
    kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
    kern$transformArgs[[2]] <- kern$options$inversePeriodWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  }
  else {
    kern$transforms <- list(list(index = c(1, 2, 3), type = "positive"))
  }
  return(kern)
}

periodicKernExtractParam <- function (kern, only.values = TRUE, 
                                      untransformed.values = TRUE) {
  params <- c(kern$inverseWidth, kern$variance, kern$invperiod)
  #params <- c(kern$inverseWidth, kern$variance)
  if (!only.values) 
    names(params) <- c("inverseWidth", "variance", "invperiod")
    #names(params) <- c("inverseWidth", "variance")
  return(params)
}

periodicKernExpandParam <- function (kern, params)  {
  if (is.list(params)) 
    params <- params$values
  kern$inverseWidth <- params[1]
  kern$variance <- params[2]
  kern$invperiod <- params[3]
  return(kern)
}

periodicKernCompute <- function (kern, x, x2 = NULL) {
  if (nargs() < 3) {
    n1 <- .dist1(x, x)
  }
  else {
    n1 <- .dist1(x, x2)
  }
  wi2 <- 2 * kern$inverseWidth
  k <- kern$variance * exp(- sin(pi * n1 / kern$invperiod)^2 * wi2)
  
  #if ("isNormalised" %in% names(kern) && kern$isNormalised) 
  #  k <- k * sqrt(kern$inverseWidth/(2 * pi))
  
  return(k)
}

periodicKernGradient <- function (kern, x, x2, covGrad) 
{
  if (nargs() == 3) {
    k <- periodicKernCompute(kern, x)
    dist1xx <- pi * .dist1(x, x)
    covGrad <- x2
  }
  else if (nargs() == 4) {
    k <- periodicKernCompute(kern, x, x2)
    dist1xx <- pi * .dist1(x, x2)
  }
  g <- array()
  
  g[1] <- - 2 * sum(covGrad * k * sin(dist1xx / kern$invperiod)^2)
  
  g[2] <- sum(covGrad * k)/kern$variance

  #g[3] <-  - 2 * kern$inverseWidth * sum(dist1xx * covGrad * k * 
  #  sin(2* kern$invperiod * dist1xx))
  g[3] <-  2 * kern$inverseWidth * sum(dist1xx * covGrad * k * 
      sin(2 * dist1xx / kern$invperiod)) / kern$invperiod^2
    

  if (any(is.nan(g))) 
    warning("g is NaN.")
  
  return(g)
}

.dist1 <- function (x, x2) 
{
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))
  xMat <- array(apply(as.matrix(x), 1, sum), c(xdim[1], 
                                                   x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2), 1, sum), c(x2dim[1], 
                                                        xdim[1])))
  if (xdim[2] != x2dim[2]) 
    stop("Data dimensions are not matched.")
  n1 <- xMat - x2Mat
  return(n1)
}

periodicKernDiagCompute <- function (kern, x) 
{
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)
  return(k)
}


