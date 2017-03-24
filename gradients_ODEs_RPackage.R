# Calculate gradients from ODE system

getODEGradient <- function(params, X, timePoints, auxVars, species) {


	XODE <- X + auxVars$constant

	if (auxVars$originalPositive==TRUE){
		XODE[XODE < 0] <- 0
	}

  if(auxVars$modelName == 'PIF45') {  
    gradient = gradientsPIF45Model(params, XODE, timePoints, auxVars$TOC1)
  } else if(auxVars$modelName == 'LV') {
    gradient = gradientsLVModel(params, XODE, timePoints)
    gradient = gradient[,species,drop=F]
  } else if(auxVars$modelName == 'VV') {
    gradient = gradientsVVModel(params, XODE, timePoints, species)
  } else if(auxVars$modelName == 'RobertaSimple') {
    gradient = gradientsRobertaSimpleModel(params, XODE, timePoints, auxVars, species)
  } else if(auxVars$modelName == 'RobertaSimpleLog') {
    gradient = gradientsRobertaSimpleLogModel(params, XODE, timePoints, auxVars, species)
  } else if(auxVars$modelName == 'RobertaSimpleStand') {
    gradient = gradientsRobertaSimpleStandModel(params, XODE, timePoints, auxVars, species) 
  } else if(auxVars$modelName == 'RobertaSimpleReparam') {
    gradient = gradientsRobertaSimpleReparamModel(params, XODE, timePoints, auxVars, species)
  } else if(auxVars$modelName == "P2011Model") {
    gradient = gradientsP2011Model(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "VGModel1Full") {
    gradient = gradientsVGModel1Full(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "VGModel2Full") {
    gradient = gradientsVGModel2Full(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "VGModel3Full") {
    gradient = gradientsVGModel3Full(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "VGModel4Full") {
    gradient = gradientsVGModel4Full(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "LotkaVolterraModel") {
    gradient = gradientsLotkaVolterraModel(params,XODE,timePoints,auxVars,species)
  } else if(auxVars$modelName == "FitzHughNagumoModel") {
    gradient = gradientsFitzHughNagumoModel(params,XODE,timePoints,auxVars,species)
  }



  
  return(gradient)
}

