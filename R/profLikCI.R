#'Plot profile likelihood for a parameter in an NBDA
#'
#'Find the profile likelihood for a specific parameter, or a difference between two parameters.
#'
#'The profile likelihood method for finding (100-X)\% confidence intervals works by finding the set of values for a
#'parameter that would not be rejected in a likelihood ratio test at the X\%  level of significance. This is equivalent to
#'finding the set of values for which the profile likelihood (-log likelihood optimized over all other parameters in the
#'model) is within C units of the -log-likelihood for the model, where C is the critical value for rejection at the X\%
#'level of significance (1.92 for 95\%  confidence intervals). The \code{plotProfLik} function can be used to plot
#'the profile likelihood for a parameter and find the approximate location of the endpoints of the confidence interval
#'after which \code{\link{profLikCI}} can be used to locate the exact endpoints. The \code{plotProfLik} function plots the
#'profile likelihood for the specified parameter with a dotted line at the point of rejection (C) at the X\% significance
#'level. The endpoints of the confidence interval are where the profile likelihood crosses the dotted line. If necessary
#'the user can reduce the range and re-plot to "zoom in" on each endpoint. Note if points are plotted in red then it means
#'the optimization algorithm did not converge when calculating the profile likelihood for that point.
#'
#'@section Warning: This function does not work when trueTies are present in an OADA. Instead use
#'\code{\link{plotProfLikTrueTies}} for the confidence intervals on a parameter, or \code{\link{plotProfLikDiffTrueTies}}
#' for the difference between two parameters.
#'
#'@section Plotting profile likelihood for a difference between two parameters: This can be achieved using the
#'\code{constraintsVect} argument. e.g. if we wish to find the confidence interval for parameter 1 - parameter 2, we
#'specify \code{which=1} and  \code{constraintsVect=c(1,1,2,3,etc.)}. This constrains parameter 1 and 2 to be the same, but adds
#'an offset to parameter 1 using the \code{\link{constrainedNBDAdata}} function. The resulting profile likelihood is for
#'parameter 1 - parameter 2. This can only be done for parameters of the same type i.e. differences must be within the s
#'parameters, asoc_ilv, int_ilv or multi_ilv categories. If the user wishes to find confidence intervals for the difference
#'between two s parameters which is thought to span zero, we advise doing this as a two step process. e.g. find the upper limit
#'for s1-s2, setting range >0, then find the upper limit for s2-s1 setting range >0. This prevents values of s1 or s2<0 being
#'condsidered in the optimization process, which may trigger errors.
#'
#'@seealso \code{\link{profLikCI}}
#'
#'@param which numeric giving the parameter for which the profile likelihood is to be plotted. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param model object of class \code{\link{oadaFit}} or \code{\link{tadaFit}}.
#'@param range numeric vector of length two, providing a range of parameter values within which the profile likelihood is
#'to be plotted.
#'@param constraintsVect optional numerical vector. This only needs to be used if the confidence interval for the difference
#'between two parameters is required (see specific section below).
#'@param resolution numeric giving the number of points to be plotted. The user is advised to start at a low resolution to
#'obtain the approproate range, then increase resolution to identify ranges for endpoints.
#'@param iterations optional numerical giving the maximum number of iterations to be used by the optimization alogorithms.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty. This simply increases the height of the dotted
#'line above the model maximum likelihood by a factor of \code{inflation}.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param retainInt logical, can be used to force the model to retain int_ilvs in an asocial model. This is used internally
#'by other functions when there is an offset on the s parameters, but can be safely ignored by the user.
#'@return A dataframe giving the plotted values, and an indicator of whether the optimization algorithms converged (0) or
#'not (1) when finding the profile likelihood for that point.


plotProfLik<-function(which,model,range,constraintsVect=NULL,resolution=20,iterations=150,inflation=1,conf=0.95,retainInt=NULL){

  whichIn<-which

  nbdadata<-model@nbdadata

  #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
  if(is.character(nbdadata)){
    newNbdaData<-list()
    for(i in 1:length(nbdadata)){
      newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
    }
    nbdadata<-newNbdaData
  }

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]]
  }else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

	type<-fitType<-model@type
	cutoff<-model@loglik+inflation*qchisq(conf,1)/2

	#If it was an asocial model- was the interaction retained?
	if(is.null(retainInt)){
	  if(is.list(nbdadata)){
	    retainInt<-FALSE
	    for (i in 1:length(nbdadata)){
	      nbdadataTemp2<-nbdadata[[i]];
	      if(sum(nbdadataTemp2@offsetCorrection[,1])>0) retainInt<-TRUE
	    }
	  }else{
	    retainInt<-sum(nbdadata@offsetCorrection[,1])>0
	  }
	}

  #If no constraints vector is specified, then just constrain the parameter being varied to 0 (the add appropriate offset later)
	if(is.null(constraintsVect)) {
		constraintsVect<-model@outputPar*0
		constraintsVect[-which]<-1:(length(constraintsVect)-1)
	}

	#If the model has fitType "asocial"
	#We need to add a one for the s parameters so the constrianed NBDA object can be created
	#And the ILV numbers need shifting up one
	#And which needs incrementing by one too
	if(fitType=="asocial"){
	  constraintsVect<-c(rep(1,noSParam),constraintsVect);
	  constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
	  which<-which+noSParam
	}

	#If the user has specified all zeroes for the s parameters, we need to change the fitType to "asocial"
	#And we need to add a one for the first s parameter so the constrianed NBDA object can be created
	#And the ILV numbers need shifting up one
	if(sum(constraintsVect[1:noSParam])==0){
	  fitType<-"asocial";
	  constraintsVect[1]<-1;
	  constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
	  if(which<=noSParam&noILVint>0) retainInt<-TRUE
	}

	#Divide up the constraints matrix into components, then add back together to account for interactions in the asocial model
	constraintsSparam<-constraintsVect[(1:noSParam)]
	if(noILVasoc==0){constraintsAsocILV<-NULL}else{constraintsAsocILV<-constraintsVect[(noSParam+1):(noSParam+noILVasoc)]}
  if(type=="asocial"&!retainInt){
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-rep(0,noILVint)}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVmulti)]}
    if (which>noSParam+noILVasoc) which<-which+noILVint
  }else{
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]}
  }
	constraintsVect<-c(constraintsSparam,constraintsAsocILV,constraintsIntILV,constraintsMultiILV)

	#Initialize offsetVect to 0s
	offsetVect<-constraintsVect*0

	xVals<-seq(range[1],range[2],length=resolution)
	profLik<-convergence<-rep(NA,length(xVals))

for(i in 1:length(xVals)){
  offsetVect<-constraintsVect*0
  offsetVect[which]<-xVals[i]

		#Create the necessary constrained data objects
		if(is.list(nbdadata)){
		  nbdadataTemp<-list()
		  for(dataset in 1:length(nbdadata)){
		    nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
		  }
		}else{
		  nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
		}

    # This blocked out bit is now being done in the oadaFit_coxme function
		#Fit the model
    #modelTemp<-oadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F,coxmeFit=F)

    #if(class(model2_RE)=="oadaFit_coxme"){
      #if the model is a coxme fit, we fit a model without random effects first, and use the MLEs as starting values for the coxme fit
    #  startVals<-modelTemp@outputPar[-((length(modelTemp@outputPar)-noILVmulti+1):length(modelTemp@outputPar))]
    #  modelTemp<-oadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F,startValue = startVals)
    #}

    if(class(model)=="oadaFit"|class(model)=="oadaFit_coxme"){
      noHazFunctPars<-0
      modelTemp<-oadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F)
    }
    if(class(model)=="tadaFit"){
      noHazFunctPars<-model@noHazFunctPars
      modelTemp<-tadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F,baseline=model@baseline,noHazFunctPars=model@noHazFunctPars,hazFunct=model@hazFunct,cumHaz=model@cumHaz)
    }

		if(is.null(modelTemp@optimisation$convergence)){
		  convergence[i]<-0
		}else{
		  convergence[i]<-modelTemp@optimisation$convergence
		}
		profLik[i]<-modelTemp@loglik
		plot(xVals[convergence==0],profLik[convergence==0],type="l",xlim=range,ylim=c(model@loglik-(max(na.omit(profLik))-model@loglik)*0.03,max(na.omit(profLik))),xlab=model@varNames[whichIn+noHazFunctPars],ylab="Profile log-likelihood")
    #Plot unconverged points in red
		points(xVals[convergence==1],profLik[convergence==1],col=2)
		abline(h= cutoff, lty=2)
    }
	converged<-rep(NA,length(convergence))
	converged[convergence==0]<-"Yes"
	converged[convergence==1]<-"No"
	return(data.frame(xVals,profLik,converged))
}

distanceFromCutoff<-function(value,which,model,constraintsVect=NULL,iterations=150,inflation=1,conf=0.95,retainInt=NULL,direction=F){

  if(model@nbdaMultiDiff[1]=="NA"){
    nbdadata<-model@nbdadata
  }else{
    nbdadata<-model@nbdaMultiDiff
    #return("Please provide specify data underlying this multi diffusion model")
  }

  #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
  if(is.character(nbdadata)){
    newNbdaData<-list()
    for(i in 1:length(nbdadata)){
      newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
    }
    nbdadata<-newNbdaData
  }

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]]
  }else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  type<-fitType<-model@type
  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  #If it was an asocial model- was the interaction retained?
  if(is.null(retainInt)){
      if(is.list(nbdadata)){
        retainInt<-FALSE
        for (i in 1:length(nbdadata)){
          nbdadataTemp2<-nbdadata[[i]];
          if(sum(nbdadataTemp2@offsetCorrection[,1])>0) retainInt<-TRUE
        }
      }else retainInt<-sum(nbdadata@offsetCorrection[,1])>0
  }

  #If no constraints vector is specified, then just constrain the parameter being varied to 0 (the add appropriate offset later)
  if(is.null(constraintsVect)) {
    constraintsVect<-model@outputPar*0
    constraintsVect[-which]<-1:(length(constraintsVect)-1)
  }

  #If the model has fitType "asocial"
  #We need to add a one for the s parameters so the constrianed NBDA object can be created
  #And the ILV numbers need shifting up one
  #And which needs incrementing by one too
  if(fitType=="asocial"){
    constraintsVect<-c(rep(1,noSParam),constraintsVect);
    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
    which<-which+noSParam
  }

  #If the user has specified all zeroes for the s parameters, we need to change the fitType to "asocial"
  #And we need to add a one for the first s parameter so the constrianed NBDA object can be created
  #And the ILV numbers need shifting up one
  if(sum(constraintsVect[1:noSParam])==0){
    fitType<-"asocial";
    constraintsVect[1]<-1;
    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
  }

  #Divide up the constraints matrix into components, then add back together to account for interactions in the asocial model
  constraintsSparam<-constraintsVect[(1:noSParam)]
  if(noILVasoc==0){constraintsAsocILV<-NULL}else{constraintsAsocILV<-constraintsVect[(noSParam+1):(noSParam+noILVasoc)]}
  if(type=="asocial"&!retainInt){
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-rep(0,noILVint)}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVmulti)]}
    if (which>noSParam+noILVasoc) which<-which+noILVint
  }else{
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]}
  }
  constraintsVect<-c(constraintsSparam,constraintsAsocILV,constraintsIntILV,constraintsMultiILV)


  #Initialize offsetVect to 0s
  offsetVect<-constraintsVect*0

    #Input specified value
  offsetVect[which]<-value

  #Create the necessary constrained data objects
  if(is.list(nbdadata)){
    nbdadataTemp<-list()
    for(dataset in 1:length(nbdadata)){
      nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
    }
  }else{
     nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
  }

  #Fit the model
  if(class(model)=="oadaFit"|class(model)=="oadaFit_coxme")    modelTemp<-oadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F)
  if(class(model)=="tadaFit")    modelTemp<-tadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F,baseline=model@baseline,noHazFunctPars=model@noHazFunctPars,hazFunct=model@hazFunct,cumHaz=model@cumHaz)

  profLik<-modelTemp@loglik
  if(direction){
    return(cutoff-profLik)
  }else{
  	return(abs(cutoff-profLik))
  }
}


#'Find confidence intervals for parameters in an NBDA
#'
#'Find confidence intervals for a specific parameter, or a difference between two parameters using the profile likelihood
#'method.
#'
#'The profile likelihood method for finding (100-X)\% confidence intervals works by finding the set of values for a
#'parameter that would not be rejected in a likelihood ratio test at the X\%  level of significance. This is equivalent to
#'finding the set of values for which the profile likelihood (-log likelihood optimized over all other parameters in the
#'model) is within C units of the -log-likelihood for the model, where C is the critical value for rejection at the X\%
#'level of significance (1.92 for 95\%  confidence intervals). The \code{\link{plotProfLik}} function can be used to plot
#'the profile likelihood for a parameter and find the approximate location of the endpoints of the confidence interval
#'after which \code{profLikCI} can be used to locate the exact endpoints.
#'
#'@section Warning: This function does not work when trueTies are present in an OADA. Instead use
#'\code{\link{profLikCITrueTies}} for the confidence intervals on a parameter, or \code{\link{profLikCIDiffTrueTies}} for
#'the difference between two parameters.
#'
#'@section Getting confidence intervals for a difference between two parameters: This can be achieved using the
#'\code{constraintsVect} argument. e.g. if we wish to find the confidence interval for parameter 1 - parameter 2, we
#'specify \code{which=1} and  \code{constraintsVect=c(1,1,2,3,etc.)}. This constrains parameter 1 and 2 to be the same, but adds
#'an offset to parameter 1 using the \code{\link{constrainedNBDAdata}} function. The resulting confidence interval is for
#'parameter 1 - parameter 2. This can only be done for parameters of the same type i.e. differences must be within the s
#'parameters, asoc_ilv, int_ilv or multi_ilv categories. If the user wishes to find confidence intervals for the difference
#'between two s parameters which is thought to span zero, we advise doing this as a two step process. e.g. find the upper limit
#'for s1-s2, setting range >0, then find the upper limit for s2-s1 setting range >0. This prevents values of s1 or s2<0 being
#'condsidered in the optimization process, which may trigger errors.
#'
#'@seealso \code{\link{plotProfLik}}
#'
#'@param which numeric giving the parameter for which the confidence interval is to be calculated. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param model object of class \code{\link{oadaFit}} or \code{\link{tadaFit}}
#'@param upperRange numeric vector of length two, providing a range within which the upper endpoint of the confidence interval
#'is known to be located. This range can be indentified using \code{\link{plotProfLik}}. \code{upperRange} can be omitted if the
#'confidence interval only has a lower endpoint- i.e. if the profile likelihood levels out below the dotted line.
#'@param lowerRange numeric vector of length two, providing a range within which the lower endpoint of the confidence interval
#'is known to be located. This range can be indentified using \code{\link{plotProfLik}}. \code{lowerRange} can be omitted if the
#'confidence interval only has a upper endpoint- e.g. if the profile likelihood is beneath the dotted line for s=0.
#'@param constraintsVect optional numeric vector. This only needs to be used if the confidence interval for the difference
#'between two parameters is required (see specific section below).
#'@param iterations optional numerical giving the maximum number of iterations to be used by the optimization alogorithms.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param retainInt logical, can be used to force the model to retain int_ilvs in an asocial model. This is used internally
#'by other functions when there is an offset on the s parameters, but can be safely ignored by the user.
#'@return A list of the form ("Lower CI","Upper CI")


profLikCI<-function(which,model,upperRange=NULL,lowerRange=NULL,constraintsVect=NULL,iterations=150,inflation=1,conf=0.95,retainInt=NULL){

	if(!is.null(upperRange)){
			temp1<-optimise(distanceFromCutoff,upperRange,which=which,model=model, constraintsVect= constraintsVect,iterations=iterations,inflation=inflation,conf=conf,retainInt=retainInt)
			upperPoint<-temp1$minimum
		}else{upperPoint<-NA}

	#Evaluate at 0
  #In order to do this we record the original values of which and constraintsVect to be used with distanceFromCutOff later

  originalWhich<-which;
  originalConstrainstsVect<-constraintsVect;

  if(model@nbdaMultiDiff[1]=="NA"){
    nbdadata<-model@nbdadata
  }else{
    nbdadata<-model@nbdaMultiDiff
    #return("Please provide specify data underlying this multi diffusion model")
  }

  #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
  if(is.character(nbdadata)){
    newNbdaData<-list()
    for(i in 1:length(nbdadata)){
      newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
    }
    nbdadata<-newNbdaData
  }

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]]
  }else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  type<-fitType<-model@type
  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  #If it was an asocial model- was the interaction retained?
  originalRetainInt<-retainInt
  if(is.null(retainInt)){
    if(is.list(nbdadata)){
      retainInt<-FALSE
      for (i in 1:length(nbdadata)){
        nbdadataTemp2<-nbdadata[[i]];
        if(sum(nbdadataTemp2@offsetCorrection[,1])>0) retainInt<-TRUE
      }
    }else{
      retainInt<-sum(nbdadata@offsetCorrection[,1])>0
    }
  }

  #If no constraints vector is specified, then just constrain the parameter being varied to 0 (the add appropriate offset later)
  if(is.null(constraintsVect)) {
    constraintsVect<-model@outputPar*0
    constraintsVect[-which]<-1:(length(constraintsVect)-1)
  }

  #If the model has fitType "asocial"
  #We need to add a one for the s parameters so the constrianed NBDA object can be created
  #And the ILV numbers need shifting up one
  #And which needs incrementing by one too
  if(fitType=="asocial"){
    constraintsVect<-c(rep(1,noSParam),constraintsVect);
    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
    which<-which+noSParam
  }

  #If the user has specified all zeroes for the s parameters, we need to change the fitType to "asocial"
  #And we need to add a one for the first s parameter so the constrianed NBDA object can be created
  #And the ILV numbers need shifting up one
  if(sum(constraintsVect[1:noSParam])==0){
    fitType<-"asocial";
    constraintsVect[1]<-1;
    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
  }

  #Divide up the constraints matrix into components, then add back together to account for interactions in the asocial model
  constraintsSparam<-constraintsVect[(1:noSParam)]
  if(noILVasoc==0){constraintsAsocILV<-NULL}else{constraintsAsocILV<-constraintsVect[(noSParam+1):(noSParam+noILVasoc)]}
  if(type=="asocial"&!retainInt){
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-rep(0,noILVint)}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVmulti)]}
    if (which>noSParam+noILVasoc) which<-which+noILVint
  }else{
    if(noILVint==0){constraintsIntILV<-NULL}else{constraintsIntILV<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]}
    if(noILVmulti==0){constraintsMultiILV<-NULL}else{constraintsMultiILV<-constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]}
  }
  constraintsVect<-c(constraintsSparam,constraintsAsocILV,constraintsIntILV,constraintsMultiILV)

  #Initialize offsetVect to 0s
  offsetVect<-constraintsVect*0
  #Input specified value
  offsetVect[which]<-0

  #Create the necessary constrained data objects
  if(is.list(nbdadata)){
    nbdadataTemp<-list()
    for(dataset in 1:length(nbdadata)){
      nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
    }
  }else{
    nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
  }

  #Fit the model
  if(class(model)=="oadaFit"|class(model)=="oadaFit_coxme")    modelTemp<-oadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F)
  if(class(model)=="tadaFit"){
    noHazFunctPars<-model@noHazFunctPars
    modelTemp<-tadaFit(nbdadata= nbdadataTemp,type=fitType,iterations=iterations,standardErrors=F,baseline=model@baseline,noHazFunctPars=model@noHazFunctPars,hazFunct=model@hazFunct,cumHaz=model@cumHaz)
  }

  profLik<-modelTemp@loglik



if(is.null(lowerRange)){
	if(profLik<cutoff){lowerPoint<-0}else{lowerPoint <-NA}
}else{
  #Put the original values of constraintsVect and which back
  originalRetainInt->retainInt
  originalWhich->which;
  originalConstrainstsVect->constraintsVect;
  temp2<-optimise(distanceFromCutoff, lowerRange,which=which,model=model, constraintsVect= constraintsVect, iterations=iterations,inflation=inflation,conf=conf,retainInt=retainInt)
  lowerPoint<-temp2$minimum
}

	output<-c(lowerPoint,upperPoint)
	names(output)<-list("Lower CI","Upper CI")
	return(output)
}


#The above functions do not work with trueTies since the likelihood cannot be corrected easily when an offsetCorrection is included (I cannot think of any way)
#I wrote the functions below to provide profile likelihood confidence intervals for data with trueTies

#Currently still being tested and not written for asocial models yet

oadaSplitParamsLikelihood<-function(parVect,which,value,nbdadata){
  newParVect<-rep(NA,length(parVect)+1)
  newParVect[-which]<-parVect
  newParVect[which]<-value
  return(oadaLikelihood(newParVect,nbdadata))
}

oadaSplitParamsGradient<-function(parVect,which,value,nbdadata){
  return(grad(func=oadaSplitParamsLikelihood,x=parVect,which=which,value=value,nbdadata=nbdadata))
}

oadaProfileLikelihood<-function(which,value,nbdadata,startValue=NULL,lower=NULL,upper=NULL){

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.character(nbdadata)){
    nbdadataTemp<-eval(as.name(nbdadata[1]));
  }
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]];
  }else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
  if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

  if((noSParam+noILVasoc+noILVint+noILVmulti)==1){
    return(oadaLikelihood(value,nbdadata))
  }else{

  if(which<=noSParam) {noSParam<-noSParam-1}else{noILVasoc<-noILVasoc-1}

  #Set lower values if not specified by the user
  if(is.null(lower)) lower<-c(rep(0,noSParam),rep(-Inf,noILVasoc+noILVint+noILVmulti))

  #Set upper values if not specified by the user
  if(is.null(upper)) upper<-c(rep(Inf,noSParam),rep(Inf,noILVasoc+noILVint+noILVmulti))

  #Set staring values if not specified by the user
  if(is.null(startValue)) startValue<-rep(0,noSParam+noILVasoc+noILVint+noILVmulti);

  model<-nlminb(start=startValue,objective=oadaSplitParamsLikelihood,gradient=oadaSplitParamsGradient,upper=upper,lower=lower,value=value,which=which,nbdadata=nbdadata)

  return(model$objective)
  }
}

#'Plot profile likelihood for a parameter in an OADA with trueTies
#'
#'Find the profile likelihood for a specific parameter in an OADA with trueTies
#'
#'An alternative simple plotting function to plot profile likelihoods, intended for data with trueTies, only works for social
#'models. Functions in much the same way as \code{\link{plotProfLik}}, except it cannot be used to find the profile likelihood
#'for the difference between two parameters. For this use \code{\link{plotProfLikDiffTrueTies}}.
#'
#'@seealso \code{\link{plotProfLik}}
#'
#'@param which numeric giving the parameter for which the profile likelihood is to be plotted. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param model object of class \code{\link{oadaFit}}.
#'@param range numeric vector of length two, providing a range of parameter values within which the profile likelihood is
#'to be plotted.
#'@param resolution numeric giving the number of points to be plotted. The user is advised to start at a low resolution to
#'obtain the approproate range, then increase resolution to identify ranges for endpoints.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty. This simply increases the height of the dotted
#'line above the model maximum likelihood by a factor of \code{inflation}.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param lower optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
#' ILVs.
#'@param upper optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be Inf for all parameters.
#'@return A dataframe giving the plotted values.

plotProfLikTrueTies<-function(which,model,range,resolution=20,inflation=1,conf=0.95,startValue=NULL,lower=NULL,upper=NULL){

  if(model@type=="asocial"){
    print("Function does not yet work for asocial models")
    return(NULL)
  }

  constraintsVect=NULL

  nbdadata<-model@nbdadata

  xVals<-seq(range[1],range[2],length=resolution)
  profLik<-rep(NA,length(xVals))

  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  for(i in 1:length(xVals)){
    profLik[i]<-oadaProfileLikelihood(which=which, value=xVals[i],nbdadata=nbdadata,startValue=startValue,lower=lower,upper=upper)
    plot(xVals,profLik,type="l",xlim=range,ylim=c(model@loglik-(max(na.omit(profLik))-model@loglik)*0.03,max(na.omit(profLik))),xlab=model@varNames[which],ylab="Profile log-likelihood")
    abline(h= cutoff, lty=2)
  }

  return(data.frame(xVals,profLik))

}

distanceFromCutoffTrueTies<-function(value,which,model,inflation=1,conf=0.95,startValue=NULL,lowerIn=NULL,upperIn=NULL){

  lower<-lowerIn
  upper<-upperIn

  if(model@nbdaMultiDiff[1]=="NA"){
    nbdadata<-model@nbdadata
  }else{
    nbdadata<-model@nbdaMultiDiff
    #return("Please provide specify data underlying this multi diffusion model")
  }

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.character(nbdadata)){
    nbdadataTemp<-eval(as.name(nbdadata[1]));
  }
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]];
  }  else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  type<-fitType<-model@type
  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  return(abs(cutoff-oadaProfileLikelihood(which=which,value=value,nbdadata=nbdadata,startValue=startValue,lower=lower,upper=upper)))
}

#'Find confidence intervals for a parameter in an OADA with trueTies
#'
#'Find the confidence intervals for a specific parameter in an OADA with trueTies
#'
#'An alternative to find confidence intervals, intended for data with trueTies, only works for socialmodels. Functions in much
#'the same way as \code{\link{profLikCI}}, except it cannot be used to find the confidence intervalfor the difference between
#'two parameters, and the upper and lower endpoints need to be found separately. To find the confidence intervalfor the
#'difference between two parameters this use \code{\link{plotProfLikDiffTrueTies}}.
#'
#'@seealso \code{\link{profLikCI}},\code{\link{plotProfLikTrueTies}}, \code{\link{plotProfLikDiffTrueTies}}
#'
#'@param which numeric giving the parameter for which the profile likelihood is to be plotted. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param model object of class \code{\link{oadaFit}}.
#'@param interval numeric vector of length two, providing a range within which the required endpoint of the confidence interval
#'is known to be located. This range can be indentified using \code{\link{plotProfLikTrueTies}}.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty. This simply increases the height of the dotted
#'line above the model maximum likelihood by a factor of \code{inflation}.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param startValue optional numeric vector giving start values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model.
#'@param lower optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
#' ILVs.
#'@param upper optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be Inf for all parameters.
#'@return numeric giving the location of the required endpoint.

profLikCITrueTies<-function(which,model,interval,inflation=1,conf=0.95,startValue=NULL,lower=NULL,upper=NULL){
  fit<-optimise(f=distanceFromCutoffTrueTies, interval=interval,which=which, model=model,inflation=inflation,conf=conf,startValue=startValue,lowerIn=lower,upperIn=upper)
  fit$minimum
}


oadaDiffParamsLikelihood<-function(parVect,which,whichBaseline,value,nbdadata){
  newParVect<-rep(NA,length(parVect)+1)
  newParVect[-which]<-parVect
  newParVect[which]<-newParVect[whichBaseline]+value
  return(oadaLikelihood(newParVect,nbdadata))
}

oadaDiffParamsGradient<-function(parVect,which,whichBaseline,value,nbdadata){
  return(grad(func=oadaDiffParamsLikelihood,x=parVect,which=which,whichBaseline=whichBaseline,value=value,nbdadata=nbdadata))
}


oadaDiffProfileLikelihood<-function(which,whichBaseline,value,nbdadata,startValue=NULL,lower=NULL,upper=NULL){

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.character(nbdadata)){
    nbdadataTemp<-eval(as.name(nbdadata[1]));
  }
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]];
  }else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
  if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

  if((noSParam+noILVasoc+noILVint+noILVmulti)==1){
    return(oadaLikelihood(value,nbdadata))
  }else{

    if(which<=noSParam) {noSParam<-noSParam-1}else{noILVasoc<-noILVasoc-1}

    #Set lower values if not specified by the user
    if(is.null(lower)) lower<-c(rep(0,noSParam),rep(-Inf,noILVasoc+noILVint+noILVmulti))

    #Set upper values if not specified by the user
    if(is.null(upper)) upper<-c(rep(Inf,noSParam),rep(Inf,noILVasoc+noILVint+noILVmulti))

    #Set staring values if not specified by the user
    if(is.null(startValue)) startValue<-rep(0,noSParam+noILVasoc+noILVint+noILVmulti);


    model<-nlminb(start=startValue,objective=oadaDiffParamsLikelihood,gradient=oadaDiffParamsGradient,upper=upper,lower=lower,value=value,which=which,whichBaseline=whichBaseline,nbdadata=nbdadata)
    return(model$objective)

  }
}


#'Plot profile likelihood for the difference between two parameters in an OADA with trueTies
#'
#'Find the profile likelihood for a specific parameter in an OADA with trueTies
#'
#'An alternative simple plotting function to plot profile likelihoods, intended for data with trueTies, only works for social
#'models. Functions in much the same way as \code{\link{plotProfLik}}.
#'
#'@seealso \code{\link{plotProfLik}}
#'
#'@param which numeric giving the first parameter for which the profile likelihood is to be plotted. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param whichBaseline numeric giving the second parameter for which the profile likelihood is to be plotted. The profile
#'likelihood is plotted for which - whichBaseline.
#'@param model object of class \code{\link{oadaFit}}.
#'@param range numeric vector of length two, providing a range of parameter values within which the profile likelihood is
#'to be plotted.
#'@param resolution numeric giving the number of points to be plotted. The user is advised to start at a low resolution to
#'obtain the approproate range, then increase resolution to identify ranges for endpoints.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty. This simply increases the height of the dotted
#'line above the model maximum likelihood by a factor of \code{inflation}.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param lower optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
#' ILVs.
#'@param upper optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be Inf for all parameters.
#'@return A dataframe giving the plotted values.


plotProfLikDiffTrueTies<-function(which,whichBaseline,model,range,constraintsVect=NULL,resolution=20,inflation=1,conf=0.95,startValue=NULL,lower=NULL,upper=NULL){

  if(model@type=="asocial"){
    print("Function does not yet work for asocial models")
    return(NULL)
  }

  nbdadata<-model@nbdadata

  xVals<-seq(range[1],range[2],length=resolution)
  profLik<-rep(NA,length(xVals))

  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  for(i in 1:length(xVals)){
    profLik[i]<-oadaDiffProfileLikelihood(which=which, whichBaseline=whichBaseline,value=xVals[i],nbdadata=nbdadata,startValue=startValue,lower=lower,upper=upper)
    plot(xVals,profLik,type="l",xlim=range,ylim=c(model@loglik-(max(na.omit(profLik))-model@loglik)*0.03,max(na.omit(profLik))),xlab=model@varNames[which],ylab="Profile log-likelihood")
    abline(h= cutoff, lty=2)
  }

  return(data.frame(xVals,profLik))

}

distanceFromCutoffDiffTrueTies<-function(value,which,whichBaseline,model,inflation=1,conf=0.95,startValue=NULL,lowerIn=NULL,upperIn=NULL){

  lower<-lowerIn
  upper<-upperIn

  if(model@nbdaMultiDiff[1]=="NA"){
    nbdadata<-model@nbdadata
  }else{
    nbdadata<-model@nbdaMultiDiff
    #return("Please provide specify data underlying this multi diffusion model")
  }

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.character(nbdadata)){
    nbdadataTemp<-eval(as.name(nbdadata[1]));
  }
  if(is.list(nbdadata)){
    nbdadataTemp<-nbdadata[[1]];
  }  else{nbdadataTemp<-nbdadata}

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  type<-fitType<-model@type
  cutoff<-model@loglik+inflation*qchisq(conf,1)/2

  return(abs(cutoff-oadaDiffProfileLikelihood(which=which,whichBaseline=whichBaseline,value=value,nbdadata=nbdadata,startValue=startValue,lower=lower,upper=upper)))
}

#'Find confidence intervals for the difference between two parameters in an OADA with trueTies
#'
#'Find the confidence intervals for the difference between two parameters in an OADA with trueTies
#'
#'An alternative to find confidence intervals, intended for data with trueTies, only works for socialmodels. Functions in much
#'the same way as \code{\link{profLikCI}}, except the upper and lower endpoints need to be found separately.
#'
#'@seealso \code{\link{profLikCI}},\code{\link{plotProfLikTrueTies}}, \code{\link{plotProfLikDiffTrueTies}}
#'
#'@param which numeric giving the first parameter for the difference. The appropriate number
#'can be identified from the fitted model, by entering \code{<modelName>@@varNames} to extract the variable names from
#'the model. Each variable name is preceded by its number.
#'@param whichBaseline numeric giving the second parameter for the difference. The confidence interval is for which -
#'whichBaseline.
#'@param model object of class \code{\link{oadaFit}}.
#'@param interval numeric vector of length two, providing a range within which the required endpoint of the confidence interval
#'is known to be located. This range can be indentified using \code{\link{plotProfLikTrueTies}}.
#'@param inflation numerical to be used if the confidence intervals are to be inflated by a specified amount, as suggested
#'by Burnham & Anderson (2000) to allow for model selection uncertainty. This simply increases the height of the dotted
#'line above the model maximum likelihood by a factor of \code{inflation}.
#'@param conf numerical giving the level of confidence required, defaulting to the traditional 0.95.
#'@param startValue optional numeric vector giving start values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model.
#'@param lower optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
#' ILVs.
#'@param upper optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be Inf for all parameters.
#'@return numeric giving the location of the required endpoint.

profLikCIDiffTrueTies<-function(which,whichBaseline,model,interval,inflation=1,conf=0.95,startValue=NULL,lower=NULL,upper=NULL){
  fit<-optimise(f=distanceFromCutoffDiffTrueTies, interval=interval,which=which, whichBaseline=whichBaseline,model=model,inflation=inflation,conf=conf,startValue=startValue,lowerIn=lower,upperIn=upper)
  fit$minimum
}
