
#This version provides a table giving the predicted probability that each event occured by social transmission via each network
#The user can just provide the names of a model- in which case the MLEs of that model are input as parameter values
#Instead the user can input their own parameter values, e.g. the upper or lower confidence intervals for a specific parameter, and the MLEs for the other parameters. In which case nbdadata must also be provided
#The function works for multiple diffusions

#Now works for both TADA and OADA, but the first function is kept as a duplicate function for backwards compatibility
oadaPropSolveByST.byevent<-function(par=NULL,nbdadata=NULL,model=NULL,type="social",retainInt=NULL){
  return (nbdaPropSolveByST.byevent(par=par,nbdadata=nbdadata,model=model,type=type,retainInt=retainInt))
}

#'Estimate the probability that each event occured by social transmission
#'
#'Generates a table giving the predicted probability that each event occured by social transmission via each network.
#'
#'The user can provide the name of a model- in which case the MLEs of that model are input as parameter values. Instead the
#'user can input their own parameter values in which case \code{nbdadata} must also be provided. The function works for
#'multiple diffusions. A column is also provided for P(S offset), which gives the probability that transmission was via
#'network effects that have been allocated to the offset fo s parameters, see \code{\link{constrainedNBDAdata}}.
#'
#'@seealso \code{\link{nbdaPropSolveByST}}
#'
#'@param par optional numeric giving the parameter values. Not necessary if \code{model} is specified.
#'@param nbdadata object of class \code{\link{nbdaData}}, list of objects of class \code{\link{nbdaData}}, object of class
#'\code{\link{dTADAData}} or list of objects of class \code{\link{dTADAData}}. Not necessary if \code{model} is specified.
#'@param model object of class \code{\link{oadaFit}} or \code{\link{tadaFit}}.
#'@param type string determining whether the estimates are for an "asocial" or "social" model. This is taken from
#'\code{model} if provided.
#'@param retainInt logical, can be used to force the model to retain int_ilvs in an asocial model. This is used internally
#'by other functions when there is an offset on the s parameters, but can usually be safely ignored by the user.
#'@return A dataframe giving the predicted probability that each event occured by social transmission via each network.
#'@aliases oadaPropSolveByST.byevent


nbdaPropSolveByST.byevent<-function(par=NULL,nbdadata=NULL,model=NULL,type="social",retainInt=NULL){
  if(is.null(model)){
    if(is.null(par)|is.null(nbdadata)){
      return("Please provide a model or input parameter values with data")
    }
#    if(is.list(nbdadata)){
#      nbdaMultiDiff<-nbdadata;
#      nbdadata<-eval(as.name(nbdaMultiDiff[1]));
#    }else{
#      nbdaMultiDiff<-"NA";
#    }
  }else{
      par<-model@outputPar;
      nbdadata<-model@nbdadata;
      type<-model@type

      if(class(model)=="tadaFit"){

        if(model@baseline=="constant") noHazFunctPars<-1
        if(model@baseline=="gamma") noHazFunctPars<-2
        if(model@baseline=="weibull") noHazFunctPars<-2

        par<-par[-(1:noHazFunctPars)]
      }
  }
  parVect<-par

  if(type=="asocial"){
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
  }

  if(is.list(nbdadata)){
    outputMatrix<-NULL
    for(i in 1:length(nbdadata)){
      subdata <- nbdadata[[i]];
      outputMatrix<-rbind(outputMatrix,oadaPropSolveByST.byevent(par=par, nbdadata=subdata,type=type,retainInt=retainInt));
    }

    }else{

    if(type=="social"){

    #Define required function
    sumWithoutNA <- function(x) sum(na.omit(x))

    #calculate the number of each type of parameter
    noSParam <- dim(nbdadata@stMetric)[2] #s parameters
    noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
    noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
    noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

    if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
    if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
    if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

    datalength <- dim(nbdadata@stMetric)[1]

    #assign different parameter values to the right vectors



    #Extract vector giving which naive individuals were present in the diffusion for each acqusition event
    presentInDiffusion<-nbdadata@ presentInDiffusion

    #assign different paramreter values to the right vectors
    sParam <- parVect[1:noSParam]
    asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
    intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]
    multiCoef<-parVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]

    if(nbdadata@asoc_ilv[1]=="ILVabsent") asocialCoef<-NULL
    if(nbdadata@int_ilv[1]=="ILVabsent") intCoef<-NULL
    if(nbdadata@multi_ilv[1]=="ILVabsent") multiCoef<-NULL

    # create a matrix of the coefficients to multiply by the observed data values, only if there are asocial variables
    if(nbdadata@asoc_ilv[1]=="ILVabsent"){
      asocialLP<-rep(0,datalength)
    }else{
      asocialCoef.mat <- matrix(data=rep(asocialCoef, datalength), nrow=datalength, byrow=T)
      asocial.sub <- nbdadata@asocILVdata
      asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
    }
    asocialLP<-asocialLP+nbdadata@offsetCorrection[,2]

    # now do the same for the interaction variables
    if(nbdadata@int_ilv[1]=="ILVabsent"){
      socialLP<-rep(0,datalength)
    }else{
      intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
      int.sub <- nbdadata@intILVdata
      socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
    }
    socialLP<-socialLP+nbdadata@offsetCorrection[,3]

    # now adjust both LPs for the variables specified to have a multiplicative effect (the same effect on asocial and social learning)
    if(nbdadata@multi_ilv[1]=="ILVabsent"){
      multiLP<-rep(0,datalength)
    }else{
      multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
      multi.sub <- nbdadata@multiILVdata
      multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
    }
    multiLP<-multiLP+nbdadata@offsetCorrection[,4]
    asocialLP<-asocialLP+multiLP
    socialLP<-socialLP+multiLP

    sParam.mat <- matrix(data=rep(sParam, datalength), nrow=datalength, byrow=T) # create a matrix of sParams
    unscaled.st <- apply(sParam.mat*nbdadata@stMetric, MARGIN=1, FUN=sum)
    unscaled.st<-unscaled.st+nbdadata@offsetCorrection[,1]

    #The totalRate is set to zero for naive individuals not in the diffusion for a given event
    totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion
    socialRate<-(exp(socialLP)*unscaled.st)* presentInDiffusion

    socialRatesPerEvent<-exp(socialLP)[nbdadata@status==1]*t(par[1:noSParam]*t(nbdadata@stMetric[nbdadata@status==1,]));
    #Quantify social learning in the offset
    socialRateInOffset<-socialRate[nbdadata@status==1]-apply(socialRatesPerEvent,1,sum)
    socialRatesPerEvent<-cbind(socialRatesPerEvent,socialRateInOffset)
    totalRatesPerEvent<-totalRate[nbdadata@status==1];
    outputMatrix<-data.frame(socialRatesPerEvent/totalRatesPerEvent)
    names(outputMatrix)<-c(paste("P(Network ",1:noSParam,")",sep=""),"P(S offset)")
    outputMatrix<-cbind(eventID=nbdadata@event.id[nbdadata@status==1],outputMatrix);

  #end social type code
    }else{
    if(type=="asocial"){

      #Define required function
      sumWithoutNA <- function(x) sum(na.omit(x))

      #calculate the number of each type of parameter
      noSParam <- dim(nbdadata@stMetric)[2] #s parameters
      noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
      noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects
      noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on social learning


      if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
      if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
      if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

      datalength <- length(nbdadata@id) #ILV effects on social transmission

      #Extract vector giving which naive individuals were present in the diffusion for each acqusition event
      presentInDiffusion<-nbdadata@ presentInDiffusion

      #Extend par to include 0 s parameters
      parVect<-c(rep(0,noSParam),parVect)

      #Allow for the fact that the user might provide offsets to the s parameters which might need to be accounted for
      if(retainInt){
        noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
        if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
        if(length(parVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
          cat("Error: parVect wrong length. \nNote a non-zero offset is provided for the s parameters. \nparVect must include values for the interaction effects\n")
          return(NA)
        }
      }else{
        noILVint<-0
        if(length(parVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
          cat("Error: parVect wrong length. \nNote a zero offset is provided for the s parameters. \nparVect must not include values for the interaction effects\n")
          return(NA)
        }
      }

      #assign different paramreter values to the right vectors
      sParam <- parVect[1:noSParam]
      asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
      multiCoef<-parVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]

      if(nbdadata@asoc_ilv[1]=="ILVabsent") asocialCoef<-NULL
      if(nbdadata@int_ilv[1]=="ILVabsent") intCoef<-NULL
      if(nbdadata@multi_ilv[1]=="ILVabsent") multiCoef<-NULL



      # create a matrix of the coefficients to multiply by the observed data values, only if there are asocial variables
      if(nbdadata@asoc_ilv[1]=="ILVabsent"){
        asocialLP<-rep(0,datalength)
      }else{
        asocialCoef.mat <- matrix(data=rep(asocialCoef, datalength), nrow=datalength, byrow=T)
        asocial.sub <- nbdadata@asocILVdata
        asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
      }
      asocialLP<-asocialLP+nbdadata@offsetCorrection[,2]

      # now calculate the multiplicative LP and add to the asocial LP
      if(nbdadata@multi_ilv[1]=="ILVabsent"){
        multiLP<-rep(0,datalength)
      }else{
        multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
        multi.sub <- nbdadata@multiILVdata
        multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
      }
      multiLP<-multiLP+nbdadata@offsetCorrection[,4]
      asocialLP<-asocialLP+multiLP

      unscaled.st<-nbdadata@offsetCorrection[,1]

      #Allow for the fact that the user might provide offsets to the s parameters which might need to be accounted for
      if(retainInt){

        #assign different paramreter values to the right vectors
        intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]

        #interaction variables
        if(nbdadata@int_ilv[1]=="ILVabsent"){
          socialLP<-rep(0,datalength)
        }else{
          intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
          int.sub <- nbdadata@intILVdata
          socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
        }
        # calculate
        socialLP<-socialLP+nbdadata@offsetCorrection[,3]+multiLP
      }else{socialLP<-rep(0,datalength)}

      #The totalRate is set to zero for naive individuals not in the diffusion for a given event
      totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion
      socialRate<-(exp(socialLP)*unscaled.st)* presentInDiffusion

      #Quantify social learning in the offset
      socialRateInOffset<-socialRate[nbdadata@status==1]
      socialRatesPerEvent<-socialRateInOffset
      totalRatesPerEvent<-totalRate[nbdadata@status==1];
      outputMatrix<-data.frame(socialRatesPerEvent/totalRatesPerEvent)
      names(outputMatrix)<-"P(S offset)"
      outputMatrix<-cbind(eventID=nbdadata@event.id[nbdadata@status==1],outputMatrix);


    }else{
      return("Invalid model type, must be social or asocial")
    }
  }
  }
  return(outputMatrix)
}

#This function gives a simplified output giving the predicted proportion of events that occurred via social transmission through each network
#By default we can exclude events that we know must have been innovations (i.e. the first individual to learn in each diffusion)
#The user can specify how many events were innovations, otherwise the code assumes one innovation for every diffusion in which there were no demostrators present at the start

oadaPropSolveByST<-function(par=NULL,nbdadata=NULL,model=NULL,type="social",exclude.innovations=T,innovations=NULL){
 return(nbdaPropSolveByST(par=par,nbdadata=nbdadata,model=model,type=type,exclude.innovations=exclude.innovations,innovations=innovations))
}

#'Estimate the proportion of events that occured by social transmission
#'
#'Generates a table giving the predicted proportion of events that occured by social transmission (propST or %ST when
#'expressed as a \%) via each network, as predicted by the fitted model.
#'
#'The user can provide the name of a model- in which case the MLEs of that model are input as parameter values. Instead the
#'user can input their own parameter values in which case \code{nbdadata} must also be provided. The function works for
#'multiple diffusions. A column is also provided for P(S offset), which gives propST for transmission via
#'network effects that have been allocated to the offset fo s parameters, see \code{\link{constrainedNBDAdata}}.
#'
#'@section Estimating propST corresponding to the upper and lower endpoints of a C.I. for s (single network):
#'Ideally, one needs to optimize the other parameters in the model while fixing s to its upper and lower limits. This can
#'be done as follows:
#' \enumerate{
#' \item Create a new \code{nbdadata} or \code{dTADAData} object using \code{\link{constrainedNBDAdata}} which has the upper
#' limit of s set as an offset for s using the \code{offsetVect} argument.
#' \item Fit an ASOCIAL model to the new constrained data, using \code{\link{oadaFit}} or \code{\link{tadaFit}}.
#' \item Use \code{nbdaPropSolveByST} to estimate propST, noting that the relevant figure will be contained in the
#' P(S offset) column.
#' \item Repeat for the lower limit for s.
#' }
#'
#'@section Estimating propST corresponding to the upper and lower endpoints of a C.I. for s (multiple networks):
#'This can be done as follows:
#' \enumerate{
#' \item Create a new \code{nbdadata} or \code{dTADAData} object using \code{\link{constrainedNBDAdata}} with the target s
#' parameter constrained to s=0 using the \code{constraintsVect} argument AND an offset equal to its upper limit using
#' the \code{offsetVect} argument.
#' \item Fit a SOCIAL model to the new constrained data, using \code{\link{oadaFit}} or \code{\link{tadaFit}}.
#' \item Use \code{nbdaPropSolveByST} to estimate propST, noting that the relevant figure will be contained in the
#' P(S offset) column.
#' \item Repeat for the lower limit for the target s parameter.
#' }
#'
#'@seealso \code{\link{nbdaPropSolveByST.byevent}}
#'
#'@param par optional numeric giving the parameter values. Not necessary if \code{model} is specified.
#'@param nbdadata object of class \code{\link{nbdaData}}, list of objects of class \code{\link{nbdaData}}, object of class
#'\code{\link{dTADAData}} or list of objects of class \code{\link{dTADAData}}. Not necessary if \code{model} is specified.
#'@param model object of class \code{\link{oadaFit}} or \code{\link{tadaFit}}.
#'@param type string determining whether the estimates are for an "asocial" or "social" model. This is taken from
#'\code{model} if provided.
#'@param retainInt logical, can be used to force the model to retain int_ilvs in an asocial model. This is used internally
#'by other functions when there is an offset on the s parameters, but can usually be safely ignored by the user.
#'@param exclude.innovations logical determining whether innovation events (the first individual to learn in each diffusion)
#'should be excluded from the calculation- since we know the innovation events must occur by asocial learning not social
#'transmission.
#'@param innovations numerical giving the number of innovations across all diffusions. By default this is assumed to be one
#'innovator per diffusion in which there were no demonstrators (see \code{demons} argument in \code{\link{nbdaData}})
#'@return A dataframe giving the predicted proportion of events that occured by social transmission via each network.
#'@aliases oadaPropSolveByST


nbdaPropSolveByST<-function(par=NULL,nbdadata=NULL,model=NULL,type="social",exclude.innovations=T,innovations=NULL){
  byEvent<-nbdaPropSolveByST.byevent(par=par, nbdadata=nbdadata, model=model, type=type)
  #The remainder calculates the number of definite innovations if not specified by the user. We assume any diffusion without demonstrators had an innovator included in the learning events
    if(is.null(model)){
    if(is.null(par)|is.null(nbdadata)){
      return("Please provide a model or input parameter values with data")
    }
      if(is.character(nbdadata)){
        newNbdaData<-list()
        for(i in 1:length(nbdadata)){
          newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
        }
        nbdadata<-newNbdaData
      }
      nbdaMultiDiff<-"NA";
    }else{
    par<-model@outputPar;
    nbdadata<-model@nbdadata;
    type<-model@type
  }
  if(is.null(innovations)){
    if(!is.list(nbdadata)){
    innovations<-is.na(nbdadata@demons[1])*1
  }else{
    innovations<-0
    for(i in 1:length(nbdadata)){
      subdata <- nbdadata[[i]];
      innovations<-innovations+is.na(subdata@demons[1])*1
    }
  }
  }
  numbers<-apply(as.matrix(byEvent[,-1]),2,sum)
  total<-dim(byEvent)[1]-innovations*exclude.innovations
  output<-numbers/total
  names(output)<-names(byEvent)[-1]
  #Return rounded to prevent numerical errors from creating non-zero terms for the offset when the offset is zero
  return(round(output,5))
}

