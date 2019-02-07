#Editted for unconstrained model v1.5 Nov 2018

#Modified by Will 11/4/17- simplified to be applied to all asocial variables in the object
#Constrained version dropped on 11/06/18

tadaGradient_fn <- function(parVect, nbdadata,baseline="constant",noHazFunctPars=NULL,hazFunct=NULL,cumHaz=NULL){

# Since this is faster when components need to be calulated using the grad function anyway
if(baseline=="gamma"|baseline=="custom"){
  gradient<-tadaGradient_fnNum(parVect, nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz)
  return(gradient)
}


if(is.character(nbdadata)){

		totalGradient <- rep(0, length(parVect));

		for(i in 1:length(nbdadata)){
			subdata <- eval(as.name(nbdadata[i]));
			totalGradient <- totalGradient + tadaGradient_fn(parVect= parVect, nbdadata=subdata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz);
			}

		return(totalGradient);

}else{


	#calculate the number of each type of parameter
	noSParam <- dim(nbdadata@stMetric)[2] #s parameters
	noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
	noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
	noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

	datalength <- dim(nbdadata@stMetric)[1]

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion


	if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
	if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
	if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

	  #assign different paramreter values to the right vectors

	if(baseline=="constant") noHazFunctPars<-1
	if(baseline=="gamma") noHazFunctPars<-2
	if(baseline=="weibull") noHazFunctPars<-2

	hazFunctPars<-parVect[1:noHazFunctPars]
	parVect<-parVect[-(1:noHazFunctPars)]


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

	  solveTimes<-nbdadata@TADAtime2[nbdadata@status==1]

	  #Plug into provided baseline hazard function
	  if(baseline=="constant"){
	    solveHazards<-(1/hazFunctPars)+0*solveTimes
	    cumHazards1<-(1/hazFunctPars)*nbdadata@TADAtime1
	    cumHazards2<-(1/hazFunctPars)*nbdadata@TADAtime2
	    cumHazDiff<-cumHazards1-cumHazards2
	    derivCumHaz1<-cbind(-nbdadata@TADAtime1/(hazFunctPars^2))
	    derivCumHaz2<-cbind(-nbdadata@TADAtime2/(hazFunctPars^2))
	    derivLogHaz<-cbind(-1/hazFunctPars+0*solveTimes)
	  }

#	  if(baseline=="gamma"){
#	    rate=1/hazFunctPars[1]
#	    shape=hazFunctPars[2]
#	    solveHazards<-dgamma(solveTimes,shape=shape,rate=rate)/pgamma(solveTimes,shape=shape,rate=rate, lower = FALSE)
#	    cumHazards1<--pgamma(nbdadata@TADAtime1,shape=shape,rate=rate, lower = FALSE, log = TRUE)
#	    cumHazards2<--pgamma(nbdadata@TADAtime2,shape=shape,rate=rate, lower = FALSE, log = TRUE)
#	    cumHazDiff<-cumHazards1-cumHazards2
#	    derivCumHaz1<-derivCumHaz2<-matrix(NA,nrow=datalength,ncol=2)
#	    cumHaz<-function(parVect,time) -pgamma(time,shape=parVect[2],rate=1/parVect[1], lower = FALSE, log = TRUE)
#	    for(i in 1:datalength){
#	      derivCumHaz1[i,]<-grad(cumHaz,hazFunctPars,time=nbdadata@TADAtime1[i])
#	      derivCumHaz2[i,]<-grad(cumHaz,hazFunctPars,time=nbdadata@TADAtime2[i])
#	    }
#	    logHazFunct<-function(parVect,time) log(dgamma(time,shape=parVect[2],rate=1/parVect[1])/pgamma(time,shape=parVect[2],rate=1/parVect[1], lower = FALSE))
#	    derivLogHaz<-matrix(NA,nrow=length(solveTimes),ncol=2)
#	    for(i in 1:length(solveTimes)){
#	      derivLogHaz[i,]<-grad(logHazFunct,hazFunctPars,time=solveTimes[i])
#	    }
	    #	  }
 #Not currently used since these options are sent to the numerical method

	  if(baseline=="weibull"){
	    scale=hazFunctPars[1]
	    shape=hazFunctPars[2]
	    solveHazards<-dweibull(solveTimes,shape=shape,scale=scale)/pweibull(solveTimes,shape=shape,scale=scale, lower = FALSE)
	    cumHazards1<--pweibull(nbdadata@TADAtime1,shape=shape,scale=scale, lower = FALSE, log = TRUE)
	    cumHazards2<--pweibull(nbdadata@TADAtime2,shape=shape,scale=scale, lower = FALSE, log = TRUE)
	    cumHazDiff<-cumHazards1-cumHazards2
	    derivCumHaz1<-derivCumHaz2<-matrix(NA,nrow=datalength,ncol=2)
	    derivCumHaz1[,1]<--(shape/scale)*(nbdadata@TADAtime1/scale)^shape
	    derivCumHaz2[,1]<--(shape/scale)*(nbdadata@TADAtime2/scale)^shape
	    derivCumHaz1[,2]<-((nbdadata@TADAtime1/scale)^shape)*log(nbdadata@TADAtime1/scale)
	    derivCumHaz2[,2]<-((nbdadata@TADAtime2/scale)^shape)*log(nbdadata@TADAtime2/scale)
      #Replace NAs for start time=0 with 0
	    derivCumHaz1[is.na(derivCumHaz1)]<-0
	    derivLogHaz<-matrix(NA,nrow=length(solveTimes),ncol=2)
	    derivLogHaz[,1]<--shape/scale
	    derivLogHaz[,2]<-(1/shape)+log(solveTimes/scale)
	  }
#	  if(baseline=="custom"){
#	    solveHazards<-hazFunct(hazFunctPars,solveTimes)
#	    cumHazards1<-cumHaz(hazFunctPars,nbdadata@TADAtime1)
#	    cumHazards2<-cumHaz(hazFunctPars,nbdadata@TADAtime2)
#	    cumHazDiff<-cumHazards1-cumHazards2
#	    derivCumHaz1<-derivCumHaz2<-matrix(NA,nrow=datalength,ncol=noHazFunctPars)
#	    for(i in 1:datalength){
#	      derivCumHaz1[i,]<-grad(cumHaz,hazFunctPars,time=nbdadata@TADAtime1[i])
#	      derivCumHaz2[i,]<-grad(cumHaz,hazFunctPars,time=nbdadata@TADAtime2[i])
#	    }
#	    logHazFunct<-function(parVect,time) log(hazFunct(parVect,time))
#	    derivHaz<-matrix(NA,nrow=length(solveTimes),ncol=noHazFunctPars)
#	    for(i in 1:length(solveTimes)){
#	      derivHaz[i,]<-grad(hazFunct,hazFunctPars,time=solveTimes[i])
#	    }
#	  }

#### BASELINE FUNCTION PARAMETERS
	  haz_grad <- vector("numeric", length=noHazFunctPars)

	  for (h in 1:noHazFunctPars){

	    haz_grad[h] <- sum((derivCumHaz1[,h]-derivCumHaz2[,h])*totalRate)+sum(derivLogHaz[,h])

	  }

#### S PARAMETERS

s_grad <- vector("numeric", length=noSParam)

	for (s in 1:noSParam){

	s_grad[s] <- sum((exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s])/totalRate[nbdadata@status==1]) +
	                   sum(exp(socialLP)* presentInDiffusion*nbdadata@stMetric[,s]*cumHazDiff)

} # closes s for loop

#### ASOCIAL PARAMETERS

if(nbdadata@asoc_ilv[1]!="ILVabsent"){

	asocial_grad <- vector("numeric", length=length(nbdadata@asoc_ilv))
	for (i in 1:length(nbdadata@asoc_ilv)){

# UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
		asocial_grad[i] <- sum((nbdadata@asocILVdata[nbdadata@status==1,i]*(exp(asocialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1]) +
		                         sum(nbdadata@asocILVdata[,i]*(exp(asocialLP))*presentInDiffusion*cumHazDiff)
			} # closes loop through asocialVar
} else {asocial_grad <- NULL} # closes if !isn.null(asocialVar)

if(nbdadata@multi_ilv[1]!="ILVabsent"){

  multi_grad <- vector("numeric", length=length(nbdadata@multi_ilv))
  for (i in 1:length(nbdadata@multi_ilv)){

    # UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
    multi_grad[i] <- sum(nbdadata@multiILVdata[nbdadata@status==1,i]) +
                           sum(nbdadata@multiILVdata[,i]*(totalRate)*cumHazDiff)

  } # closes loop through asocialVar
} else {multi_grad <- NULL} # closes if !isn.null(asocialVar)

#### SOCIAL PARAMETERS

if(nbdadata@int_ilv[1]!="ILVabsent"){

	social_grad <- vector("numeric", length=length(nbdadata@int_ilv))
	for (i in 1:length(nbdadata@int_ilv)){


		social_grad[i] <- sum((nbdadata@intILVdata[nbdadata@status==1,i]*(unscaled.st[nbdadata@status==1]*exp(socialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1]) +
		                        sum(nbdadata@intILVdata[,i]*unscaled.st*(exp(socialLP))*presentInDiffusion*cumHazDiff)


	} # closes loop through social var
} else {social_grad <- NULL} # closes if !is.null(nbdadata@asoc)


gradient <- c(haz_grad,s_grad, asocial_grad, social_grad, multi_grad)
return(-gradient)
}
} # end function


# For the gamma function- since numerical approximations are used,it is quicker to just numerically approximate the whole thing thus:
tadaGradient_fnNum <- function(parVect, nbdadata,baseline="constant",noHazFunctPars=NULL,hazFunct=NULL,cumHaz=NULL){
  grad(tadaLikelihood,parVect,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz)
}
