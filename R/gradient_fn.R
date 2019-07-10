#Editted for unconstrained model v1.5 Nov 2018

#Modified by Will 11/4/17- simplified to be applied to all asocial variables in the object
#Constrained version dropped on 11/06/18

gradient_fn <- function(parVect, nbdadata){

  if(is.list(nbdadata)){

    totalGradient <- rep(0, length(parVect));

    for(i in 1:length(nbdadata)){
      subdata <- nbdadata[[i]];
      totalGradient <- totalGradient + gradient_fn(parVect= parVect, nbdadata=subdata);
    }

    return(totalGradient);

  }else{

  if(!is.null(nbdadata@trueTies[[1]])){
    return(grad(oadaLikelihood,parVect, nbdadata=nbdadata))
  }

	#calculate the number of each type of parameter
	noSParam <- dim(nbdadata@stMetric)[2] #s parameters
	noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
	noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
	noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

	includeInOADA<-nbdadata@event.id%in% nbdadata@event.id[nbdadata@status==1]
	#Exclude the lines of data corresponding to the final period to endtime, if necessary
	datalength <- sum(includeInOADA)

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion[includeInOADA]


	if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
	if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
	if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

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
	    asocial.sub <- nbdadata@asocILVdata[includeInOADA,]
	    asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
	  }
	  asocialLP<-asocialLP+nbdadata@offsetCorrection[includeInOADA,2]

	  # now do the same for the interaction variables
	  if(nbdadata@int_ilv[1]=="ILVabsent"){
	    socialLP<-rep(0,datalength)
	  }else{
	    intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
	    int.sub <- nbdadata@intILVdata[includeInOADA,]
	    socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
	  }
	  socialLP<-socialLP+nbdadata@offsetCorrection[includeInOADA,3]

	  # now adjust both LPs for the variables specified to have a multiplicative effect (the same effect on asocial and social learning)
	  if(nbdadata@multi_ilv[1]=="ILVabsent"){
	    multiLP<-rep(0,datalength)
	  }else{
	    multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
	    multi.sub <- nbdadata@multiILVdata[includeInOADA,]
	    multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
	  }
	  multiLP<-multiLP+nbdadata@offsetCorrection[includeInOADA,4]
	  asocialLP<-asocialLP+multiLP
	  socialLP<-socialLP+multiLP

	  sParam.mat <- matrix(data=rep(sParam, datalength), nrow=datalength, byrow=T) # create a matrix of sParams
	  unscaled.st <- apply(sParam.mat*nbdadata@stMetric[includeInOADA,], MARGIN=1, FUN=sum)
	  unscaled.st<-unscaled.st+nbdadata@offsetCorrection[includeInOADA,1]

	  #The totalRate is set to zero for naive individuals not in the diffusion for a given event
	  totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion



# gradient for any s parameter - make sure you apply it to the relevant column of stMetric matrix

#### S PARAMETERS

s_grad <- vector("numeric", length=noSParam)

	for (s in 1:noSParam){

	s_grad[s] <- sum((exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s])/totalRate[nbdadata@status==1] - tapply(exp(socialLP)* presentInDiffusion*nbdadata@stMetric[includeInOADA,s], INDEX=nbdadata@event.id[includeInOADA], FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id[includeInOADA], FUN=sum))
	# NUM: solver social rate/s
	# DENOM: solver total rate # solver total rate

} # closes s for loop

#### ASOCIAL PARAMETERS

if(nbdadata@asoc_ilv[1]!="ILVabsent"){

	asocial_grad <- vector("numeric", length=length(nbdadata@asoc_ilv))
	for (i in 1:length(nbdadata@asoc_ilv)){

# UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
		asocial_grad[i] <- sum((nbdadata@asocILVdata[nbdadata@status==1,i]*(exp(asocialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1] -tapply(nbdadata@asocILVdata[includeInOADA,i]*(exp(asocialLP))*presentInDiffusion, INDEX=nbdadata@event.id[includeInOADA], FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id[includeInOADA], FUN=sum))
	# NUM: variable for solver * solver asocial rate / solver total rate
	# DENOM: variable for all individiduals * asocial rate, summed over all acquisition events / total naive rate
	} # closes loop through asocialVar
} else {asocial_grad <- NULL} # closes if !isn.null(asocialVar)

if(nbdadata@multi_ilv[1]!="ILVabsent"){

  multi_grad <- vector("numeric", length=length(nbdadata@multi_ilv))
  for (i in 1:length(nbdadata@multi_ilv)){

    # UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
    multi_grad[i] <- sum(nbdadata@multiILVdata[nbdadata@status==1,i] - tapply(nbdadata@multiILVdata[includeInOADA,i]*(totalRate), INDEX=nbdadata@event.id[includeInOADA], FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id[includeInOADA], FUN=sum))
    # NUM: variable for solver * solver asocial rate / solver total rate
    # DENOM: variable for all individiduals * asocial rate, summed over all acquisition events / total naive rate
  } # closes loop through asocialVar
} else {multi_grad <- NULL} # closes if !isn.null(asocialVar)

#### SOCIAL PARAMETERS

if(nbdadata@int_ilv[1]!="ILVabsent"){

	social_grad <- vector("numeric", length=length(nbdadata@int_ilv))
	for (i in 1:length(nbdadata@int_ilv)){


		social_grad[i] <- sum((nbdadata@intILVdata[nbdadata@status==1,i]*(unscaled.st[nbdadata@status==1]*exp(socialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1] - tapply(nbdadata@intILVdata[includeInOADA,i]*unscaled.st*(exp(socialLP))*presentInDiffusion, INDEX=nbdadata@event.id[includeInOADA], FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id[includeInOADA], FUN=sum))
	# variable for solver * solver social rate / solver total rate
	# variable for all individiduals * social rate, summed over all acquisition events / total naive rate
	} # closes loop through social var
} else {social_grad <- NULL} # closes if !is.null(nbdadata@asoc)


gradient <- c(s_grad, asocial_grad, social_grad, multi_grad)
return(-gradient)
}
} # end function






