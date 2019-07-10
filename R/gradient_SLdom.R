#Editted for unconstrained model v1.5 Nov 2018

#Modified by Will 11/4/17- simplified to be applied to all asocial variables in the object
#Constrained version dropped on 11/06/18

gradient_SLdom <- function(parVect, nbdadata){

if(is.list(nbdadata)){

    totalGradient <- rep(0, length(parVect));

    for(i in 1:length(nbdadata)){
      subdata <- nbdadata[[i]];
      totalGradient <- totalGradient + gradient_SLdom(parVect= parVect, nbdadata=subdata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz);
    }

    return(totalGradient);

  }else{


	#calculate the number of each type of parameter
  noSParam <- dim(nbdadata@stMetric)[2] -1#s parameters
  #MInus 1 to account for the fixed reference, s1
  noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
	noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
	noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

	# extract the length of the data as the sum of naive individuals over all acquisition events
	datalength <- length(nbdadata@id)

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion


	if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
	if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
	if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

	  #assign different paramreter values to the right vectors
	#assign different paramreter values to the right vectors
	if(noSParam==0){sParam <- 1}else{sParam <- c(1,parVect[1:noSParam])}
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

	# create a matrix of s parameters
	sParam.mat <- matrix(data=rep(parVect[1:(noSParam+1)],datalength), nrow=datalength, byrow=T)
	# multiply the matrix of s parameters, by the matrix of observed strength of associations (stMetric), and sum the rows of the resulting matrix to get the unscaled strength of association data
	unscaled.st <- apply(sParam.mat*nbdadata@stMetric, MARGIN=1, FUN=sum)
	unscaled.st<-unscaled.st+nbdadata@offsetCorrection[,1]

	# calculate the total rate of learning (of naive individuals) by taking the exponentials of the linear predictors, and multiplying the socialLP by the unscaled association data
	#Individuals not present in the diffusion have their rate set to zero

	#The totalRate is set to zero for naive individuals not in the diffusion for a given event
	asocialRate <- exp(asocialLP)* presentInDiffusion
	socialRate<- exp(socialLP)*unscaled.st*presentInDiffusion

	#Assuming social transmission is dominant, i.e. individuals with non-zero connections always learn before individuals with 0 connections

	solverTotalRate<-asocialRate[nbdadata@status==1]*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)+socialRate[nbdadata@status==1]*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)

	summedNaiveTotalRate <-tapply(asocialRate, INDEX=nbdadata@event.id, FUN=sum)*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)+tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)


# gradient for any s parameter - make sure you apply it to the relevant column of stMetric matrix

#### S PARAMETERS
if(noSParam==0){s_grad<-NULL}else{
s_grad <- vector("numeric", length=noSParam-1)

	for (s in 2:(noSParam+1)){

	s_grad[(s-1)] <- sum((tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)*((exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s])/solverTotalRate - tapply(exp(socialLP)* presentInDiffusion*nbdadata@stMetric[,s], INDEX=nbdadata@event.id, FUN=sum)/summedNaiveTotalRate))
	# NUM: solver social rate/s
	# DENOM: solver total rate # solver total rate

} # closes s for loop
}

#### ASOCIAL PARAMETERS

if(nbdadata@asoc_ilv[1]!="ILVabsent"){

	asocial_grad <- vector("numeric", length=length(nbdadata@asoc_ilv))
	for (i in 1:length(nbdadata@asoc_ilv)){

# UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
	  temp_grad<-rep(NA,sum(nbdadata@status==1))
	  temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)]<-
		                         ((nbdadata@asocILVdata[nbdadata@status==1,i]*(exp(asocialLP[nbdadata@status==1])))/asocialRate[nbdadata@status==1] -
		                                                                                    tapply(nbdadata@asocILVdata[,i]*(exp(asocialLP))*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/
		                            tapply(asocialRate, INDEX=nbdadata@event.id, FUN=sum))[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)]
	  temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)]<-0

	    asocial_grad[i] <-sum(temp_grad)
	# NUM: variable for solver * solver asocial rate / solver total rate
	# DENOM: variable for all individiduals * asocial rate, summed over all acquisition events / total naive rate
	} # closes loop through asocialVar
} else {asocial_grad <- NULL} # closes if !isn.null(asocialVar)

if(nbdadata@multi_ilv[1]!="ILVabsent"){

  multi_grad <- vector("numeric", length=length(nbdadata@multi_ilv))
  for (i in 1:length(nbdadata@multi_ilv)){

    # UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables

    temp_grad<-rep(NA,sum(nbdadata@status==1))
    temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)]<-(
                           ((nbdadata@multiILVdata[nbdadata@status==1,i]*(exp(asocialLP[nbdadata@status==1])))/asocialRate[nbdadata@status==1] -
                              tapply(nbdadata@multiILVdata[,i]*(exp(asocialLP))*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/
                              tapply(asocialRate, INDEX=nbdadata@event.id, FUN=sum)))[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)]
    temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)]<-(
                           (nbdadata@multiILVdata[nbdadata@status==1,i] -
                                                                                       tapply(nbdadata@multiILVdata[,i]*(socialRate), INDEX=nbdadata@event.id, FUN=sum)/
                                                                                       tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)))[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)]
    multi_grad[i] <-sum(temp_grad)
    # NUM: variable for solver * solver asocial rate / solver total rate
    # DENOM: variable for all individiduals * asocial rate, summed over all acquisition events / total naive rate
  } # closes loop through asocialVar
} else {multi_grad <- NULL} # closes if !isn.null(asocialVar)

#### SOCIAL PARAMETERS

if(nbdadata@int_ilv[1]!="ILVabsent"){

	social_grad <- vector("numeric", length=length(nbdadata@int_ilv))
	for (i in 1:length(nbdadata@int_ilv)){


		social_grad[i] <- sum()

		temp_grad<-rep(NA,sum(nbdadata@status==1))
		temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)]<-0
		temp_grad[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)]<-(
		  nbdadata@intILVdata[nbdadata@status==1,i] -
		    tapply(nbdadata@intILVdata[,i]*unscaled.st*(exp(socialLP))*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/
		    tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum))[(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)]

		social_grad[i] <- sum(temp_grad)

	# variable for solver * solver social rate / solver total rate
	# variable for all individiduals * social rate, summed over all acquisition events / total naive rate
	} # closes loop through social var
} else {social_grad <- NULL} # closes if !is.null(nbdadata@asoc)


gradient <- c(s_grad, asocial_grad, social_grad, multi_grad)
return(-gradient)
}
} # end function






