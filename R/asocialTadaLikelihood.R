#Editted for constrained model
#No need for sOffset type anymore- since the ILVs are provided separately for additive, interactive and multiplicative effects

asocialTadaLikelihood <- function(parVect, nbdadata, retainInt=NULL,baseline="constant",noHazFunctPars=NULL,hazFunct=NULL,cumHaz=NULL){
  #We need to know whether to remove the interaction variables. This depends on whether an offset is included for any of the s parameters in any of the diffusions.
  #This will be passed on by the model fitting function, but if the function is called independently we need to calculate this here
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

  if(is.list(nbdadata)){

    totalLikelihood <- 0;

    for(i in 1:length(nbdadata)){
      subdata <- nbdadata[[i]];
      totalLikelihood <- totalLikelihood+ asocialTadaLikelihood(parVect= parVect, nbdadata=subdata,retainInt=retainInt,baseline=baseline,hazFunct=hazFunct,cumHaz=cumHaz,noHazFunctPars=noHazFunctPars);
    }

    return(totalLikelihood);

  }else{

  #If the object is a dTADAData object return the likelihood for the discrete time of acquisition diffusion analysis
  if(class(nbdadata)=="dTADAData"){
    return(asocialDisTadaLikelihood(parVect=parVect, nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz))
  }else{

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

	datalength <- dim(nbdadata@stMetric)[1]

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion

	#assign different parameter values to the right vectors


	if(baseline=="constant") noHazFunctPars<-1
	if(baseline=="gamma") noHazFunctPars<-2
	if(baseline=="weibull") noHazFunctPars<-2

	hazFunctPars<-parVect[1:noHazFunctPars]
	parVect<-parVect[-(1:noHazFunctPars)]

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
	  asocial.sub <- nbdadata@asocILVdata[,]
	  asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
	}
	asocialLP<-asocialLP+nbdadata@offsetCorrection[,2]


	# now calculate the multiplicative LP and add to the asocial LP
	if(nbdadata@multi_ilv[1]=="ILVabsent"){
	  multiLP<-rep(0,datalength)
	}else{
	  multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
	  multi.sub <- nbdadata@multiILVdata[,]
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
	  int.sub <- nbdadata@intILVdata[,]
	  socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
	}
  # calculate
	socialLP<-socialLP+nbdadata@offsetCorrection[,3]+multiLP
}else{socialLP<-rep(0,datalength)}


	#The totalRate is set to zero for naive individuals not in the diffusion for a given event
	totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion


	#Take logs and add across acquisition events
	lComp1 <- sum(log(totalRate[nbdadata@status==1])) # group by skilled

	solveTimes<-nbdadata@TADAtime2[nbdadata@status==1]

	#Plug into provided baseline hazard function
	if(baseline=="constant"){
	  solveHazards<-(1/hazFunctPars)+0*solveTimes
	}
	if(baseline=="gamma"){
	  rate=1/hazFunctPars[1]
	  shape=hazFunctPars[2]
	  solveHazards<-dgamma(solveTimes,shape=shape,rate=rate)/pgamma(solveTimes,shape=shape,rate=rate, lower = FALSE)
	}
	if(baseline=="weibull"){
	  scale=hazFunctPars[1]
	  shape=hazFunctPars[2]
	  solveHazards<-dweibull(solveTimes,shape=shape,scale=scale)/pweibull(solveTimes,shape=shape,scale=scale, lower = FALSE)
	}
	if(baseline=="custom"){
	  solveHazards<-hazFunct(hazFunctPars,solveTimes)
	}

	lComp2<-sum(log(solveHazards))
	#log baseline hazards for solvers at times of solving

	#lComp3
	#Summed relative rates x difference in cumumative baseline hazards

	if(baseline=="constant"){
	  cumHazards1<-(1/hazFunctPars)*nbdadata@TADAtime1
	  cumHazards2<-(1/hazFunctPars)*nbdadata@TADAtime2
	  cumHazDiff<-cumHazards1-cumHazards2
	}
	if(baseline=="gamma"){
	  cumHazards1<--pgamma(nbdadata@TADAtime1,shape=shape,rate=rate, lower = FALSE, log = TRUE)
	  cumHazards2<--pgamma(nbdadata@TADAtime2,shape=shape,rate=rate, lower = FALSE, log = TRUE)
	  cumHazDiff<-cumHazards1-cumHazards2
	}
	if(baseline=="weibull"){
	  cumHazards1<--pweibull(nbdadata@TADAtime1,shape=shape,scale=scale, lower = FALSE, log = TRUE)
	  cumHazards2<--pweibull(nbdadata@TADAtime2,shape=shape,scale=scale, lower = FALSE, log = TRUE)
	  cumHazDiff<-cumHazards1-cumHazards2
	}
	if(baseline=="custom"){
	  cumHazards1<-cumHaz(hazFunctPars,nbdadata@TADAtime1)
	  cumHazards2<-cumHaz(hazFunctPars,nbdadata@TADAtime2)
	  cumHazDiff<-cumHazards1-cumHazards2
	}

	lComp3.1 <- tapply(totalRate*cumHazDiff, INDEX=nbdadata@event.id, FUN=sum) # check this works. this is total rate per event across all naive id
	lComp3.2 <- sum(lComp3.1)
	#Need to calculate

	negloglik <- -lComp1-lComp2-lComp3.2

	return(negloglik)
	}
}}


