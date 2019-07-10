#for dTADA model

dTadaLikelihood <- function(parVect, nbdadata,baseline="constant",noHazFunctPars=NULL,hazFunct=NULL,cumHaz=NULL){

if(is.list(nbdadata)){

		totalLikelihood <- 0;

		for(i in 1:length(nbdadata)){
			subdata <- nbdadata[[i]];
			totalLikelihood <- totalLikelihood+ dTadaLikelihood(parVect= parVect, nbdadata=subdata,baseline=baseline,hazFunct=hazFunct,cumHaz=cumHaz,noHazFunctPars=noHazFunctPars);
			}

		return(totalLikelihood);

}else{

  if(class(nbdadata)!="dTADAData"){
    return("dTADALikelihood for use with dTADAData only")
  }

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

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion

	#assign different parameter values to the right vectors

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

	#The totalRate is set to zero for naive individuals not in the diffusion for a given period
	#This is the total relative rate from OADA- the beaseline hazard is dealt with in a separate component of the logLik here
	totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion

	#Plug into provided baseline hazard function to get cumlative hazard at start-end

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

	#For non-solvers in each time period
  lComp1<- sum((totalRate*cumHazDiff)[nbdadata@status==0])
  #For solvers in each time period
  lComp2<- sum(log(1-exp((totalRate*cumHazDiff)[nbdadata@status==1])))


	negloglik <- -lComp1-lComp2

	return(negloglik)
	}
}


