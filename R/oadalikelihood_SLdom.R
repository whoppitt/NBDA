#Editted for constrained model

oadaLikelihood_SLdom <- function(parVect, nbdadata){
 
if(is.character(nbdadata)){
	
		totalLikelihood <- 0;
			
		for(i in 1:length(nbdadata)){
			subdata <- eval(as.name(nbdadata[i]));
			totalLikelihood <- totalLikelihood+ oadaLikelihood_SLdom(parVect= parVect, nbdadata=subdata);
			}
					
		return(totalLikelihood);
					
}else{

	#Define required function
	sumWithoutNA <- function(x) sum(na.omit(x))

	#calculate the number of each type of parameter
	noSParam <- dim(nbdadata@stMetric)[2] -1#s parameters
	#MInus 1 to account for the fixed reference, s1
	noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
	noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
	noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects
	
	if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
	if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0
	
	datalength <- length(nbdadata@id) #ILV effects on social transmission
	
	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion
	
	#assign different paramreter values to the right vectors
	asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
	intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]
	multiCoef<-parVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]
	if(noSParam==0){sParam <- 1}else{sParam <- c(1,parVect[1:noSParam])}
	#Extra 1 added to sParam for the reference s1 level
	
		
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

	sParam.mat <- matrix(data=rep(sParam+1, datalength), nrow=datalength, byrow=T) # create a matrix of sParams 
	unscaled.st <- apply(sParam.mat*nbdadata@stMetric, MARGIN=1, FUN=sum)
	unscaled.st<-unscaled.st+nbdadata@offsetCorrection[,1]
	
	#The totalRate is set to zero for naive individuals not in the diffusion for a given event
	asocialRate <- exp(asocialLP)* presentInDiffusion
	socialRate<- exp(socialLP)*unscaled.st*presentInDiffusion
	
	#Assuming social transmission is dominant, i.e. individuals with non-zero connections always learn before individuals with 0 connections
	
	solverTotalRate<-asocialRate[nbdadata@status==1]*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)+socialRate[nbdadata@status==1]*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)>0)
	
	#Take logs and add across acquisition events
	lComp1 <- sum(log(solverTotalRate)) # group by skilled
	
	lComp2.1 <-tapply(asocialRate, INDEX=nbdadata@event.id, FUN=sum)*(tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)==0)+tapply(socialRate, INDEX=nbdadata@event.id, FUN=sum)
	lComp2.2 <- sum(log(lComp2.1))
	
	negloglik <- lComp2.2 - lComp1
	
	return(negloglik)
	}
}


