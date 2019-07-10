#Now corrected for trueTies

asocialLikelihood <- function(parVect, nbdadata, retainInt=NULL){
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
      totalLikelihood <- totalLikelihood + asocialLikelihood(parVect= parVect, nbdadata=subdata,retainInt = retainInt);
    }

    return(totalLikelihood);

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

	includeInOADA<-nbdadata@event.id %in% nbdadata@event.id[nbdadata@status==1]
	#Exclude the lines of data corresponding to the final period to endtime, if necedssary
	datalength <- sum(includeInOADA)

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion[includeInOADA]


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
	  asocial.sub <- nbdadata@asocILVdata[includeInOADA,]
	  asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
	}
	asocialLP<-asocialLP+nbdadata@offsetCorrection[includeInOADA,2]


	# now calculate the multiplicative LP and add to the asocial LP
	if(nbdadata@multi_ilv[1]=="ILVabsent"){
	  multiLP<-rep(0,datalength)
	}else{
	  multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
	  multi.sub <- nbdadata@multiILVdata[includeInOADA,]
	  multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
	}
	multiLP<-multiLP+nbdadata@offsetCorrection[includeInOADA,4]
	asocialLP<-asocialLP+multiLP

	unscaled.st<-nbdadata@offsetCorrection[includeInOADA,1]


#Allow for the fact that the user might provide offsets to the s parameters which might need to be accounted for
if(retainInt){

  #assign different paramreter values to the right vectors
  intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]

	#interaction variables
	if(nbdadata@int_ilv[1]=="ILVabsent"){
	  socialLP<-rep(0,datalength)
	}else{
	  intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
	  int.sub <- nbdadata@intILVdata[includeInOADA,]
	  socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
	}
  # calculate
	socialLP<-socialLP+nbdadata@offsetCorrection[includeInOADA,3]+multiLP
}else{socialLP<-rep(0,datalength)}


	#The totalRate is set to zero for naive individuals not in the diffusion for a given event
	totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion

	#Take logs and add across acquisition events
	lComp1 <- sum(log(totalRate[nbdadata@status==1])) # group by skilled

	lComp2.1 <- tapply(totalRate, INDEX=nbdadata@event.id[includeInOADA], FUN=sum) # check this works. this is total rate per event across all naive id
	lComp2.2 <- sum(log(lComp2.1))

	negloglik <- lComp2.2 - lComp1

	if(!is.null(nbdadata@trueTies[[1]])){
	  negloglik<-negloglik+asocialCorrectTrueTies(parVect[-(1:noSParam)],nbdadata)
	}
	return(negloglik)
	}
}



