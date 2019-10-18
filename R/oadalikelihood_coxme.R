# Returns the likelihood for an nbda model based around the coxme model. The parameter values are provided
# for s parameters, asocial ILVS and social (int) ILVs, but multi ILV parameters are optimized by the coxme function
# Random effects can be taken into account (and are by default.)


createCoxmeData<-function(parVect,nbdadata,retainInt=TRUE){

  #Define required function
	sumWithoutNA <- function(x) sum(na.omit(x))

	#calculate the number of each type of parameter
	noSParam <- dim(nbdadata@stMetric)[2] #s parameters
	noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
	noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
	noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

	if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
	if(nbdadata@int_ilv[1]=="ILVabsent"|!retainInt) noILVint<-0
  if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

	datalength <- length(nbdadata@id) #ILV effects on social transmission

	#Extract vector giving which naive individuals were present in the diffusion for each acqusition event
	presentInDiffusion<-nbdadata@ presentInDiffusion

	#assign different paramreter values to the right vectors
	sParam <- parVect[1:noSParam]
	asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
	if(noILVint>0) {intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]}else(intCoef<-NULL)

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
	if(nbdadata@int_ilv[1]=="ILVabsent"|!retainInt){
	  socialLP<-rep(0,datalength)
	}else{
	  intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
	  int.sub <- nbdadata@intILVdata
	  socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
	}
	socialLP<-socialLP+nbdadata@offsetCorrection[,3]

	# now add in the multiplicative effect offset(the same effect on asocial and social learning)
	multioffset<-nbdadata@offsetCorrection[,4]

	sParam.mat <- matrix(data=rep(sParam, datalength), nrow=datalength, byrow=T) # create a matrix of sParams
	unscaled.st <- apply(sParam.mat*nbdadata@stMetric, MARGIN=1, FUN=sum)
	unscaled.st<-unscaled.st+nbdadata@offsetCorrection[,1]

	totalOffset<-log((unscaled.st/exp(asocialLP))+exp(-socialLP))+asocialLP+socialLP+multioffset


	coxModelData<-data.frame(time1=nbdadata@time1,time2=nbdadata@time2, status=nbdadata@status,stratum=nbdadata@label,offset=totalOffset,nbdadata@multiILVdata,nbdadata@randomEffectdata)

	return(coxModelData)

}

oadalikelihood_coxme<-function(parVect, nbdadata, formula=NULL){

if(is.character(nbdadata)){
  subdata <- nbdadatatemp<-eval(as.name(nbdadata[1]));
  coxmeData<-createCoxmeData(parVect,subdata)
  if (length(nbdadata)>1){
    for(i in 2:length(nbdadata)){
      subdata <- eval(as.name(nbdadata[i]));
      coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata))
    }
  }
}else{
  if(is.list(nbdadata)){
    subdata <- nbdadatatemp<-nbdadata[[1]];
    coxmeData<-createCoxmeData(parVect,subdata)
    if (length(nbdadata)>1){
      for(i in 2:length(nbdadata)){
        subdata <- nbdadata[[i]];
        coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata))
      }
    }
  }else{
  coxmeData<-createCoxmeData(parVect,nbdadata);
  nbdadatatemp<-nbdadata
}}
  if(is.null(formula)){
    formula<-"Surv(time1,time2,status)~strata(stratum)"
    if(var(coxmeData$offset)>0) {formula<-paste(formula,"+offset(offset)")}
    if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){
      formula<-paste(formula,"+1")
    }else{
      formula<-paste(formula,paste("+",nbdadatatemp@multi_ilv, collapse=""))
    }
#    if(nbdadatatemp@random_effects[1]=="REabsent"){
#      formula<-formula
#    }else{
      formula<-paste(formula,paste("+(1|",nbdadatatemp@random_effects,")", collapse=""))
#    }
    formula<-as.formula(formula)
  }
  model<-coxme(formula=formula,data=coxmeData)
  return(-model$loglik[2])
}

asocialLikelihood_coxme<-function(parVect, nbdadata, retainInt=NULL, formula=NULL){

  if(is.null(retainInt)){
    if(is.character(nbdadata)){
      retainInt<-FALSE
      for (i in 1:length(nbdadata)){
        nbdadatatemp2<-eval(as.name(nbdadata[i]));
        if(sum(nbdadatatemp2@offsetCorrection[,1])>0) retainInt<-TRUE
      }
    }else{
      if(is.list(nbdadata)){
        retainInt<-FALSE
        for (i in 1:length(nbdadata)){
          nbdadatatemp2<-nbdadata[[i]];
          if(sum(nbdadatatemp2@offsetCorrection[,1])>0) retainInt<-TRUE
        }
      }else{
        retainInt<-sum(nbdadata@offsetCorrection[,1])>0
      }
    }
  }


    if(is.list(nbdadata)){
      subdata <- nbdadatatemp<-nbdadata[[1]];
      noSParam<-dim(nbdadatatemp@stMetric)[2] #number of s parameters
      #Append 0s to parVect for the s parameters
      parVect<-c(rep(0,noSParam),parVect)
      coxmeData<-createCoxmeData(parVect,subdata,retainInt=retainInt)
      if (length(nbdadata)>1){
        for(i in 2:length(nbdadata)){
          subdata <- nbdadata[[i]];
          coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata,retainInt=retainInt))
        }
      }
    }else{
    nbdadatatemp<-nbdadata
    noSParam<-dim(nbdadatatemp@stMetric)[2] #number of s parameters
    #Append 0s to parVect for the s parameters
    parVect<-c(rep(0,noSParam),parVect)
    coxmeData<-createCoxmeData(parVect,nbdadata);
    }

  if(is.null(formula)){
    formula<-"Surv(time1,time2,status)~strata(stratum)"
    if(var(coxmeData$offset)>0) {formula<-paste(formula,"+offset(offset)")}
    if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){
      formula<-paste(formula,"+1")
    }else{
      formula<-paste(formula,paste("+",nbdadatatemp@multi_ilv, collapse=""))
    }
#    if(nbdadatatemp@random_effects[1]=="REabsent"){
#      formula<-formula
#    }else{
      formula<-paste(formula,paste("+(1|",nbdadatatemp@random_effects,")", collapse=""))
#    }
    formula<-as.formula(formula)
  }
  model<-coxme(formula=formula,data=coxmeData)
  return(-model$loglik[2])
}
