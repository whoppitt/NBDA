#CANNOT GET THIS TO WORK WITH RANDOME EFFECTS
#TRIED EXTRACTING THE LP FROM THE COXME AND USING THAT BUT CANNOT GET IT TO WORK
#ALSO TRIED EXTRACTING THE RE FROM THE COXME BUT CANNOT GET THAT TO WORK EITHER

gradient_coxme <- function(parVect, nbdadata,coxmeModel=NULL,coxmeData=NULL){

if(is.character(nbdadata)){
	
		totalGradient <- rep(0, length(parVect));
    
		subdata <- nbdadatatemp<-eval(as.name(nbdadata[1]));
		coxmeData<-createCoxmeData(parVect,subdata)
		if (length(nbdadata)>1){
		  for(i in 2:length(nbdadata)){
		    subdata <- eval(as.name(nbdadata[i]));
		    coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata))
		  }
		}  

		coxmeModel<-return_coxme(parVect=parVect,nbdadata =nbdadata)
							
		for(i in 1:length(nbdadata)){
			subdata <- eval(as.name(nbdadata[i]));

			totalGradient <- totalGradient + gradient_coxme(parVect= parVect, nbdadata=subdata, coxmeModel=coxmeModel,coxmeData=coxmeData);
			}
					
		return(totalGradient);
					
}else{


	#calculate the number of each type of parameter
	noSParam <- dim(nbdadata@stMetric)[2] #s parameters
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
	  sParam <- parVect[1:noSParam]
	  asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
	  intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]
	  
	  #Now get the coxme model and data so the values for the multiILVs can be input, and the LP including REs extracted
	  
	  if(is.null(coxmeModel)){
	    coxmeModel<-return_coxme(parVect=parVect,nbdadata =nbdadata)
	  }
	  if(is.null(coxmeData)){
	    coxmeData<-createCoxmeData(parVect=parVect,nbdadata =nbdadata)
	  }

	  if(nbdadata@asoc_ilv[1]=="ILVabsent") asocialCoef<-NULL
	  if(nbdadata@int_ilv[1]=="ILVabsent") intCoef<-NULL
	  if(nbdadata@multi_ilv[1]=="ILVabsent") {multiCoef<-NULL}else{multiCoef<-coxmeModel$coefficients}
	  
	  
	  # create a matrix of the coefficients to multiply by the observed data values, only if there are asocial variables 
	  if(nbdadata@asoc_ilv[1]=="ILVabsent"){
	    asocialLP<-rep(0,datalength)
	  }else{
	    asocialCoef.mat <- matrix(data=rep(asocialCoef, datalength), nrow=datalength, byrow=T)
	    asocial.sub <- nbdadata@asocILVdata
	    asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
	  }
	  asocialLP<-asocialLP+nbdadata@offsetCorrection[,2]+coxmeModel$frail$Ind[nbdadata@randomEffectdata[,1]]

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
	  multiLP<-multiLP+nbdadata@offsetCorrection[,4]+coxmeModel$frail$Ind[nbdadata@randomEffectdata[,1]]

	  asocialLP<-asocialLP+multiLP
	  socialLP<-socialLP+multiLP

	# create a matrix of s parameters
	sParam.mat <- matrix(data=rep(parVect[1:noSParam],datalength), nrow=datalength, byrow=T) 
	# multiply the matrix of s parameters, by the matrix of observed strength of associations (stMetric), and sum the rows of the resulting matrix to get the unscaled strength of association data
	unscaled.st <- apply(sParam.mat*nbdadata@stMetric, MARGIN=1, FUN=sum)
	unscaled.st<-unscaled.st+nbdadata@offsetCorrection[,1]
		
	# calculate the total rate of learning (of naive individuals) by taking the exponentials of the linear predictors, and multiplying the socialLP by the unscaled association data
	#Individuals not present in the diffusion have their rate set to zero
	totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)*presentInDiffusion
 # totalRate<-exp(overallLP)

# gradient for any s parameter - make sure you apply it to the relevant column of stMetric matrix

#### S PARAMETERS

s_grad <- vector("numeric", length=noSParam)

	for (s in 1:noSParam){
		
	s_grad[s] <- sum((exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s])/totalRate[nbdadata@status==1] - tapply(exp(socialLP)* presentInDiffusion*nbdadata@stMetric[,s], INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)) 
	# NUM: solver social rate/s
	# DENOM: solver total rate # solver total rate

} # closes s for loop

#### ASOCIAL PARAMETERS

if(nbdadata@asoc_ilv[1]!="ILVabsent"){
	
	asocial_grad <- vector("numeric", length=length(nbdadata@asoc_ilv))
	for (i in 1:length(nbdadata@asoc_ilv)){

# UNCONSTRAINED OR ADDITIVE - first derivative of the likelihood function for asocial variables
		asocial_grad[i] <- sum((nbdadata@asocILVdata[nbdadata@status==1,i]*(exp(asocialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1] -tapply(nbdadata@asocILVdata[,i]*(exp(asocialLP))*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)) 
	# NUM: variable for solver * solver asocial rate / solver total rate
	# DENOM: variable for all individiduals * asocial rate, summed over all acquisition events / total naive rate
	} # closes loop through asocialVar	
} else {asocial_grad <- NULL} # closes if !isn.null(asocialVar)


#NO MULTI ILV GRADIENT REQUIRED SINCE THESE ARE OPTIMIZED WITHIN THE COXME MODEL

#### SOCIAL PARAMETERS

if(nbdadata@int_ilv[1]!="ILVabsent"){
	
	social_grad <- vector("numeric", length=length(nbdadata@int_ilv))
	for (i in 1:length(nbdadata@int_ilv)){
			

		social_grad[i] <- sum((nbdadata@intILVdata[nbdadata@status==1,i]*(unscaled.st[nbdadata@status==1]*exp(socialLP[nbdadata@status==1])))/totalRate[nbdadata@status==1] - tapply(nbdadata@intILVdata[,i]*unscaled.st*(exp(socialLP))*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))  
	# variable for solver * solver social rate / solver total rate 
	# variable for all individiduals * social rate, summed over all acquisition events / total naive rate 
	} # closes loop through social var
} else {social_grad <- NULL} # closes if !is.null(nbdadata@asoc)


gradient <- c(s_grad, asocial_grad, social_grad)
return(-gradient)
}
} # end function






