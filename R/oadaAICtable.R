
#Define class of object for the fitted additive model
setClass("oadaAICtable",representation(nbdaMultiDiff
="character",nbdadata="list",models='list',convergence="logical",loglik="numeric",aic="numeric",aicc="numeric",constraintsVectMatrix="matrix", offsetVectMatrix="matrix",
MLEs="matrix",SEs="matrix",MLEilv="matrix",SEilv="matrix",MLEint="matrix",SEint="matrix",
typeVect="character",deltaAIC="numeric",RelSupport="numeric",AkaikeWeight="numeric",printTable="data.frame"));


#Method for initializing oadaAICtable object- including model fitting
setMethod("initialize",
    signature(.Object = "oadaAICtable"),
    function (.Object, nbdadata,typeVect,constraintsVectMatrix,offsetVectMatrix,startValue,method,gradient,iterations,aicUse,lowerList,upperList,writeProgressFile,combineTables=F,
              MLEs,SEs,MLEilv,SEilv,MLEint,SEint,
              convergence,loglik,aic,aicc,netComboModifierVect,statusBar,saveModels,stripData,models,...)
    {

    #Save any provided list of models
    modelsIn<-models


    if(is.null(typeVect)){typeVect<-rep("social",dim(constraintsVectMatrix)[1])}

      #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
      if(is.character(nbdadata)){
        newNbdaData<-list()
        for(i in 1:length(nbdadata)){
          newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
        }
        nbdadata<-newNbdaData
      }

      #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
      if(is.list(nbdadata)){
        nbdadataTemp1<-nbdadata[[1]]
      }else{nbdadataTemp1<-nbdadata}


		#if offset matrix is null set it up to contain zeroes
		if(is.null(offsetVectMatrix)) offsetVectMatrix<-constraintsVectMatrix*0

		noModels<-dim(constraintsVectMatrix)[1]

		#Calculate the number of different s parameters, ILVs and models to be fitted
		noSParam<-dim(nbdadataTemp1@stMetric)[2]
		noILVasoc<- dim(nbdadataTemp1@asocILVdata)[2] #ILV effects on asocial learning
		noILVint<- dim(nbdadataTemp1@intILVdata)[2] #ILV effects on interation (social learning)
		noILVmulti<- dim(nbdadataTemp1@multiILVdata)[2] #ILV multiplicative model effects
		if(nbdadataTemp1@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
		if(nbdadataTemp1@int_ilv[1]=="ILVabsent") noILVint<-0
		if(nbdadataTemp1@multi_ilv[1]=="ILVabsent") noILVmulti<-0

		#Record asocialVar names
		asocialVarNames<-unique(c(nbdadataTemp1@asoc_ilv,nbdadataTemp1@int_ilv,nbdadataTemp1@multi_ilv))
		asocialVarNames<-asocialVarNames[asocialVarNames!="ILVabsent"]
		if(is.null(asocialVarNames)){noILVs<-0}else{noILVs<-length(asocialVarNames)}

		#Unless we are combining AICtables already fitted, loop through the constrainstVectMatrix and fit all the models
		if(!combineTables){

		models<-NULL

		#set up progress bar
		if(statusBar) pb <- txtProgressBar(min=0, max=noModels, style=3)

		#Set up matrices to record maximum likelihood estimators and SEs
		MLEs<-matrix(NA,nrow=noModels,ncol=noSParam,dimnames=list(1:noModels, paste("s",1:noSParam,sep="")))
		SEs<-matrix(NA,nrow=noModels,ncol=noSParam,dimnames=list(1:noModels, paste("SEs",1:noSParam,sep="")))
		if(noILVasoc==0){
		  MLEadd<-SEadd<-rep(NA,noModels)
		}else{
		  MLEadd<-matrix(NA,nrow=noModels,ncol= noILVasoc, dimnames=list(1:noModels, nbdadataTemp1@asoc_ilv))
		  SEadd<-matrix(NA,nrow=noModels,ncol= noILVasoc, dimnames=list(1:noModels, nbdadataTemp1@asoc_ilv))
		}
    if(noILVint==0){
      MLEintUC<-SEintUC<-rep(NA,noModels)
    }else{
  		MLEintUC<-matrix(NA,nrow=noModels,ncol= noILVint, dimnames=list(1:noModels, nbdadataTemp1@int_ilv))
  		SEintUC<-matrix(NA,nrow=noModels,ncol= noILVint, dimnames=list(1:noModels,nbdadataTemp1@int_ilv))
    }
		if(noILVmulti==0){
		  MLEmulti<-SEmulti<-rep(NA,noModels)
		}else{
		  MLEmulti<-matrix(NA,nrow=noModels,ncol= noILVmulti, dimnames=list(1:noModels, nbdadataTemp1@multi_ilv))
		  SEmulti<-matrix(NA,nrow=noModels,ncol= noILVmulti, dimnames=list(1:noModels, nbdadataTemp1@multi_ilv))
		}


		#Set up various vectors to record things about each model
		convergence<-rep(NA,noModels)
		loglik<-aic<-aicc<-seApprox<-rep(NaN,noModels)

		#Loop through the rows of the constrainstsVectMatrix creating the constrained objects and thus fitting the specified model each time
		for (i in 1:noModels){

		  #Update progress bar
		  if(statusBar)setTxtProgressBar(pb, i)
		  #Write file to working directory saying what model we are on
      if(writeProgressFile){write.csv(paste("Currently fitting model",i, "out of", noModels),file=paste("oadaTableProgressFile",nbdadataTemp1@label[1],".txt",sep=""),row.names =F)}


		  constraintsVect<-constraintsVectMatrix[i,]
		  offsetVect <-offsetVectMatrix[i,]

		  #If the user has specified all zeroes for the s parameters, we need to change it to an "asocial" type
		  #And we need to add a one for the first s parameter so the constrained NBDA object can be created
		  #And the ILV numbers need shifting up one, to be shifted down later
		  if(sum(constraintsVect[1:noSParam])==0){
		    typeVect[i]<-"asocial";
		    constraintsVect[1]<-1;
		    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
		  }

		  if(is.null(startValue)) {
		    newStartValue<-NULL
		  }else{
		    newStartValue<-startValue[constraintsVect!=0]
		  }

		  if(is.null(lowerList)) {
		    lower<-NULL
		  }else{
		    lower<-lowerList[i,]
		    lower<-lower[constraintsVect!=0]
		  }
		  if(is.null(upperList)) {
		    upper<-NULL
		  }else{
		    upper<-upperList[i,]
		    upper<-upper[constraintsVect!=0]
		  }
		  #Create the necessary constrained data objects
		  if(is.list(nbdadata)){
		    nbdadataTemp<-list()
		    for(dataset in 1:length(nbdadata)){
		      nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
		      }
		  }else{
		    nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
		  }

		  #See if a list of fitted models has been provided and, if so, take the appropriate one
		  if(!is.null(modelsIn)){
		    model<-modelsIn[[i]]
		  }else{
  			#Fit the model
  		  model<-NULL
  			try(model<-oadaFit(nbdadata= nbdadataTemp,type=typeVect[i],lower=lower,upper=upper,startValue=newStartValue,method=method,gradient=gradient,iterations=iterations))
		  }

		  #If the fitted models are to be saved, add this one to the list
  		if(saveModels){
  		  #The user can choose to strip the data from the saved models to save on memory which can be massive with big datasets
  		  if(stripData&!is.null(model)){
  		    model@nbdadata<-list(NA)
  		  }
  		  models<-c(models,list(model))
  		}
      if(!is.null(model)){

			#If it is an asocial model, set constraints to 0 for all s parameters and adjust those for ILVs so they start at 1
			if(typeVect[i]=="asocial"){
			  constraintsVect[1:noSParam]<-0;
			  tempCV<-constraintsVect[-(1:noSParam)]
			  if(max(tempCV)>0) constraintsVect[-(1:noSParam)]<-(tempCV-min(tempCV[tempCV>0])+1)*(tempCV>0)
			}

			#Did the model converge?
			if(is.null(unlist(model@optimisation)[1])){
			  convergence[i]<-T
			}else{
			  if(is.na(unlist(model@optimisation)[1])){convergence[i]<-T}else{convergence[i]<-model@optimisation$convergence==0}
			}

			#Record loglik AIC and AICc
			loglik[i]<-model@loglik
			aic[i]<-model@aic
			aicc[i]<-model@aicc

			#Record MLE and SE for s parameters
			for(j in unique(constraintsVect[1:noSParam])){
			  if(j==0){
			    MLEs[i,constraintsVect[1:noSParam]==j]<-0
			    SEs[i,constraintsVect[1:noSParam]==j]<-0
			  }else{
  			  MLEs[i,constraintsVect[1:noSParam]==j]<-model@outputPar[j]
  			  SEs[i,constraintsVect[1:noSParam]==j]<-model@se[j]
			  }
 			}

			#Record MLE and SE for the  effect of additive ILVs on asocial learning
			if(noILVasoc>0){
			for(j in unique(constraintsVect[(noSParam+1):(noSParam+ noILVasoc)])){
			  if(j==0){
			    MLEadd[i,constraintsVect[(noSParam+1):(noSParam+ noILVasoc)]==j]<-0
			    SEadd[i,constraintsVect[(noSParam+1):(noSParam+ noILVasoc)]==j]<-0
			  }else{
			    MLEadd[i,constraintsVect[(noSParam+1):(noSParam+ noILVasoc)]==j]<-model@outputPar[j]
			    SEadd[i,constraintsVect[(noSParam+1):(noSParam+ noILVasoc)]==j]<-model@se[j]
			  }
			}
			}

			#Record MLE and SE for the  effect of interactive ILVs on social learning
      if(noILVint>0){
			for(j in unique(constraintsVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)])){
			  if(j==0){
			    MLEintUC[i,constraintsVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]==j]<-0
			    SEintUC[i,constraintsVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]==j]<-0
			  }else{
			    MLEintUC[i,constraintsVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]==j]<-model@outputPar[j]
			    SEintUC[i,constraintsVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]==j]<-model@se[j]
			  }
			}
      }

			#Record MLE and SE for the  effect of multiplicative ILVs on social and asocial learning
			if(noILVmulti>0){
			for(j in unique(constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)])){
			 if(j==0){
			    MLEmulti[i,constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]==j]<-0
			    SEmulti[i,constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]==j]<-0
			  }else{
			    MLEmulti[i,constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]==j]<-model@outputPar[j]
			    SEmulti[i,constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]==j]<-model@se[j]
			  }
			}
		}
      }
		}
		if(statusBar)close(pb)

		#We can now sum up the effects on asocial and social learning for each variable
		MLEilv<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("ASOCIAL",asocialVarNames,sep="")))
		SEilv<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SEasocial",asocialVarNames,sep="")))
		MLEint<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SOCIAL",asocialVarNames,sep="")))
		SEint<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SEsocial",asocialVarNames,sep="")))



		for(variable in 1:length(asocialVarNames)){
		  if(sum(unlist(dimnames(MLEadd)[2])==asocialVarNames[variable])>0){
		    MLEilv[,variable]<-MLEilv[,variable]+MLEadd[,unlist(dimnames(MLEadd)[2])==asocialVarNames[variable]]
		    SEilv[,variable]<-SEilv[,variable]+SEadd[,unlist(dimnames(SEadd)[2])==asocialVarNames[variable]]
		  }
		  if(sum(unlist(dimnames(MLEmulti)[2])==asocialVarNames[variable])>0){
		    MLEilv[,variable]<-MLEilv[,variable]+MLEmulti[,unlist(dimnames(MLEmulti)[2])==asocialVarNames[variable]]
		    SEilv[,variable]<-SEilv[,variable]+SEmulti[,unlist(dimnames(SEmulti)[2])==asocialVarNames[variable]]
	      MLEint[,variable]<-MLEint[,variable]+MLEmulti[,unlist(dimnames(MLEmulti)[2])==asocialVarNames[variable]]*(typeVect!="asocial")
		    SEint[,variable]<-SEint[,variable]+SEmulti[,unlist(dimnames(SEmulti)[2])==asocialVarNames[variable]]*(typeVect!="asocial")

		  }
		  if(sum(unlist(dimnames(MLEintUC)[2])==asocialVarNames[variable])>0){
		    MLEint[,variable]<-MLEint[,variable]+MLEintUC[,unlist(dimnames(MLEintUC)[2])==asocialVarNames[variable]]
		    SEint[,variable]<-SEint[,variable]+SEintUC[,unlist(dimnames(SEintUC)[2])==asocialVarNames[variable]]
		  }
		}
		}

		#calculate deltaAIC based on AICc unless user specifies AIC
		if(aicUse=="aic") {deltaAIC<-aic-min(aic)}else{deltaAIC<-aicc-min(aicc)}
		RelSupport<-exp(-0.5*deltaAIC)
		AkaikeWeight<-RelSupport/sum(RelSupport)

		#Give some dimnames to constraints and offset matrices
		varNames<-c(unlist(dimnames(MLEs)[2]))
		if(noILVasoc>0) varNames<-c(varNames,paste("ASOC:",nbdadataTemp1@asoc_ilv, sep=""))
		if(noILVint>0) varNames<-c(varNames,paste("SOCIAL:",nbdadataTemp1@int_ilv, sep=""))
		if(noILVmulti>0) varNames<-c(varNames,paste("A&S:",nbdadataTemp1@multi_ilv,sep=""))

		dimnames(constraintsVectMatrix)=list(1:noModels,paste("CONS:", varNames    ,sep=""))
		dimnames(offsetVectMatrix)=list(1:noModels,paste("OFF", varNames  ,sep=""))

		#Identify which models are asocial models and assign them that type
		typeVect[apply(cbind(constraintsVectMatrix[,1:noSParam]),1,sum)==0]<-"asocial"

		#Classify model types according to ILV effects fitted
		newType<-rep(NA,length(typeVect))

		if(noILVasoc>0&noILVint>0&noILVmulti>0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"additive"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"unconstrained"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"mixed"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"multiplicative"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"socialEffectsOnly"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"noILVs"
		}
		if(noILVasoc==0&noILVint>0&noILVmulti>0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"socialEffectsOnly"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"mixed"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"multiplicative"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"noILVs"
		}
		if(noILVasoc>0&noILVint==0&noILVmulti>0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"additive"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"mixed"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"multiplicative"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"noILVs"
		}
		if(noILVasoc>0&noILVint>0&noILVmulti==0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0]<-"additive"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0]<-"unconstrained"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0]<-"noILVs"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0]<-"socialEffectsOnly"

		}
		if(noILVasoc==0&noILVint==0&noILVmulti>0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)>0]<-"multiplicative"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]),1,sum)==0]<-"noILVs"
		}
		if(noILVasoc>0&noILVint==0&noILVmulti==0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)>0]<-"additive"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+1):(noSParam+noILVasoc)]),1,sum)==0]<-"noILVs"
		}
		if(noILVasoc==0&noILVint>0&noILVmulti==0){
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)>0]<-"socialEffectsOnly"
		  newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)>0&apply(as.matrix(constraintsVectMatrix[,(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]),1,sum)==0]<-"noILVs"
		}
		newType[apply(as.matrix(constraintsVectMatrix[,1:noSParam]),1,sum)==0]<-"asocial"
		newType[typeVect=="asocial"]<-"asocial"

		#Classify model types according to combination of network constraints used
		netCombo<-rep(NA,length(typeVect))
		for(i in 1:length(typeVect)){
		  netCombo[i]<- paste(constraintsVectMatrix[i,1:noSParam],collapse=":")
		}
		#This modifies the netCombo vector if multiple tables are being combined with different network properties
		#most likey weighted versus non-weighted
		netCombo<-paste(netComboModifierVect,netCombo,sep="")

		if(is.null(models)) models<-list(NA)

		if(aicUse=="aic"){
		  printTable<-data.frame(model=1:noModels,type=newType,netCombo=netCombo,baseline=NA,constraintsVectMatrix, offsetVectMatrix,convergence,loglik,MLEs,MLEilv,MLEint,
		                        SEs,SEilv,SEint,aic,aicc,deltaAIC,RelSupport,AkaikeWeight)
		  printTable <-printTable[order(aic),]
		}else{
		  printTable<-data.frame(model=1:noModels,type=newType,netCombo=netCombo,baseline=NA,constraintsVectMatrix, offsetVectMatrix,convergence,loglik,MLEs,MLEilv,MLEint,
		                         SEs,SEilv,SEint,aic,aicc, deltaAICc=deltaAIC,RelSupport,AkaikeWeight)
		  printTable <-printTable[order(aicc),]
		}


		if(is.list(nbdadata)){
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata,models=models,convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)
		}else{
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), models=models,convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)

		}

    }
)


#'Fit a set of OADA models for multi-model inference
#'
#'\code{oadaAICtable} takes diffusion data in the form of an nbdaData object (\code{\link{nbdaData}}) or a list of nbdaData
#'objects (for multiple diffusions). It then fits a set of models using \code{\link{oadaFit}} and return them in an object
#'of class \code{oadaAICtable}. Arguments not described below are used internally by \code{\link{combineOadaAICtables}} when
#'making calls to \code{oadaAICtable} and can be ignored by the user.
#'
#'Each row of \code{constraintsVectMatrix} and \code{offsetVectMatrix} determines a model to
#'be fitted. For each row, constrained \code{nbdaData} objects are created using (\code{\link{constrainedNBDAdata}}) with
#'\code{constraintsVect=constraintsVectMatrix[i,]} and \code{offsetVect=offsetVectMatrix[i,]}. A model is then fitted to
#'the \code{nbdaData} object(s) using \code{\link{oadaFit}}.
#'
#'@seealso \code{\link{typeSupport}}, \code{\link{networksSupport}}, \code{\link{typeByNetworksSupport}}, \code{\link{modelAverageEstimates}},
#' \code{\link{variableSupport}}, \code{\link{unconditionalStdErr}}, \code{\link{combineOadaAICtables}}. For TADA models
#'use \code{\link{tadaAICtable}}.
#'
#'@param nbdadata an object of class (\code{\link{nbdaData}}) to fit models to a single diffusion or a list of
#'nbdaData objects to fit a model to multiple diffusions.
#'@param constraintsVectMatrix a numerical matrix specifying the constraints, with each row specifiying a model to be fitted.
#'The number of columns is equal to the number of parameters in the (\code{\link{nbdaData}}) object(s) input to argument
#'\code{nbdadata}. Each row then specifies the \code{constraintsVect} to be passed to the \code{\link{constrainedNBDAdata}}
#' when creating the data object(s) to which that model is fitted.
#'@param typeVect optional character vector specifying if each model is "asocial" or "social". However, it is not usually
#'necessary to specify, since models with all s parameters constrained =0 are automatically classified as "asocial", and
#'others are assumed to be "social".
#'@param offsetVectMatrix an optional numerical matrix specifying the offsets, with each row specifiying the offsets for each
#'model to be fitted. The number of columns is equal to the number of parameters in the (\code{\link{nbdaData}}) object(s)
#'input to argument \code{nbdadata}. Each row then specifies the \code{offsetVect} to be passed to the
#'\code{\link{constrainedNBDAdata}} when creating the data object(s) to which that model is fitted.
#'@param cores numerical giving the number of computer cores to be used in parallel to fit the models in the set, thus
#'speeding up the process. By default set to 1. For a standard desktop computer at the time of writing 4-6 is advised.
#'@param modelsPerCorePerSet  optional numerical. If specified the models can be fit in sets, and a progress file written
#'after each set is completed. This means progress is not completely lost in the case of a crash/ powercut etc. For example,
#'if we have 400 models, we can specify \code{cores=4} and \code{modelsPerCorePerSet=10}. This means 10 models are fitted on
#'each core, then progress is saved with the first 40 models, then the next 40 and so on.
#'@param writeProgressFile logical. If set to T, a file is written to the working directory when each set of models have
#'been completed with the oadaAICtable for the models fitted so far. In the event of a crash, the remining models can be
#'fitted as a separate set, then combined using \code{\link{combineOadaAICtables}}.
#'@param saveTableList logical. If set to T, a file is written to the working directory containing a list of oadaAICtable
#'objects- one for each set of models per core (see \code{modelsPerCorePerSet}). This is useful if there is an error in one
#'or more of the models preventing the oadaAICtable from being calculated properly. In this case the set(s) containing the
#'faulty model can be corrected, refitted using \code{oadaAICtable}, and replaced in the list. The oadaAICtable can then
#'be reformed using \code{\link{combineOadaAICtables}}.
#'@param statusBar optional logical. Status bar only works when \code{cores=1}.
#'@param startValue optional numeric vector giving start values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the full model.
#'@param lowerList optional numeric matrix giving lower values for the maximum likelihood optimization for each model.
#'Columns to match the number of parameters fitted in the full model, rows matched to the number of models. Can be used if
#'some models have convergence problems or trigger errors.
#'@param upperList optional numeric matrix giving upper values for the maximum likelihood optimization for each model.
#'Columns to match the number of parameters fitted in the full model, rows matched to the number of models. Can be used if
#'some models have convergence problems or trigger errors.
#'@param method optional character string passed to \code{\link{oadaFit}}.
#'@param gradient optional logical passed to \code{\link{oadaFit}}.
#'@param iterations optional numerical passed to \code{\link{oadaFit}}.
#'@param aicUse string specifying whether to use "aicc" or "aic".
#'@param combineTables logical used internally by \code{\link{combineOadaAICtables}} when making calls to \code{oadaAICtable}.
#'@param saveModels logical determining whether the fitted model objects should be saved as a list as a component in the
#'\code{oadaAICtable} object. This can be useful for accessing individual models without re-fitting them (e.g. to get
#'confidence intervals using \code{\link{profLikCI}}). Alternatively individual problem models can be refitted, replaced in
#'the list, and the \code{oadaAICtable} re-generated using \code{oadaAICtable} using the \code{models} argument (see below).
#'@param stripData logical. By default, when \code{saveModels=T} each model is saved with its constrained version of the
#'data. Whilst this can be useful, if there are lots of models and/or a lot of data, this can take up too much memory.
#'setting \code{saveModels=T} and \code{stripData=T} enables the models to be saved without their data.
#'@param models an optional list of model objects of class \code{\link{oadaFit}}. If provided, the models are not fitted
#'but taken from the list when generating the oadaAICtable. This is useful if individual models have been corrected and the
#'user needs to re-form the oadaAICtable.
#'
#'@return An object of class \code{oadaAICtable}.
#'@section print(oadaAICtable)  components:
#'A data.frame giving a summary of models ordered by AIC can be obtained using \code{print(<name of oadaAICtable>)}. This has
#'the following columns, listed in order: \describe{
#'   \item{model}{Model number, i.e. the row of \code{constraintsVectMatrix} used to generate the model.}
#'   \item{type}{Type of model. noILVs, additive (ILV effects on asocial learning only), multiplicative (all ILVs have same
#'   effect on asocial and social learning), unconstrained (differing effects on asocial and social learning for at least one
#'   ILV), socialEffectsOnly (ILV effects only on social learning),asocial}
#'   \item{netcombo}{A representaion of the network effects present in the model, i.e. the constraints on the s
#'   parameters. See \code{\link{constrainedNBDAdata}}.}
#'   \item{baseline}{Used for TADA models only, set to NA here.}
#'   \item{CONS.}{The constraint on each parameter, as taken from \code{constraintsVectMatrix}}
#'   \item{OFF}{The offset on each parameter, as taken from \code{offsetVectMatrix}}
#'   \item{convergence}{Was convergence reported by the optimization algorithm?}
#'   \item{loglik}{-log-likelihood for the model.}
#'   \item{s...}{Maximum likelihood estimates for s parameters.}
#'   \item{ASOCIAL...}{Maximum likelihood estimates for effects of ILVs on asocial learning.}
#'   \item{SOCIAL...}{Maximum likelihood estimates for effects of ILVs on social learning.}
#'   \item{ASOCIAL:SOCIAL...}{Maximum likelihood estimates for multiplicative effects of ILVs see \code{\link{oadaFit}}.}
#'   \item{SE...}{Standard errors for each parameter, set to 0 when a parameter was constrained. See \code{\link{oadaFit}}.}
#'   \item{aic}{AIC for the model.}
#'   \item{aicc}{AICc for the model.}
#'   \item{deltaAICc}{Difference in AICc or AIC from the best model.}
#'   \item{RelSupport}{Relative support for the model compared to the best model, calculated as exp(-0.5*deltaAICc).}
#'   \item{AkaikeWeight}{Akaike weight for the model. Can be interpretted as the probability that model has the highest
#'   predictive power (K-L information) out of the set of models considered.}
#'   }
#'@section oadaAICtable components:\describe{
#'   \item{@@nbdadata}{The unconstrained data the model is fitted to, as a list of nbdaData objects.}
#'   \item{@@models}{A list of fitted model objects of class \code{\link{oadaFit}}, if saved.}
#'   \item{@@convergence}{Was convergence reported by the optimization algorithm? (Ordered by \code{constraintsVectMatrix}).}
#'   \item{@@logLik}{-log-likeihood for models (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@aicc}{AICc for models (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@aic}{AIC for models (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@constraintsVectMatrix}{\code{constraintsVectMatrix} input to the function.}
#'   \item{@@offsetVectMatrix}{\code{offsetVectMatrix} input to the function.}
#'   \item{@@MLEs}{Maximum likelihood estimates for s parameters (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@SEs}{Standard errors for s parameters (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@MLEilv}{Maximum likelihood estimates for effect of ILs on asocial learning (ordered by
#'   \code{constraintsVectMatrix}).}
#'   \item{@@SEilv}{Standard errors for effect of ILs on asocial learning  (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@MLEint}{Maximum likelihood estimates for effect of ILs on social learning (ordered by
#'   \code{constraintsVectMatrix}).}
#'   \item{@@SEint}{Standard errors for effect of ILs on social learning  (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@typeVect}{Type of models (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@deltaAICc}{Difference in AICc or AIC from the best model. (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@RelSupport}{Relative support for the model compared to the best model, calculated as exp(-0.5*deltaAICc).
#'   (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@AkaikeWeight}{Akaike weight for the model. (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@printTable}{data.frame to be output by the print method for \code{oadaAICtable} (see above).}
#'}


#Function for implementing the initialization
oadaAICtable <-function(nbdadata,  constraintsVectMatrix,typeVect=NULL, offsetVectMatrix = NULL, cores=1, modelsPerCorePerSet=NULL,writeProgressFile=F,saveTableList=F,statusBar=NULL,startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,upperList=NULL,combineTables=F,
                        MLEs=NULL,SEs=NULL,MLEilv=NULL,SEilv=NULL,MLEint=NULL,SEint=NULL,
                        convergence=NULL,loglik=NULL,aic=NULL,aicc=NULL,netComboModifierVect="",saveModels=F,stripData=F,models=NULL){

  #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
  if(is.character(nbdadata)){
    newNbdaData<-list()
    for(i in 1:length(nbdadata)){
      newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
    }
    nbdadata<-newNbdaData
  }

  if(cores>1){
    if(is.null(statusBar))statusBar<-F
    oadaAICtable_multiCore(nbdadata=nbdadata,constraintsVectMatrix=constraintsVectMatrix,cores=cores,typeVect=typeVect, offsetVectMatrix = offsetVectMatrix,
                           modelsPerCorePerSet=modelsPerCorePerSet,writeProgressFile=writeProgressFile,statusBar=statusBar, startValue=startValue,method=method, gradient=gradient,iterations=iterations,
                           aicUse=aicUse,lowerList=lowerList,upperList=upperList,saveTableList=saveTableList,saveModels=saveModels,stripData=stripData)

  }else{
    if(is.null(statusBar))statusBar<-T
    return(new("oadaAICtable",nbdadata= nbdadata, typeVect= typeVect, constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix = offsetVectMatrix, startValue= startValue,method= method, gradient= gradient,
	           iterations= iterations,aicUse= aicUse,lowerList=lowerList,upperList=upperList,writeProgressFile=writeProgressFile,combineTables=combineTables,statusBar=statusBar,
	           MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint,
	           convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect,saveModels=saveModels,stripData=stripData,models=models))
  }
}

#Method for printing oadaAICtable
print.oadaAICtable<-function (oadaAICtable)
    {
		oadaAICtable@printTable
	}


#'Get the support for each type of model from an oadaAICtable or tadaAICtable
#'
#'Calculates the support for each type of model fitted in a set of models fitted using \code{\link{oadaAICtable}} or
#'\code{\link{tadaAICtable}}.
#'
#'Models are classified as: \describe{
#'   \item{noILVs;}{}
#'   \item{additive: ILV effects on asocial learning only;}{}
#'   \item{multiplicative: all ILVs have same effect on asocial and social learning;}{}
#'   \item{unconstrained: differing effects on asocial and social learning for at least one ILV;}{}
#'   \item{socialEffectsOnly: ILV effects only on social learning;}{}
#'   \item{asocial: no social transmission (all s=0).}{}
#'   }
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{networksSupport}},
#'\code{\link{typeByNetworksSupport}}, \code{\link{modelAverageEstimates}}, \code{\link{variableSupport}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'
#'@return dataframe giving the support (total Akaike Weight) and number of models of each type.

typeSupport<-function(nbdaAICtable){
	#Calculate support for each type of model in the table
	support<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$type,sum)
	numbers<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$type,length)
	return(data.frame(support=support,numberOfModels=numbers))
}


#'Get the support for each network combination from an oadaAICtable or tadaAICtable
#'
#'Calculates the support for network combination fitted in a set of models fitted using \code{\link{oadaAICtable}} or
#'\code{\link{tadaAICtable}}.
#'
#'Network combinations are represented by the set of constraints on the s parameters used to generate the model, using the
#'\code{constraintsVectMatrix} argument. e.g. 1:1:1 means all three networks in the data were constrained to have the same
#'effect, 1:2:0 means the first two networks have different effects, whereas the third is constrained to have no effect, and
#'0:0:0 represents asocial models.
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{networksSupport}}, \code{\link{modelAverageEstimates}}, \code{\link{variableSupport}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'
#'@return dataframe giving the support (total Akaike Weight) and number of models fitted for each combination.


networksSupport<-function(nbdaAICtable){
  #Calculate support for each combination of network constraints in the table
  support<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$netCombo,sum)
  numbers<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$netCombo,length)
  return(data.frame(support=support,numberOfModels=numbers))
}

#'Get the support for each network x type combination from an oadaAICtable or tadaAICtable
#'
#'Calculates the support for network x type of model combination fitted in a set of models fitted using
#'\code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'
#'Network combinations are represented by the set of constraints on the s parameters used to generate the model, using the
#'\code{constraintsVectMatrix} argument. e.g. 1:1:1 means all three networks in the data were constrained to have the same
#'effect, 1:2:0 means the first two networks have different effects, whereas the third is constrained to have no effect, and
#'0:0:0 represents asocial models. For type, models are classified as: \describe{
#'   \item{noILVs;}{}
#'   \item{additive: ILV effects on asocial learning only;}{}
#'   \item{multiplicative: all ILVs have same effect on asocial and social learning;}{}
#'   \item{unconstrained: differing effects on asocial and social learning for at least one ILV;}{}
#'   \item{socialEffectsOnly: ILV effects only on social learning;}{}
#'   \item{asocial: no social transmission (all s=0).}{}
#'   }
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{typeByNetworksSupport}}, \code{\link{modelAverageEstimates}}, \code{\link{variableSupport}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'
#'@return array giving the support (total Akaike Weight) and number of models fitted for each network combination.


typeByNetworksSupport<-function(nbdaAICtable){
  if(class(nbdaAICtable)=="oadaAICtable"){
    typesList<-levels(nbdaAICtable@printTable$type)
    netComboList<-levels(nbdaAICtable@printTable$netCombo)
    output<-array(NA,dim=c(length(netComboList),nrow=length(typesList),2))
    for(i in 1:length(typesList)){
      output[,i,1]<-tapply(nbdaAICtable@printTable$AkaikeWeight[nbdaAICtable@printTable$type==typesList[i]], nbdaAICtable@printTable$netCombo[nbdaAICtable@printTable$type==typesList[i]],sum)
      output[,i,2]<-tapply(nbdaAICtable@printTable$AkaikeWeight[nbdaAICtable@printTable$type==typesList[i]], nbdaAICtable@printTable$netCombo[nbdaAICtable@printTable$type==typesList[i]],length)
    }
    output[is.na(output[,,2])]<-0
    dimnames(output)<-list(netComboList,typesList,c("Support","NumberOfModels"))
    return(output)
  }
  if(class(nbdaAICtable)=="tadaAICtable"){
    typesList<-levels(nbdaAICtable@printTable$type)
    netComboList<-levels(nbdaAICtable@printTable$netCombo)
    baselineList<-levels(nbdaAICtable@printTable$baseline)
    output<-array(NA,dim=c(length(netComboList),length(typesList),length(baselineList),2))
    for(j in 1:length(baselineList)){
      for(i in 1:length(typesList)){
        output[,i,j,1]<-tapply(nbdaAICtable@printTable$AkaikeWeight[nbdaAICtable@printTable$type==typesList[i]&nbdaAICtable@printTable$baseline==baselineList[j]], nbdaAICtable@printTable$netCombo[nbdaAICtable@printTable$type==typesList[i]&nbdaAICtable@printTable$baseline==baselineList[j]],sum)
        output[,i,j,2]<-tapply(nbdaAICtable@printTable$AkaikeWeight[nbdaAICtable@printTable$type==typesList[i]&nbdaAICtable@printTable$baseline==baselineList[j]], nbdaAICtable@printTable$netCombo[nbdaAICtable@printTable$type==typesList[i]&nbdaAICtable@printTable$baseline==baselineList[j]],length)
      }
    }
    output[is.na(output[,,,2])]<-0
    dimnames(output)<-list(netComboList,typesList,paste("baseline=",baselineList),c("Support","NumberOfModels"))
    return(output)
  }
}

#'Get the support for each variable from an oadaAICtable or tadaAICtable
#'
#'Calculates the support for each variable fitted in a set of models fitted using \code{\link{oadaAICtable}} or
#'code{\link{tadaAICtable}}.
#'
#'The support for each variable is the total Akaike weight for the models in which that variable is present, i.e. not
#'constrained to =0.
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{typeByNetworksSupport}}, \code{\link{networksSupport}}, \code{\link{modelAverageEstimates}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'@param typeFilter  an optional string allowing the user to get the support for variables within a specified type of models.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models. See \code{\link{typeSupport}}
#'for an explanation of the different model types.
#'@param baselineFilter  an optional string allowing the user to get the support for variables within the subset of models
#'with a specific baseline function, e.g. \code{typeFilter="gamma"} gets the variable support in the subset of models with a
#'gamma baseline rate function.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models.
#'@param includeAsocial logical indicating whether asocial models should be included.
#'
#'@return matrix giving the support (total Akaike Weight) for each variable.

variableSupport<-function(nbdaAICtable,typeFilter=NULL,baselineFilter=NULL,includeAsocial=TRUE){
  #Extract the printTable and correct type to include asocial labels
  printTable<-nbdaAICtable@printTable

  #Extract number of s parameters
  noSParam<-dim(nbdaAICtable@MLEs)[2]
  #Extract number of ILVs
  noILVs <-dim(nbdaAICtable@MLEilv)[2]


  #Filter as requested by the user
  if(!is.null(typeFilter)) {
    if(includeAsocial){
      printTable<-printTable[printTable$type==typeFilter|printTable$type=="asocial",]
    }else{
      printTable<-printTable[printTable$type==typeFilter,]
    }
  }
  if(!is.null(baselineFilter)&class(nbdaAICtable)=="tadaAICtable") {
    printTable<-printTable[printTable$baseline==baselineFilter,]

  }

  #Correct Akaike Weights for the new subset of models
  printTable$AkaikeWeight<-printTable$AkaikeWeight/sum(printTable$AkaikeWeight)

  #Set up a vector to record the support for each variable
  support<-rep(NA,dim(nbdaAICtable@constraintsVectMatrix)[2])
  for(i in 1:dim(nbdaAICtable@constraintsVectMatrix)[2]){
    support[i]<-sum(printTable$AkaikeWeight[printTable[,(i+4)]!=0])
  }
  #Convert support into a matrix so I can add dimnames
  support<-rbind(support)
  #Give the support vector some dimension names
  tempNames<-unlist(dimnames(nbdaAICtable@constraintsVectMatrix)[2])
  dimnames(support)[2]=list(gsub("CONS:","",tempNames))

  return(support)

}

#'Get the model averaged estimate for each variable from an oadaAICtable or tadaAICtable
#'
#'Calculates the model averaged estimate for each variable fitted in a set of models fitted using \code{\link{oadaAICtable}}
#'or code{\link{tadaAICtable}}.
#'
#'By default the model averaged estimate is calculated as the weighted mean of maximum likelihood estimates in each model,
#'weighted by Akaike weights. This can be changed to a weighted median using the \code{averageType} argument. This could be
#'useful in cases where a few models (potentially of low weight) have an essentially infinite estimate for a model parameter,
#'e.g. this can occur for s parameters in an OADA when the diffusion follows the network closely. In such cases a Akaike
#'weighted median may better reflect the overall findings of the multi-model analysis. Models where a variable does not appear
#'are treated as models in which the variable has a MLE=0.
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{typeByNetworksSupport}}, \code{\link{networksSupport}}, \code{\link{modelAverageEstimates}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'@param typeFilter  an optional string allowing the user to get the support for variables within a specified type of models.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models. See \code{\link{typeSupport}}
#'for an explanation of the different model types.
#'@param netFilter  an optional string allowing the user to get the support for variables from models with a specific
#'combination of network effects, e.g. \code{netFilter = "0:1:0"} gets the variable support in the subset of models containing
#'only network 2. See \code{\link{networksSupport}} for an explanation of network combinations.
#'@param baselineFilter  an optional string allowing the user to get the support for variables within the subset of models
#'with a specific baseline function, e.g. \code{typeFilter="gamma"} gets the variable support in the subset of models with a
#'gamma baseline rate function.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models.
#'@param includeAsocial logical indicating whether asocial models should be included.
#'@param averageType string indicating whether a "mean" or "median" should be calculated.
#'
#'@return numeric vector giving the model averaged estimate for each variable.

modelAverageEstimates<-function(nbdaAICtable,typeFilter=NULL,netFilter=NULL,baselineFilter=NULL,includeAsocial=TRUE,averageType="mean"){
  #Extract the printTable and correct type to include asocial labels
  printTable<-nbdaAICtable@printTable

  #Extract number of s parameters
  noSParam<-dim(nbdaAICtable@MLEs)[2]
  #Extract number of ILVs
  noILVs <-dim(nbdaAICtable@MLEilv)[2]

  AkaikeWeight <-nbdaAICtable@AkaikeWeight[order(-nbdaAICtable@AkaikeWeight)]
  MLEs<-as.matrix(nbdaAICtable@MLEs[order(-nbdaAICtable@AkaikeWeight),])
  MLEilv<-as.matrix(nbdaAICtable@MLEilv[order(-nbdaAICtable@AkaikeWeight),])
  MLEint<-as.matrix(nbdaAICtable@MLEint[order(-nbdaAICtable@AkaikeWeight),])


  #Filter as requested by the user
  if(!is.null(typeFilter)) {
    if(includeAsocial){
      MLEs <-cbind(MLEs[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      MLEilv <-cbind(MLEilv[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      MLEint <-cbind(MLEint[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      AkaikeWeight<-AkaikeWeight[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial"]
      netCombo<-printTable$netCombo[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial"]
    }else{
      MLEs <-cbind(MLEs[printTable$type==typeFilter|printTable$type=="noILVs",])
      MLEilv <-cbind(MLEilv[printTable$type==typeFilter|printTable$type=="noILVs",])
      MLEint <-cbind(MLEint[printTable$type==typeFilter|printTable$type=="noILVs",])
      AkaikeWeight<-AkaikeWeight[printTable$type==typeFilter|printTable$type=="noILVs"]
      netCombo<-printTable$netCombo[printTable$type==typeFilter|printTable$type=="noILVs"]
    }
  }else{netCombo<-printTable$netCombo}
  #Filter by network as requested by the user
  if(!is.null(netFilter)) {
    MLEs <-MLEs[netCombo==netFilter,]
    MLEilv <-MLEilv[netCombo==netFilter,]
    MLEint <-MLEint[netCombo==netFilter,]
    AkaikeWeight<-AkaikeWeight[netCombo==netFilter]
  }
  if(!is.null(baselineFilter)&class(nbdaAICtable)=="tadaAICtable") {
    MLEs <-cbind(MLEs[printTable$baseline==baselineFilter,])
    MLEilv <-cbind(MLEilv[printTable$baseline==baselineFilter,])
    MLEint <-cbind(MLEint[printTable$baseline==baselineFilter,])
    AkaikeWeight<-AkaikeWeight[printTable$baseline==baselineFilter]
  }

  #Correct Akaike Weights for the new subset of models
  AkaikeWeight<-AkaikeWeight/sum(AkaikeWeight)

  #Do means first and then replace with model weighted medians if requested
  MAvs<-apply(as.matrix(MLEs*AkaikeWeight),2,sum)
  MAvilv<-apply(as.matrix(MLEilv*AkaikeWeight),2,sum)
  MAvint<-apply(as.matrix(MLEint*AkaikeWeight),2,sum)

  if(averageType=="mean"){

  }else{
    if(averageType=="median"){
      for(i in 1:dim(MLEs)[2]){
        tempMLE<-MLEs[,i]
        tempMLEordered<-tempMLE[order(tempMLE)]
        tempAW<-AkaikeWeight[order(tempMLE)]
        cumulAW<-cumsum(tempAW)
        MAvs[i]<-tempMLEordered[min(which(cumulAW>0.5))]
      }
      for(i in 1:dim(MLEilv)[2]){
        tempMLE<-MLEilv[,i]
        tempMLEordered<-tempMLE[order(tempMLE)]
        tempAW<-AkaikeWeight[order(tempMLE)]
        cumulAW<-cumsum(tempAW)
        MAvilv[i]<-tempMLEordered[min(which(cumulAW>0.5))]
      }
      for(i in 1:dim(MLEint)[2]){
        tempMLE<-MLEint[,i]
        tempMLEordered<-tempMLE[order(tempMLE)]
        tempAW<-AkaikeWeight[order(tempMLE)]
        cumulAW<-cumsum(tempAW)
        MAvint[i]<-tempMLEordered[min(which(cumulAW>0.5))]
      }

    }else{
      print("Invalid averageType, please select 'mean' or 'median'");
      return(NULL)
    }
  }


  return(c(MAvs, MAvilv, MAvint))
}

#'Get the unconditional standard error for each variable from an oadaAICtable or tadaAICtable
#'
#'Calculates the unconditional standard error for each variable fitted in a set of models fitted using \code{\link{oadaAICtable}}
#'or code{\link{tadaAICtable}}.
#'
#'A standard error (SE) from an individual model is conditional on that model. The unconditional standard error (USE) is not
#'really  "unconditional" but is rather conditional on the whole set of models fitted rather than an individual model. It takes
#'into account the SEs within each model and the variation in maximum likelihood estimates among models. SEs cannot be always be
#'calulated in an NBDA, meaning we are left with a situation in which USEs cannot be calculated if 1 or more models have SE=NaN-
#'even if there are only a few models of low weight for which this is the case. Here we offer a pragmatic
#'solution to such cases- by replacing the SEs in such models with a model-averaged SE across all models in which the SE can be
#'calculated, we are able to estimate a USE that is likely to be approximately correct. However, we recommend that this solution
#'only be used if SEs are absent for only a few models of low weight. In other cases we recommend that the user use confidence intervals
#'from \code{\link{profLikCI}} conditional on the best model as a measure of uncertainty in parameter estimates.
#'For s parameters \code{\link{multiModelLowerLimits}} and \code{\link{multiModelPropST}} then provide a means to assess the
#'robustness of the findings to model selection uncertainty.
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{typeByNetworksSupport}}, \code{\link{networksSupport}}, \code{\link{modelAverageEstimates}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'@param typeFilter  an optional string allowing the user to get the support for variables within a specified type of models.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models. See \code{\link{typeSupport}}
#'for an explanation of the different model types.
#'@param netFilter  an optional string allowing the user to get the support for variables from models with a specific
#'combination of network effects, e.g. \code{netFilter = "0:1:0"} gets the variable support in the subset of models containing
#'only network 2. See \code{\link{networksSupport}} for an explanation of network combinations.
#'@param baselineFilter  an optional string allowing the user to get the support for variables within the subset of models
#'with a specific baseline function, e.g. \code{typeFilter="gamma"} gets the variable support in the subset of models with a
#'gamma baseline rate function.
#'e.g. \code{typeFilter="additive"} gets the variable support in the subset of additive models.
#'@param includeAsocial logical indicating whether asocial models should be included.
#'@param includeNoILVs logical indicating whether models with no ILVs should be included.
#'@param nanReplace logical indicating whether standard errors recorded as NaNs should be replaced with the model averaged mean
#'SE across all models for which the SE for that variable could be calculated.
#'
#'@return numeric vector giving the model averaged estimate for each variable.


unconditionalStdErr<-function(nbdaAICtable,typeFilter=NULL,netFilter=NULL,baselineFilter=NULL,includeAsocial=TRUE,includeNoILVs=TRUE,nanReplace=FALSE){
  #Extract the printTable and correct type to include asocial labels
  #Extract the printTable and correct type to include asocial labels
  printTable<-nbdaAICtable@printTable

  #Extract number of s parameters
  noSParam<-dim(nbdaAICtable@MLEs)[2]
  #Extract number of ILVs
  noILVs <-dim(nbdaAICtable@MLEilv)[2]

  AkaikeWeight <-nbdaAICtable@AkaikeWeight[order(-nbdaAICtable@AkaikeWeight)]
  MLEs<-as.matrix(nbdaAICtable@MLEs[order(-nbdaAICtable@AkaikeWeight),])
  MLEilv<-as.matrix(nbdaAICtable@MLEilv[order(-nbdaAICtable@AkaikeWeight),])
  MLEint<-as.matrix(nbdaAICtable@MLEint[order(-nbdaAICtable@AkaikeWeight),])
  SEs<-as.matrix(nbdaAICtable@SEs[order(-nbdaAICtable@AkaikeWeight),])
  SEilv<-as.matrix(nbdaAICtable@SEilv[order(-nbdaAICtable@AkaikeWeight),])
  SEint<-as.matrix(nbdaAICtable@SEint[order(-nbdaAICtable@AkaikeWeight),])


  #Filter as requested by the user
  if(!is.null(typeFilter)) {
    if(includeAsocial){
      MLEs <-cbind(MLEs[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      MLEilv <-cbind(MLEilv[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      MLEint <-cbind(MLEint[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      AkaikeWeight<-AkaikeWeight[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial"]
      netCombo<-printTable$netCombo[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial"]
    }else{
      MLEs <-cbind(MLEs[printTable$type==typeFilter|printTable$type=="noILVs",])
      MLEilv <-cbind(MLEilv[printTable$type==typeFilter|printTable$type=="noILVs",])
      MLEint <-cbind(MLEint[printTable$type==typeFilter|printTable$type=="noILVs",])
      AkaikeWeight<-AkaikeWeight[printTable$type==typeFilter|printTable$type=="noILVs"]
      netCombo<-printTable$netCombo[printTable$type==typeFilter|printTable$type=="noILVs"]
    }
  }else{netCombo<-printTable$netCombo}
  #Filter by network as requested by the user
  if(!is.null(netFilter)) {
    MLEs <-MLEs[netCombo==netFilter,]
    MLEilv <-MLEilv[netCombo==netFilter,]
    MLEint <-MLEint[netCombo==netFilter,]
    AkaikeWeight<-AkaikeWeight[netCombo==netFilter]
  }
  if(!is.null(baselineFilter)&class(nbdaAICtable)=="tadaAICtable") {
    MLEs <-MLEs[printTable$baseline==baselineFilter,]
    MLEilv <-MLEilv[printTable$baseline==baselineFilter,]
    MLEint <-MLEint[printTable$baseline==baselineFilter,]
    AkaikeWeight<-AkaikeWeight[printTable$baseline==baselineFilter]
  }

  #Filter as requested by the user
  if(!is.null(typeFilter)) {
    if(includeAsocial){
      SEs <-cbind(SEs[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      SEilv <-cbind(SEilv[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
      SEint <-cbind(SEint[printTable$type==typeFilter|printTable$type=="noILVs"|printTable$type=="asocial",])
    }else{
      SEs <-cbind(SEs[printTable$type==typeFilter|printTable$type=="noILVs",])
      SEilv <-cbind(SEilv[printTable$type==typeFilter|printTable$type=="noILVs",])
      SEint <-cbind(SEint[printTable$type==typeFilter|printTable$type=="noILVs",])
    }
  }else{netCombo<-printTable$netCombo}
  #Filter by network as requested by the user
  if(!is.null(netFilter)) {
    SEs <-SEs[netCombo==netFilter,]
    SEilv <-SEilv[netCombo==netFilter,]
    SEint <-SEint[netCombo==netFilter,]
  }
  if(!is.null(baselineFilter)&class(nbdaAICtable)=="tadaAICtable") {
    SEs <-SEs[printTable$baseline==baselineFilter,]
    SEilv <-SEilv[printTable$baseline==baselineFilter,]
    SEint <-SEint[printTable$baseline==baselineFilter,]
  }

  #If nanReplace option is used, then nan for individual model SEs are replaced with a weighted average across all model containing that parameter
  if(nanReplace){
    for(i in 1:dim(SEs)[2]){
      SEs[is.nan(SEs[,i]),i]<-sum(SEs[!is.nan(SEs[,i])&(SEs[,i]>0),i]*AkaikeWeight[!is.nan(SEs[,i])&(SEs[,i]>0)])/sum(AkaikeWeight[!is.nan(SEs[,i])&(SEs[,i]>0)])
    }
    for(i in 1:dim(SEilv)[2]){
      SEilv[is.nan(SEilv[,i]),i]<-sum(SEilv[!is.nan(SEilv[,i])&(SEilv[,i]>0),i]*AkaikeWeight[!is.nan(SEilv[,i])&(SEilv[,i]>0)])/sum(AkaikeWeight[!is.nan(SEilv[,i])&(SEilv[,i]>0)])
    }
    for(i in 1:dim(SEint)[2]){
      SEint[is.nan(SEint[,i]),i]<-sum(SEint[!is.nan(SEint[,i])&(SEint[,i]>0),i]*AkaikeWeight[!is.nan(SEint[,i])&(SEint[,i]>0)])/sum(AkaikeWeight[!is.nan(SEint[,i])&(SEint[,i]>0)])
    }
  }

  #Correct Akaike Weights for the new subset of models
  AkaikeWeight<-AkaikeWeight/sum(AkaikeWeight)

  MAvs<-apply(as.matrix(MLEs*AkaikeWeight),2,sum)
  MAvilv<-apply(as.matrix(MLEilv*AkaikeWeight),2,sum)
  MAvint<-apply(as.matrix(MLEint*AkaikeWeight),2,sum)

  modelContributionToSEs<-AkaikeWeight*(SEs^2 + t((t(MLEs)-MAvs))^2)
  modelContributionToSEilv<-AkaikeWeight*(SEilv^2 + t((t(MLEilv)-MAvilv))^2)
  modelContributionToSEint<-AkaikeWeight*(SEint^2 + t((t(MLEint)-MAvint))^2)

  UCSEs<-apply(as.matrix(modelContributionToSEs),2,sum)
  UCSEilv<-apply(as.matrix(modelContributionToSEilv),2,sum)
  UCSEint<-apply(as.matrix(modelContributionToSEint),2,sum)

  return(c(UCSEs, UCSEilv, UCSEint))
}

#'Combine two or more oadaAICtables into a single oadaAICtable
#'
#'Takes two or more \code{\link{oadaAICtable}} objects, containing different models including the same number of networks and
#'the same ILVs, and combine them into a single table.
#'
#'This function can be used for a variety of practical reasons. If the \code{\link{oadaAICtable}} is interrupted and a partial
#'oadaAICtable recovered, the remaining models can be run as a separate set and combined with the first partial set. Alternatively,
#'the user might want to add more models to a set without re-running the entire set, or run different parts of a model set on
#'different computers to speed up computation time. It is also possible to run model sets with different networks and then combine
#'them using \code{combineOadaAICtables}, so long as the number of networks matches (though it is possible to have dummy networks
#'that are constrained to have s=0 in all models, in order to match network number). In such cases a modifier can be added to the
#'netcombo codes from each oadaAICtable so they can be distinguished by functions such as \code{\link{networksSupport}}. e.g.
#'if we have two oadaAICtable objects each with the netcombos "1:0","0:1","1:2". If we use
#'\code{netComboModifier= c("NetworksA:","NetworksB:")} we get netcombos "NetworksA:1:0","NetworksA:0:1","NetworksA:1:2",
#'"NetworksB:1:0","NetworksB:0:1","NetworksB:1:2" in the resulting oadaAICtable object.
#'
#'@seealso \code{\link{oadaAICtable}}
#'
#'@param oadaAICtableList a lift of \code{\link{oadaAICtable}} objects to be combined into a single table.
#'@param aicUse string specifying whether to use "aicc" or "aic".
#'@param netComboModifier optional character vector with length matching the length of \code{oadaAICtableList} to modify the
#'netcombo strings recorded from each individual \code{\link{oadaAICtable}} object (see below).
#'@return An object of class \code{oadaAICtable}.


combineOadaAICtables<-function(oadaAICtableList,aicUse="aicc",netComboModifier=rep("",length(oadaAICtableList))){

  i<-1
  oadaAICtableTemp<-oadaAICtableList[[i]]
  nbdadata<-oadaAICtableTemp@nbdadata

  writeProgressFile=F;
  startValue=NULL;
  method=NULL;
  gradient=NULL;
  iterations=NULL;
  lowerList=NULL;
  upperList=NULL;

  typeVect<-oadaAICtableTemp@typeVect;
  constraintsVectMatrix<-oadaAICtableTemp@constraintsVectMatrix;
  offsetVectMatrix<-oadaAICtableTemp@offsetVectMatrix;
  netComboModifierVect<-rep(netComboModifier[i],dim(oadaAICtableTemp@constraintsVectMatrix)[1])
  models<-oadaAICtableTemp@models

  #Loop through the oadaAICtableList and create combined versions of typeVect, constraintsVectMatrix, and offsetVectMatrix
  for(i in 2:length(oadaAICtableList)){
    #Read in the next oadaAICtable
    oadaAICtableTemp<-oadaAICtableList[[i]]
    typeVect<-c(typeVect,oadaAICtableTemp@typeVect)
    constraintsVectMatrix<-rbind(constraintsVectMatrix,oadaAICtableTemp@constraintsVectMatrix)
    offsetVectMatrix<-rbind(offsetVectMatrix,oadaAICtableTemp@offsetVectMatrix)
    netComboModifierVect<-c(netComboModifierVect,rep(netComboModifier[i],dim(oadaAICtableTemp@constraintsVectMatrix)[1]))
    models<-c(models,oadaAICtableTemp@models)
  }

  typeVect[typeVect!="asocial"]<-"social"
  typeVect[is.na(typeVect)]<-"social"


  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    nbdadataTemp1<-nbdadata[[1]];
  }else{nbdadataTemp1<-nbdadata}

  noModels<-dim(constraintsVectMatrix)[1]
  dimnames(constraintsVectMatrix)[[1]]<-dimnames(offsetVectMatrix)[[1]]<-1:noModels


  #Calculate the number of different s parameters, ILVs and models to be fitted
  noSParam<-dim(nbdadataTemp1@stMetric)[2]
  noILVasoc<- dim(nbdadataTemp1@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp1@intILVdata)[2] #ILV effects on interation (social learning)
  noILVmulti<- dim(nbdadataTemp1@multiILVdata)[2] #ILV multiplicative model effects
  if(nbdadataTemp1@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp1@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp1@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  #Record asocialVar names
  asocialVarNames<-unique(c(nbdadataTemp1@asoc_ilv,nbdadataTemp1@int_ilv,nbdadataTemp1@multi_ilv))
  asocialVarNames<-asocialVarNames[asocialVarNames!="ILVabsent"]
  if(is.null(asocialVarNames)){noILVs<-0}else{noILVs<-length(asocialVarNames)}

  #Set up matrices to record maximum likelihood estimators and SEs
  MLEs<-matrix(NA,nrow=noModels,ncol=noSParam,dimnames=list(1:noModels, paste("s",1:noSParam,sep="")))
  SEs<-matrix(NA,nrow=noModels,ncol=noSParam,dimnames=list(1:noModels, paste("SEs",1:noSParam,sep="")))
  MLEilv<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("ASOCIAL",asocialVarNames,sep="")))
  SEilv<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SEasocial",asocialVarNames,sep="")))
  MLEint<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SOCIAL",asocialVarNames,sep="")))
  SEint<-matrix(0,nrow=noModels,ncol= noILVs, dimnames=list(1:noModels, paste("SEsocial",asocialVarNames,sep="")))

  #Set up various vectors to record things about each model
  convergence<-loglik<-aic<-aicc<-seApprox<-rep(NA,noModels)

  #Loop through the oadaAICtableList and fill in the various matrices and vectors
  #Track the start point to fill in each time
  startPoint<-1
  for(i in 1:length(oadaAICtableList)){
    #Read in the next oadaAICtable
    oadaAICtableTemp<-oadaAICtableList[[i]]
    endPoint<-dim(oadaAICtableTemp@constraintsVectMatrix)[1]+startPoint-1

    MLEs[startPoint:endPoint,]<-oadaAICtableTemp@MLEs
    SEs[startPoint:endPoint,]<-oadaAICtableTemp@SEs
    if(noILVasoc!=0){
      MLEilv[startPoint:endPoint,]<-oadaAICtableTemp@MLEilv
      SEilv[startPoint:endPoint,]<-oadaAICtableTemp@SEilv
      MLEint[startPoint:endPoint,]<-oadaAICtableTemp@MLEint
      SEint[startPoint:endPoint,]<-oadaAICtableTemp@SEint
    }
    convergence[startPoint:endPoint]<-oadaAICtableTemp@convergence
    loglik[startPoint:endPoint]<-oadaAICtableTemp@loglik
    aic[startPoint:endPoint]<-oadaAICtableTemp@aic
    aicc[startPoint:endPoint]<-oadaAICtableTemp@aicc
    startPoint<-endPoint+1
  }

  tableTemp<-oadaAICtable(nbdadata=nbdadata,  constraintsVectMatrix=constraintsVectMatrix,typeVect=typeVect, offsetVectMatrix = offsetVectMatrix,
               startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,upperList=upperList,writeProgressFile=F,
               combineTables=T, models=models,
               MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint,
               convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect)

  return(tableTemp)

}

#This function runs an oadaAICtable with multiple cores, but is called via oadaAICtable with the cores argument
#so does not need calling by the user
oadaAICtable_multiCore<-function(nbdadata,constraintsVectMatrix,cores,typeVect=NULL, offsetVectMatrix = NULL, modelsPerCorePerSet=NULL,writeProgressFile=F,statusBar=F,startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,upperList=NULL,saveTableList=F,saveModels=F,stripData=F){
  noModels<-dim(constraintsVectMatrix)[1];
  #If cores is more than the number of models, reduce so one model is fit per cors
  if(cores>noModels) cores<-noModels
  #How many models will be run on each core?
  modelsPerCore<-noModels%/%cores
  #How many will be left over to run at the end?
  remainderModels<-noModels%%cores

  #I will set the models to run in sets of modelsPerCorePerSet per core (modelsPerCore by default) then write the oadaAICtable so far to a
  #file in the working directory if writeProgressFile=T (the default is F)
  #This means in a long model run, if there is a powercut or crash the progress will be saved
  #This slows things down so only worth it for big models that take a long time to fit
  if(is.null(modelsPerCorePerSet)){modelsPerCorePerSet<-modelsPerCore}
  if(modelsPerCorePerSet>modelsPerCore)modelsPerCorePerSet<-modelsPerCore
  numberOfInitialSets<-modelsPerCore%/%modelsPerCorePerSet
  remainderModelsPerCore<-modelsPerCore%%modelsPerCorePerSet

  #aicTable funcions do not work with a single model, so adjust to make sure a single model is not fitted
  while(remainderModelsPerCore==1|remainderModels==1){
    cores<-cores-1
    #How many models will be run on each core?
    modelsPerCore<-noModels%/%cores
    #How many will be left over to run at the end?
    remainderModels<-noModels%%cores
    #How many models will be run on each core?
    numberOfInitialSets<-modelsPerCore%/%modelsPerCorePerSet
    remainderModelsPerCore<-modelsPerCore%%modelsPerCorePerSet
    if(cores==1) break
  }

  #Do the remainder set first
  if(remainderModels>0){
    constraintsVectMatrixTemp<- constraintsVectMatrix[1:remainderModels,]
    if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[1:remainderModels]}
    if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[1:remainderModels,]}

  cumulativeAICtable<-oadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,
                                   startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,upperList=upperList,writeProgressFile=F,statusBar = F,saveModels=saveModels,stripData=stripData)
  aicTableList<-list(cumulativeAICtable)
  }else{cumulativeAICtable<-aicTableList<-NULL}

  if(statusBar) pb <- txtProgressBar(min=0, max=numberOfInitialSets, style=3)

  #Loop through the initial sets
  for(set in 1:numberOfInitialSets){

    #Set up for parallel processing as detailed in http://www.parallelr.com/r-with-parallel-computing/
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)

    tablesFromSet <- foreach(i=1:cores) %dopar%
    {
      #I think we need to reload the NBDA package into each thread
      library(NBDA)
      library(coxme)
      #Identify which models need to be fitted in this core x set combination
      modelSet<-((set-1)*modelsPerCorePerSet*cores + (i-1)*modelsPerCorePerSet+remainderModels)+1:modelsPerCorePerSet
      #Cut down constraintsVectMatrix and, if necessary typeVect and offsetVectMatrix
      constraintsVectMatrixTemp<- constraintsVectMatrix[modelSet,]
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-oadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,
                                       startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,upperList=upperList,writeProgressFile=F,saveModels=saveModels,stripData=stripData)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative oadaAICtable so far
    cumulativeAICtable<-combineOadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
    aicTableList<-c(aicTableList,tablesFromSet)
    if(writeProgressFile) save(cumulativeAICtable,file="cumulativeAICtable_parallel.Rdata")
    if(statusBar)setTxtProgressBar(pb, set)

  }
  #Do the remainder for each core
  if(remainderModelsPerCore>0){
    #Set up for parallel processing as detailed in http://www.parallelr.com/r-with-parallel-computing/
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)

    tablesFromSet <- foreach(i=1:cores) %dopar%
    {
      #I think we need to reload the NBDA package into each thread
      library(NBDA)
      library(coxme)
      #Identify which models need to be fitted in this core x set combination
      modelSet<-(numberOfInitialSets*modelsPerCorePerSet*cores + (i-1)*remainderModelsPerCore+remainderModels)+1:remainderModelsPerCore
      #Cut down constraintsVectMatrix and, if necessary typeVect and offsetVectMatrix
      constraintsVectMatrixTemp<- constraintsVectMatrix[modelSet,]
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-oadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,
                                 startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,upperList=upperList,writeProgressFile=F,saveModels=saveModels,stripData=stripData)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative oadaAICtable so far
    cumulativeAICtable<-combineOadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
    aicTableList<-c(aicTableList,tablesFromSet)
  }
  if(saveTableList) save(aicTableList,file="finalAICtable_parallel_list.Rdata")
  if(writeProgressFile) save(cumulativeAICtable,file="finalAICtable_parallel.Rdata")
  if(statusBar)close(pb)

  return(cumulativeAICtable)
}


#'Test the sensitivity of the lower limit of an s parameter to model selection uncertainty
#'
#'\code{multiModelLowerLimits} returns the lower endpoint of the confidence interval for a specified s parameter for all models
#'(or just the top models) in a set of models fitted using \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}. It also
#'returns the estimated proportion of events that occurred by social transmission via the corresponding network, if social
#'tranmssion occurred at this lowest plausible rate, for each model. Only models including the chosen s parameter are included.
#'Since this can take a long time, \code{\link{multiModelLowerLimits_multicore}} is available to run the function in parallel
#'on multiple computer cores.
#'
#'The goal of this function is to test if conclusions about social transmission are robust to model selection
#'uncertainty. Often unconditional standard errors (USEs) \code{\link{unconditionalStdErr}} are used to allow for model selection
#'uncertainty. However, these are often inappropriate for s parameters in an NBDA, due to the high asymmetry in the profile
#'likelihood for s parameters. In other words, we can have a lot of information about the lower plausible limit for social
#'transmission rate, but little information about the upper plausible limit. Standard errors, or USEs only reflect overall levels
#'of information, so can make it appear like there is little evidence for social transmission when in fact there is a strong
#'evidence. The solution is to obtain confidence intervals using the profile likelihood method for s parameters (at least) using
#'\code{\link{unconditionalStdErr}}.
#'
#'However, since these confidence intervals are conditional on a single model (usually the
#'top model by AICc), it makes sense to test whether our conclusions are robust to model selection uncertainty. This function
#'returns the lower endpoint of the 95\% (by default) confidence interval: if this endpoint is well away from zero for all models,
#'it indicates that our conclusion is highly robust to model selection uncertainty. If the endpoint is well away from zero for
#'most models with a sizeable Akaike weight, it indicates that our conclusion is moderately robust to model selection
#'uncertainty, and so on.
#'
#'Since it can be difficult to interpret whether an s parameter is "far from zero" the function also provides the corresponding
#'propST, the proportion of events estimated to have occurred by social transmission via that network, see
#'\code{\link{nbdaPropSolveByST}}.
#'
#'
#'@seealso \code{\link{multiModelLowerLimits_multicore}}, \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}},
#'\code{\link{multiModelPropST}}, \code{\link{plotProfLik}}, \code{\link{nbdaPropSolveByST}}.
#'
#'@param aicTable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'@param which numerical giving the s parameter for which the lower confidence interval endpoints are to be calculated.
#'@param deltaThreshold optional numerical determining the threshold difference in AICc/AIC for a model to be included in the
#'output. e.g. \code{deltaThreshold=10} includes all models within 10 AICc units of the best model.
#'@param modelIndex optional numeric vector specifiying which models to include in the output, subject to \code{deltaThreshold}.
#'@param searchRange optional numeric vector of length two, giving the range within which to search for the lower endpoint. If
#'omitted, the function searches between 0 and the MLE for s in each model.
#'@inheritParams oadaAICtable
#'@inheritParams profLikCI
#'@inheritParams nbdaPropSolveByST
#'
#'@return data.frame giving: \describe{
#'   \item{model: model number;}{}
#'   \item{netCombo: the network combination for the model;}{}
#'   \item{lowerCI: the lower endpoint of the confidence interval for the chosen s parameter;}{}
#'   \item{propST: the proportion of social transmission events estimated to have occurred by social transmission, corresponding
#'   to lowerCI;}{}
#'   \item{deltaAICc: difference in AICc from the best model in the original set;}{}
#'   \item{akaikeWeight: Akaike weight in the original model set, see \code{\link{oadaAICtable}};}{}
#'   \item{adjAkWeight: Akaike weight adjusted to the set of models considered here;}{}
#'   \item{cumulAdjAkWeight: Cumulative adjusted Akaike weight.}{}
#'   }


multiModelLowerLimits<-function(which,aicTable,deltaThreshold=Inf,conf=0.95,modelIndex=NULL,searchRange=NULL,exclude.innovations=T,innovations=NULL,startValue=NULL,lowerList=NULL,
                                upperList=NULL,
                                method="nlminb", gradient=T,iterations=150){

  #Cut down to the models specificied in modelIndex (this is mostly included to allow multicore implementation)
  if(is.null(modelIndex)) modelIndex<-1:dim(aicTable@printTable)[1]



  if(class(aicTable)=="oadaAICtable"){

  #get the model set in which AIC or AICc was < the threshold and in which the parameter was present
  modelsIncluded<-aicTable@printTable$deltaAIC<deltaThreshold&aicTable@printTable[,which+4]
  modelSet<-aicTable@printTable$model[modelsIncluded]
  modelSet<-modelSet[modelSet%in%modelIndex]

  if(length(modelSet)==0) return(NULL)

  lowerLimitRecord<-lowerLimitPropST<-rep(NA,length(modelSet))

  nbdadata<-aicTable@nbdadata

  #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
  if(is.list(nbdadata)){
    nbdadataTemp1<-nbdadata[[1]]
  }else{nbdadataTemp1<-nbdadata}

  #Calculate the number of different s parameters, ILVs and models to be fitted
  noSParam<-dim(nbdadataTemp1@stMetric)[2]
  noILVasoc<- dim(nbdadataTemp1@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadataTemp1@intILVdata)[2] #ILV effects on interation (social learning)
  noILVmulti<- dim(nbdadataTemp1@multiILVdata)[2] #ILV multiplicative model effects
  if(nbdadataTemp1@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadataTemp1@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadataTemp1@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  pb <- txtProgressBar(min=0, max=length(modelSet), style=3)

  for(i in 1:length(modelSet)){
    setTxtProgressBar(pb, i)

    model<-modelSet[i]
    constraintsVect<-aicTable@constraintsVectMatrix[model,]
    offsetVect<-aicTable@offsetVectMatrix[model,]

    if(is.null(startValue)) {
      newStartValue<-NULL
    }else{
      newStartValue<-startValue[constraintsVect!=0]
    }

    if(is.null(lowerList)) {
      lower<-NULL
    }else{
      lower<-lowerList[i,]
      lower<-lower[constraintsVect!=0]
    }
    if(is.null(upperList)) {
      upper<-NULL
    }else{
      upper<-upperList[i,]
      upper<-upper[constraintsVect!=0]
    }
    #Create the necessary constrained data objects
    if(is.list(nbdadata)){
      nbdadataTemp<-list()
      for(dataset in 1:length(nbdadata)){
        nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
      }
    }else{
      nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
    }

    #Fit the model
    modelFit<-NULL
    try(modelFit<-oadaFit(nbdadata= nbdadataTemp,startValue=newStartValue,method=method,gradient=gradient,iterations=iterations))
    if(distanceFromCutoff(value=0,which=constraintsVect[which],model=modelFit,conf=conf,direction=T)>0){
      lowerLimitRecord[i]<-0
    }else{
      #Get the lower 95% CI
      if(is.null(searchRange)){
        lowerLimitRecord[i]<-profLikCI(which=constraintsVect[which],model=modelFit,conf=conf,lowerRange=c(0,modelFit@outputPar[constraintsVect[which]]))[1]
      }else{
        lowerLimitRecord[i]<-profLikCI(which=constraintsVect[which],model=modelFit,conf=conf,lowerRange=searchRange)[1]
      }
      }


    #Fit lower limit model

    #Create the necessary constrained data objects
    newConstraintsVect<-constraintsVect
    newConstraintsVect[which]<-0
    newConstraintsVect[newConstraintsVect>0]<-as.numeric(factor(newConstraintsVect[newConstraintsVect>0]))
    newOffsetVect<-offsetVect
    newOffsetVect[which]<-newOffsetVect[constraintsVect[which]]+lowerLimitRecord[i]

    type="social"
    if(sum(newConstraintsVect[1:noSParam])==0){
      type<-"asocial";
      newConstraintsVect[1]<-1;
      newConstraintsVect[-(1:noSParam)]<-(newConstraintsVect[-(1:noSParam)]+1)*(newConstraintsVect[-(1:noSParam)]>0);
    }

    if(is.list(nbdadata)){
      nbdadataTemp<-list()
      for(dataset in 1:length(nbdadata)){
        nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=newConstraintsVect,offsetVect=newOffsetVect))
      }
    }else{
      nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=newConstraintsVect,offsetVect=newOffsetVect)
    }

    lowerLimModel<-oadaFit(nbdadataTemp,type=type)
    if(lowerLimModel@varNames[1]=="No variables"){
      propST<-oadaPropSolveByST(par=0,nbdadata=lowerLimModel@nbdadata,exclude.innovations=exclude.innovations,innovations=innovations)
    }else{
      propST<-oadaPropSolveByST(model=lowerLimModel,exclude.innovations=exclude.innovations,innovations=innovations)
    }
    lowerLimitPropST[i]<-propST[length(propST)]
  }




  lowerLimitRecordFull<-lowerLimitPropSTFull<-rep(NA,dim(aicTable@printTable)[1])
  lowerLimitRecordFull[modelsIncluded]<-lowerLimitRecord
  lowerLimitPropSTFull[modelsIncluded]<-lowerLimitPropST

  output1<-data.frame(model=aicTable@printTable$model,netCombo=aicTable@printTable$netCombo,lowerCI=lowerLimitRecordFull,propST=lowerLimitPropSTFull,deltaAICc=aicTable@printTable$deltaAICc,
                      akaikeWeight=aicTable@printTable$AkaikeWeight
  )[!is.na(lowerLimitRecordFull),]

  output<-cbind(output1,adjAkWeight=output1$akaikeWeight/sum(output1$akaikeWeight),
                cumulAdjAkWeight=cumsum(output1$akaikeWeight/sum(output1$akaikeWeight)))

  return(output)
  }

  if(class(aicTable)=="tadaAICtable"){
    #get the model set in which AIC or AICc was < the threshold and in which the parameter was present
    modelsIncluded<-aicTable@printTable$deltaAIC<deltaThreshold&aicTable@printTable[,which+4]
    modelSet<-aicTable@printTable$model[modelsIncluded]
    modelSet<-modelSet[modelSet%in%modelIndex]

    if(length(modelSet)==0) return(NULL)


    lowerLimitRecord<-lowerLimitPropST<-rep(NA,length(modelSet))

    nbdadata<-aicTable@nbdadata

    #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
    if(is.list(nbdadata)){
      nbdadataTemp1<-nbdadata[[1]]
    }else{nbdadataTemp1<-nbdadata}

    #Calculate the number of different s parameters, ILVs and models to be fitted
    noSParam<-dim(nbdadataTemp1@stMetric)[2]
    noILVasoc<- dim(nbdadataTemp1@asocILVdata)[2] #ILV effects on asocial learning
    noILVint<- dim(nbdadataTemp1@intILVdata)[2] #ILV effects on interation (social learning)
    noILVmulti<- dim(nbdadataTemp1@multiILVdata)[2] #ILV multiplicative model effects
    if(nbdadataTemp1@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
    if(nbdadataTemp1@int_ilv[1]=="ILVabsent") noILVint<-0
    if(nbdadataTemp1@multi_ilv[1]=="ILVabsent") noILVmulti<-0

    pb <- txtProgressBar(min=0, max=length(modelSet), style=3)

    for(i in 1:length(modelSet)){
      setTxtProgressBar(pb, i)

      model<-modelSet[i]
      constraintsVect<-aicTable@constraintsVectMatrix[model,]
      offsetVect<-aicTable@offsetVectMatrix[model,]
      baseline<-aicTable@baselineVect[i]

      if(is.null(startValue)) {
        newStartValue<-NULL
      }else{
        newStartValue<-startValue[constraintsVect!=0]
      }

      if(is.null(lowerList)) {
        lower<-NULL
      }else{
        lower<-lowerList[i,]
        lower<-lower[constraintsVect!=0]
      }
      #Create the necessary constrained data objects
      if(is.list(nbdadata)){
        nbdadataTemp<-list()
        for(dataset in 1:length(nbdadata)){
          nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=constraintsVect,offsetVect=offsetVect))
        }
      }else{
        nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=constraintsVect,offsetVect=offsetVect)
      }

      #Fit the model

      if(baseline=="constant") noHazFunctPars<-1
      if(baseline=="gamma") noHazFunctPars<-2
      if(baseline=="weibull") noHazFunctPars<-2

      modelFit<-NULL
      try(modelFit<-tadaFit(nbdadata= nbdadataTemp,baseline=baseline,startValue=newStartValue,method=method,gradient=gradient,iterations=iterations))
      if(distanceFromCutoff(value=0,which=constraintsVect[which],model=modelFit,conf=conf,direction=T)>0){
        lowerLimitRecord[i]<-0
      }else{
        #Get the lower 95% CI
        lowerLimitRecord[i]<-profLikCI(which=constraintsVect[which],model=modelFit,conf=conf,lowerRange=c(0,modelFit@outputPar[constraintsVect[which]+noHazFunctPars]))[1]
      }
      #Fit lower limit model

      #Create the necessary constrained data objects
      newConstraintsVect<-constraintsVect
      newConstraintsVect[which]<-0
      newConstraintsVect[newConstraintsVect>0]<-as.numeric(factor(newConstraintsVect[newConstraintsVect>0]))
      newOffsetVect<-offsetVect
      newOffsetVect[which]<-newOffsetVect[constraintsVect[which]]+lowerLimitRecord[i]

      type="social"
      if(sum(newConstraintsVect[1:noSParam])==0){
        type<-"asocial";
        newConstraintsVect[1]<-1;
        newConstraintsVect[-(1:noSParam)]<-(newConstraintsVect[-(1:noSParam)]+1)*(newConstraintsVect[-(1:noSParam)]>0);
      }

      if(is.list(nbdadata)){
        nbdadataTemp<-list()
        for(dataset in 1:length(nbdadata)){
          nbdadataTemp<-c(nbdadataTemp,constrainedNBDAdata(nbdadata=nbdadata[[dataset]],constraintsVect=newConstraintsVect,offsetVect=newOffsetVect))
        }
      }else{
        nbdadataTemp<-constrainedNBDAdata(nbdadata=nbdadata,constraintsVect=newConstraintsVect,offsetVect=newOffsetVect)
      }

      lowerLimModel<-tadaFit(nbdadataTemp,type=type,baseline=baseline)
      if(lowerLimModel@varNames[1]=="No variables"){
        propST<-oadaPropSolveByST(par=0,nbdadata=lowerLimModel@nbdadata,exclude.innovations=exclude.innovations,innovations=innovations)
      }else{
        propST<-oadaPropSolveByST(model=lowerLimModel,exclude.innovations=exclude.innovations,innovations=innovations)
      }
      lowerLimitPropST[i]<-propST[length(propST)]
    }

    lowerLimitRecordFull<-lowerLimitPropSTFull<-rep(NA,dim(aicTable@printTable)[1])
    lowerLimitRecordFull[modelsIncluded]<-lowerLimitRecord
    lowerLimitPropSTFull[modelsIncluded]<-lowerLimitPropST

    output1<-data.frame(model=aicTable@printTable$model,netCombo=aicTable@printTable$netCombo,lowerCI=lowerLimitRecordFull,propST=lowerLimitPropSTFull,deltaAICc=aicTable@printTable$deltaAICc,
                        akaikeWeight=aicTable@printTable$AkaikeWeight
    )[!is.na(lowerLimitRecordFull),]

    output<-cbind(output1,adjAkWeight=output1$akaikeWeight/sum(output1$akaikeWeight),
                  cumulAdjAkWeight=cumsum(output1$akaikeWeight/sum(output1$akaikeWeight)))

    return(output)
  }


}

#'Test the sensitivity of the lower limit of an s parameter to model selection uncertainty
#'
#'Runs \code{\link{multiModelLowerLimits}} across multiple computer cores. see documentation for
#'\code{\link{multiModelLowerLimits}} for details.
#'
#'@seealso \code{\link{multiModelLowerLimits}}
#'
#'@param cores numerical giving the number of computer cores to be used in parallel to fit the models in the set, thus speeding
#'up the process. By default set to 2. For a standard desktop computer at the time of writing 4-6 is advised.
#'@inheritParams oadaAICtable
#'@inheritParams profLikCI
#'@inheritParams nbdaPropSolveByST
#'@inheritParams multiModelLowerLimits
#'
#'@return data.frame. See \code{\link{multiModelLowerLimits}} for details.


multiModelLowerLimits_multicore<-function(which,aicTable,cores=2,deltaThreshold=Inf,conf=0.95,
                                          modelIndex=NULL,searchRange=NULL,exclude.innovations=T,innovations=NULL,startValue=NULL,
                                          lowerList=NULL,upperList=NULL,method="nlminb", gradient=T,iterations=150){

  noModels<-dim(aicTable@printTable)[1];
  #If cores is more than the number of models, reduce so two models are fit per core
  if(cores>noModels) cores<-noModels%/%2
  #How many models will be run on each core?
  modelsPerCore<-noModels%/%cores
  #How many will be left over to run at the end?
  remainderModels<-noModels%%cores

  #Do the remainder set first
  if(remainderModels>0){
    output<-multiModelLowerLimits(which=which,aicTable=aicTable,deltaThreshold=deltaThreshold,conf=conf,modelIndex=1:remainderModels,searchRange=searchRange,exclude.innovations=exclude.innovations,innovations=innovations,
                                  startValue=startValue,lowerList=lowerList,upperList=upperList,
                                            method=method, gradient=gradient,iterations=iterations)
  }else{output<-NULL}

  #Set up for parallel processing as detailed in http://www.parallelr.com/r-with-parallel-computing/
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)

    fromCores <- foreach(i=1:cores,.combine=rbind) %dopar%
      {
        #I think we need to reload the NBDA package into each thread
        library(NBDA)
        library(coxme)
        #Identify which models need to be fitted by this core
        modelSet<-(1:modelsPerCore)+modelsPerCore*(i-1)+remainderModels

        output<-multiModelLowerLimits(which=which,aicTable=aicTable,deltaThreshold=deltaThreshold,conf=conf,modelIndex=modelSet,searchRange=searchRange,exclude.innovations=exclude.innovations,innovations=innovations,
                                      startValue=startValue,lowerList=lowerList,upperList=upperList,
                                      method=method, gradient=gradient,iterations=iterations)

        output
      }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the table for the multiple cores with the remainder set fitted first

    output<-rbind(output,fromCores)
    return(output)

}


#'Obtain estimates of the proportion of events that occured by social transmission for a set of models
#'
#'Takes a set of models fitted by \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}} and calculates the estimated
#'proportion of events occuring by social transmssion via each network (propST) in each
#'model using \code{\link{nbdaPropSolveByST}}. This, alongside the \code{\link{multiModelLowerLimits}} function, allows the user to test the robustness of conclusions
#'about social transmission to model selection uncertainty. The function also calculates a model averaged estimate of propST for
#'each network.
#'
#'
#'@seealso \code{\link{multiModelLowerLimits}}, \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}},
#'\code{\link{plotProfLik}}, \code{\link{nbdaPropSolveByST}}, \code{\link{multiModelLowerLimits}}.
#'
#'@param aicTable an object of class \code{\link{oadaAICtable}} or \code{\link{tadaAICtable}}.
#'@param subset optional numerical vector specifying a subset of model numbers to be included.
#'@param statusBar logical indictaing whether a status/ progress bar should be produced
#'@return list with the following components:\describe{
#'   \item{$propSTtable}{propST via each network for each model in the set}
#'   \item{$modelAverages}{model averaged propST for each network}
#'   }


multiModelPropST<-function(aicTable,subset=NULL,statusBar=T){

  if(aicTable@nbdadata[[1]]@multi_ilv[1]=="ILVabsent"){
    noMultiILVs<-0
  }else{
    noMultiILVs<-length(aicTable@nbdadata[[1]]@multi_ilv)
  }
  coefficients<-cbind(aicTable@MLEs,aicTable@MLEilv,aicTable@MLEint,matrix(0,ncol=noMultiILVs,nrow=dim(aicTable@MLEs)[1]))

  if(is.null(subset)) subset<-1:dim(coefficients)[1]
  output<-NULL
  if(statusBar) pb <- txtProgressBar(min=0, max=length(subset), style=3)
  for(i in subset){
    output<-rbind(output,nbdaPropSolveByST(coefficients[i,],nbdadata =aicTable@nbdadata))
    if(statusBar)setTxtProgressBar(pb, i)
  }
  if(statusBar)close(pb)

  newAkaikeWeight<-aicTable@AkaikeWeight[subset]/sum(aicTable@AkaikeWeight[subset])

  modelAverages<-apply(output*newAkaikeWeight,2,sum)

  propSTtable<-cbind(model=subset,output,weight=newAkaikeWeight)
  propSTtable<-propSTtable[order(-newAkaikeWeight),]

  return(list(propSTtable=propSTtable,modelAverages=modelAverages))
}

