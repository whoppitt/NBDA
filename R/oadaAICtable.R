#At a future date: Make it possible to use multiple Cores on this. Need to make sure the temp objects are being written to different environments but should work automatically.
#Update- now added a function to combine oadaAICtables, so the user can at least run parts of the constraintsVectMatrix on different cores then manually
#combine them using this fuction
#I should also be able to use this to create a multicore version later.

#Define class of object for the fitted additive model
setClass("oadaAICtable",representation(nbdaMultiDiff
="character",nbdadata="list",convergence="logical",loglik="numeric",aic="numeric",aicc="numeric",constraintsVectMatrix="matrix", offsetVectMatrix="matrix",
MLEs="matrix",SEs="matrix",MLEilv="matrix",SEilv="matrix",MLEint="matrix",SEint="matrix",
typeVect="character",deltaAIC="numeric",RelSupport="numeric",AkaikeWeight="numeric",printTable="data.frame"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
    signature(.Object = "oadaAICtable"),
    function (.Object, nbdadata,typeVect,constraintsVectMatrix,offsetVectMatrix,startValue,method,gradient,iterations,aicUse,lowerList,writeProgressFile,combineTables=F,
              MLEs,SEs,MLEilv,SEilv,MLEint,SEint,
              convergence,loglik,aic,aicc,netComboModifierVect,statusBar,...)
    {

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
		convergence<-loglik<-aic<-aicc<-seApprox<-rep(NA,noModels)

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
		  model<-NULL
			try(model<-oadaFit(nbdadata= nbdadataTemp,type=typeVect[i],startValue=newStartValue,method=method,gradient=gradient,iterations=iterations))
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
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata,convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)
		}else{
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)

		}

    }
)



#Function for implementing the initialization
oadaAICtable <-function(nbdadata,  constraintsVectMatrix,typeVect=NULL, offsetVectMatrix = NULL, cores=1, modelsPerCorePerSet=NULL,writeProgressFile=F,statusBar=NULL,startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,combineTables=F,
                        MLEs=NULL,SEs=NULL,MLEilv=NULL,SEilv=NULL,MLEint=NULL,SEint=NULL,
                        convergence=NULL,loglik=NULL,aic=NULL,aicc=NULL,netComboModifierVect=""){

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
                           aicUse=aicUse,lowerList=lowerList)

  }else{
    if(is.null(statusBar))statusBar<-T
    return(new("oadaAICtable",nbdadata= nbdadata, typeVect= typeVect, constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix = offsetVectMatrix, startValue= startValue,method= method, gradient= gradient,
	           iterations= iterations,aicUse= aicUse,lowerList=lowerList,writeProgressFile=writeProgressFile,combineTables=combineTables,statusBar=statusBar,
	           MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint,
	           convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect))
  }
}

#Method for initializing addFit object- including model fitting
print.oadaAICtable<-function (oadaAICtable)
    {
		oadaAICtable@printTable
	}

typeSupport<-function(nbdaAICtable){
	#Calculate support for each type of model in the table
	support<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$type,sum)
	numbers<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$type,length)
	return(data.frame(support=support,numberOfModels=numbers))
}

networksSupport<-function(nbdaAICtable){
  #Calculate support for each combination of network constraints in the table
  support<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$netCombo,sum)
  numbers<-tapply(nbdaAICtable@printTable$AkaikeWeight, nbdaAICtable@printTable$netCombo,length)
  return(data.frame(support=support,numberOfModels=numbers))
}

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

#To be modified from model averaged estimates function
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

  typeVect<-oadaAICtableTemp@typeVect;
  constraintsVectMatrix<-oadaAICtableTemp@constraintsVectMatrix;
  offsetVectMatrix<-oadaAICtableTemp@offsetVectMatrix;
  netComboModifierVect<-rep(netComboModifier[i],dim(oadaAICtableTemp@constraintsVectMatrix)[1])

  #Loop through the oadaAICtableList and create combined versions of typeVect, constraintsVectMatrix, and offsetVectMatrix
  for(i in 2:length(oadaAICtableList)){
    #Read in the next oadaAICtable
    oadaAICtableTemp<-oadaAICtableList[[i]]
    typeVect<-c(typeVect,oadaAICtableTemp@typeVect)
    constraintsVectMatrix<-rbind(constraintsVectMatrix,oadaAICtableTemp@constraintsVectMatrix)
    offsetVectMatrix<-rbind(offsetVectMatrix,oadaAICtableTemp@offsetVectMatrix)
    netComboModifierVect<-c(netComboModifierVect,rep(netComboModifier[i],dim(oadaAICtableTemp@constraintsVectMatrix)[1]))

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
               startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,
               combineTables=T,
               MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint,
               convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect)

  return(tableTemp)

}

#This function runs an oadaAICtable with multiple cores, but is called via oadaAICtable with the cores argument
#so does not need calling by the user
oadaAICtable_multiCore<-function(nbdadata,constraintsVectMatrix,cores,typeVect=NULL, offsetVectMatrix = NULL, modelsPerCorePerSet=NULL,writeProgressFile=F,statusBar=F,startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL){
  noModels<-dim(constraintsVectMatrix)[1];
  #If cores is more than the number of models, reduce so one model is fit per cors
  if(cores>noModels) cores<-noModels
  #How many models will be run on each core?
  modelsPerCore<-noModels%/%cores
  #How many will be left over to run at the end?
  remainderModels<-noModels%%cores

  #I will set the models to run in sets of 10 per core (default) then write the oadaAICtable so far to a
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
                                   startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,statusBar = F)
  }else{cumulativeAICtable<-NULL}

  if(statusBar) pb <- txtProgressBar(min=0, max=numberOfInitialSets, style=3)

  #Loop through the initial sets
  for(set in 1:numberOfInitialSets){

    #Set up for parallel processing as detailed in http://www.parallelr.com/r-with-parallel-computing/
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)

    tablesFromSet <- foreach(i=1:cores) %dopar%
    {
      #I think we need to reload the NBDA package into each thread
     # library(NBDA)
      #Identify which models need to be fitted in this core x set combination
      modelSet<-((set-1)*modelsPerCorePerSet*cores + (i-1)*modelsPerCorePerSet+remainderModels)+1:modelsPerCorePerSet
      #Cut down constraintsVectMatrix and, if necessary typeVect and offsetVectMatrix
      constraintsVectMatrixTemp<- constraintsVectMatrix[modelSet,]
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-oadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,
                                       startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative oadaAICtable so far
    cumulativeAICtable<-combineOadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
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
      #Identify which models need to be fitted in this core x set combination
      modelSet<-(numberOfInitialSets*modelsPerCorePerSet*cores + (i-1)*remainderModelsPerCore+remainderModels)+1:remainderModelsPerCore
      #Cut down constraintsVectMatrix and, if necessary typeVect and offsetVectMatrix
      constraintsVectMatrixTemp<- constraintsVectMatrix[modelSet,]
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-oadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,
                                 startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative oadaAICtable so far
    cumulativeAICtable<-combineOadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
  }
  if(writeProgressFile) save(cumulativeAICtable,file="finalAICtable_parallel.Rdata")
  if(statusBar)close(pb)

  return(cumulativeAICtable)
}


#This function allows one to get the lower limits of the C.I. for a chosen s parameter for all models in the set containing that s parameter
#or within deltaThreshold AIC/AICc units. Works for OADA
multiModelLowerLimits<-function(which,aicTable,deltaThreshold=Inf,conf=0.95,exclude.innovations=T,innovations=NULL,startValue=NULL,lowerList=NULL,method="nlminb", gradient=T,iterations=150){
  if(class(aicTable)=="oadaAICtable"){
  #get the model set in which AIC or AICc was < the threshold and in which the parameter was present
  modelsIncluded<-aicTable@printTable$deltaAIC<deltaThreshold&aicTable@printTable[,which+4]
  modelSet<-aicTable@printTable$model[modelsIncluded]

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
      lowerLimitRecord[i]<-profLikCI(which=constraintsVect[which],model=modelFit,conf=conf,lowerRange=c(0,modelFit@outputPar[constraintsVect[which]]))[1]
    }


    #Fit lower limit model

    #Create the necessary constrained data objects
    newConstraintsVect<-constraintsVect
    newConstraintsVect[constraintsVect[which]]<-0
    newConstraintsVect[newConstraintsVect>0]<-as.numeric(factor(newConstraintsVect[newConstraintsVect>0]))
    newOffsetVect<-offsetVect
    newOffsetVect[constraintsVect[which]]<-newOffsetVect[constraintsVect[which]]+lowerLimitRecord[i]

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
      newConstraintsVect[constraintsVect[which]]<-0
      newConstraintsVect[newConstraintsVect>0]<-as.numeric(factor(newConstraintsVect[newConstraintsVect>0]))
      newOffsetVect<-offsetVect
      newOffsetVect[constraintsVect[which]]<-newOffsetVect[constraintsVect[which]]+lowerLimitRecord[i]

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


#This function takes the output of multiModelLowerLimits and returns the minimunm lowerCI and propST within a specified confidence set of models
confSetLowerLimit<-function(multiModelLowerLimitsSet,confSet=.95){
  if(confSet>=1){
    newSet<-multiModelLowerLimitsSet
  }else{
    newLimit<-min(multiModelLowerLimitsSet$cumulAdjAkWeight[multiModelLowerLimitsSet$cumulAdjAkWeight>=confSet])
    newSet<-multiModelLowerLimitsSet[multiModelLowerLimitsSet$cumulAdjAkWeight<newLimit,]
  }
  lowerCI<-min(newSet$lowerCI)
  propST<-min(newSet$propST)
  list(lowerCI=lowerCI,propST=propST)
}


