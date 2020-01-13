#At a future date: Make it possible to use multiple Cores on this. Need to make sure the temp objects are being written to different environments but should work automatically.

#I now have a combineOadaAICtable function, I need to write an equivalent one for TADA

#Define class of object for the fitted additive model
setClass("tadaAICtable",representation(nbdaMultiDiff
="character",nbdadata="list",convergence="logical",loglik="numeric",aic="numeric",aicc="numeric",constraintsVectMatrix="matrix", offsetVectMatrix="matrix", MLEs="matrix",SEs="matrix",MLEilv="matrix",SEilv="matrix",MLEint="matrix",SEint="matrix",MLEhaz="matrix",SEhaz="matrix",typeVect="character",baselineVect="character",deltaAIC="numeric",RelSupport="numeric",AkaikeWeight="numeric",printTable="data.frame"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
    signature(.Object = "tadaAICtable"),
    function (.Object, nbdadata,typeVect,baselineVect,constraintsVectMatrix,offsetVectMatrix,noHazFunctParsCustom,hazFunct,cumHaz,startValue,method,gradient,iterations,aicUse,lowerList,writeProgressFile,combineTables=F,
              MLEs,SEs,MLEilv,SEilv,MLEint,SEint,MLEhaz,SEhaz,
              convergence,loglik,aic,aicc,netComboModifierVect,statusBar,...)
          {


    if(is.null(typeVect)){typeVect<-rep("social",dim(constraintsVectMatrix)[1])}
    if(is.null(baselineVect)){baselineVect<-rep("constant",dim(constraintsVectMatrix)[1])}

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

      if(is.na(nbdadataTemp1@TADAtime1[1])){
        print("TADA times absent from nbdaData object")
        return(NULL)
      }

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
		  if(statusBar)  pb <- txtProgressBar(min=0, max=noModels, style=3)

		noHazFunctParsVect<-rep(NA,length(baselineVect))
		noHazFunctParsVect[baselineVect=="constant"]<-1
		noHazFunctParsVect[baselineVect=="weibull"]<-2
		noHazFunctParsVect[baselineVect=="gamma"]<-2
		noHazFunctParsVect[baselineVect=="custom"]<-noHazFunctParsCustom

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
		MLEhaz<-matrix(NA,nrow=noModels,ncol=max(noHazFunctParsVect),dimnames=list(1:noModels, paste("Baseline parameter",1:max(noHazFunctParsVect),sep="")))
		SEhaz<-matrix(NA,nrow=noModels,ncol=max(noHazFunctParsVect),dimnames=list(1:noModels, paste("SE Baseline parameter",1:max(noHazFunctParsVect),sep="")))

		#Set up various vectors to record things about each model
		convergence<-loglik<-aic<-aicc<-seApprox<-rep(NA,noModels)

		#Loop through the rows of the constrainstsVectMatrix creating the constrained objects and thus fitting the specified model each time
		for (i in 1:noModels){

		  #Update progress bar
		  if(statusBar) setTxtProgressBar(pb, i)
		  #Write file to working directory saying what model we are on
      if(writeProgressFile){write.csv(paste("Currently fitting model",i, "out of", noModels),file=paste("oadaTableProgressFile",nbdadataTemp1@label[1],".txt",sep=""),row.names =F)}


		  constraintsVect<-constraintsVectMatrix[i,]
		  offsetVect <-offsetVectMatrix[i,]


		  if(is.null(startValue)) {
		    newStartValue<-NULL
		  }else{
		    newStartValue<-c(startValue[1:noHazFunctParsVect[i]],tapply(startValue[-(1:max(noHazFunctParsVect))],constraintsVect,mean)[-1])
		  }

		  if(is.null(lowerList)) {
		    lower<-NULL
		  }else{
		    lower<-lowerList[i,]
		    lower<-lower[constraintsVect!=0]
		  }

		  #If the user has specified all zeroes for the s parameters, we need to change it to an "asocial" type
		  #And we need to add a one for the first s parameter so the constrained NBDA object can be created
		  #And the ILV numbers need shifting up one, to be shifted down later
		  if(sum(constraintsVect[1:noSParam])==0){
		    typeVect[i]<-"asocial";
		    constraintsVect[1]<-1;
		    constraintsVect[-(1:noSParam)]<-(constraintsVect[-(1:noSParam)]+1)*(constraintsVect[-(1:noSParam)]>0);
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
			try(model<-tadaFit(nbdadata= nbdadataTemp,type=typeVect[i],startValue=newStartValue,method=method,gradient=gradient,iterations=iterations,standardErrors=T,baseline=baselineVect[i],noHazFunctPars=c,hazFunct=hazFunct,cumHaz=cumHaz))
      if(!is.null(model)){

      #Record baseline hazard rate parameters, then remove them from the model object so the OADA code works
      MLEhaz[i,1:noHazFunctParsVect[i]]<-model@outputPar[1:noHazFunctParsVect[i]]
      SEhaz[i,1:noHazFunctParsVect[i]]<-model@se[1:noHazFunctParsVect[i]]
      model@outputPar<-model@outputPar[-(1:noHazFunctParsVect[i])]
      model@se<-model@se[-(1:noHazFunctParsVect[i])]


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

			#Record MLE and SE for the  effect of  ILVs on social learning
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
		if(statusBar) close(pb)
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
		  printTable<-data.frame(model=1:noModels,type=newType,netCombo=netCombo,baseline=baselineVect,constraintsVectMatrix, offsetVectMatrix,convergence,loglik,MLEhaz,MLEs,MLEilv,MLEint,SEhaz,SEs,SEilv,SEint,
		                         aic,aicc,deltaAIC,RelSupport,AkaikeWeight)
		  printTable <-printTable[order(aic),]
		}else{
		  printTable<-data.frame(model=1:noModels,type=newType,netCombo=netCombo,baseline=baselineVect,constraintsVectMatrix, offsetVectMatrix,convergence,loglik,MLEhaz,MLEs,MLEilv,MLEint,SEhaz,SEs,SEilv,SEint,
		                         aic,aicc, deltaAICc=deltaAIC,RelSupport,AkaikeWeight)
		  printTable <-printTable[order(aicc),]
		}


		if(is.list(nbdadata)){
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata,convergence= convergence, loglik= as.numeric(loglik),aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 baselineVect=baselineVect, MLEhaz= MLEhaz,SEhaz= SEhaz,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)
		}else{
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), convergence= convergence, loglik= as.numeric(loglik),aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 baselineVect=baselineVect, MLEhaz= MLEhaz,SEhaz= SEhaz,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)

		}
    }
)


#'Fit a set of TADA models for multi-model inference
#'
#'\code{tadaAICtable} takes diffusion data in the form of an nbdaData object (\code{\link{nbdaData}}) or (\code{\link{dTADAData}})
#' or a list of nbdaData/ dTADAData objects (for multiple diffusions). It then fits a set of models using \code{\link{tadaFit}}
#' and return them in an object of class \code{tadaAICtable}. Arguments not listed below are used internally by
#' \code{\link{combineTadaAICtables}} when making calls to \code{tadaAICtable} and can be ignored by the user.
#'
#'Each row of \code{constraintsVectMatrix}, \code{offsetVectMatrix}, and \code{baselineVect} determines a model to
#'be fitted. For each row, constrained \code{nbdaData} objects are created using (\code{\link{constrainedNBDAdata}}) with
#'\code{constraintsVect=constraintsVectMatrix[i,]} and \code{offsetVect=offsetVectMatrix[i,]}. A model is then fitted to
#'the \code{nbdaData} object(s) using \code{\link{tadaFit}} with a baseline determined by \code{offsetVect=baselineVect[i]}.
#'
#'@seealso \code{\link{networksSupport}}, \code{\link{typeByNetworksSupport}}, \code{\link{modelAverageEstimates}},
#' \code{\link{variableSupport}}, \code{\link{unconditionalStdErr}},\code{\link{baselineSupport}},
#' \code{\link{combineTadaAICtables}}. For OADA models use \code{\link{oadaAICtable}}.
#'
#'@param nbdadata an object of class (\code{\link{nbdaData}}) or (\code{\link{dTADAData}}) to fit models to a single diffusion or a list of
#'nbdaData/dTADAData objects to fit a model to multiple diffusions.
#'@param constraintsVectMatrix a numerical matrix specifying the constraints, with each row specifiying a model to be fitted.
#'The number of columns is equal to the number of parameters in the (\code{\link{nbdaData}}) or (\code{\link{dTADAData}}) object(s) input to argument
#'\code{nbdadata}. Each row then specifies the \code{constraintsVect} to be passed to the \code{\link{constrainedNBDAdata}}
#' when creating the data object(s) to which that model is fitted.
#'@param typeVect optional character vector specifying if each model is "asocial" or "social". However, it is not usually
#'necessary to specify, since models with all s parameters constrained =0 are automatically classified as "asocial", and
#'others are assumed to be "social".
#'@param baselineVect optional character vector specifying the baseline function for each model (see \code{\link{tadaFit}}).
#'@param offsetVectMatrix an optional numerical matrix specifying the offsets, with each row specifiying the offsets for each
#'model to be fitted. The number of columns is equal to the number of parameters in the (\code{\link{nbdaData}}) or (\code{\link{dTADAData}}) object(s)
#'input to argument \code{nbdadata}. Each row then specifies the \code{offsetVect} to be passed to the
#'\code{\link{constrainedNBDAdata}} when creating the data object(s) to which that model is fitted.
#'@param cores numerical giving the number of computer cores to be used in parallel to fit the models in the set, thus
#'speeding up the process. By default set to 1. For a standard desktop computer at the time of writing 4-6 is advised.
#'@param modelsPerCorePerSet  optional numerical. If specified the models can be fit in sets, and a progress file written
#'after each set is completed. This means progress is not completely lost in the case of a crash/ powercut etc. For example,
#'if we have 400 models, we can specify \code{cores=4} and \code{modelsPerCorePerSet=10}. This means 10 models are fitted on
#'each core, then progress is saved with the first 40 models, then the next 40 and so on.
#'@param writeProgressFile logical. If set to T, a file is written to the working directory when each set of models have
#'been completed with the tadaAICtable for the models fitted so far. In the event of a crash, the remining models can be
#'fitted as a separate set, then combined using \code{\link{combineTadaAICtables}}.
#'@param statusBar optional logical. Status bar only works when \code{cores=1}.
#'@param noHazFunctParsCustom  optional numerical necessary if "custom" is specified for any models in \code{baselineVect}.
#'See \code{\link{tadaFit}} for details.
#'@param hazFunct  optional function necessary if "custom" is specified for any models in \code{baselineVect}.
#'See \code{\link{tadaFit}} for details.
#'@param cumHaz  optional function necessary if "custom" is specified for any models in \code{baselineVect}.
#'See \code{\link{tadaFit}} for details.
#'@param startValue optional numeric vector giving start values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the full model.
#'@param lowerList optional numeric matrix giving lower values for the maximum likelihood optimization for each model.
#'Columns to match the number of parameters fitted in the full model, rows matched to the number of models. Can be used if
#'some models have convergence problems or trigger errors.
#'@param upperList optional numeric matrix giving upper values for the maximum likelihood optimization for each model.
#'Columns to match the number of parameters fitted in the full model, rows matched to the number of models. Can be used if
#'some models have convergence problems or trigger errors.
#'@param method optional character string passed to \code{\link{tadaFit}}.
#'@param gradient optional logical passed to \code{\link{tadaFit}}.
#'@param iterations optional numerical passed to \code{\link{tadaFit}}.
#'@param aicUse string specifying whether to use "aicc" or "aic".
#'@param combineTables logical used internally by \code{\link{combineTadaAICtables}} when making calls to \code{tadaAICtable}.
#'
#'@return An object of class \code{tadaAICtable}.
#'@section print(tadaAICtable)  components:
#'A data.frame giving a summary of models ordered by AIC can be obtained using \code{print(<name of tadaAICtable>)}. This has
#'the following columns, listed in order: \describe{
#'   \item{model}{Model number, i.e. the row of \code{constraintsVectMatrix} used to generate the model.}
#'   \item{type}{Type of model. noILVs, additive (ILV effects on asocial learning only), multiplicative (all ILVs have same
#'   effect on asocial and social learning), unconstrained (differing effects on asocial and social learning for at least one
#'   ILV), asocial}
#'   \item{netcombo}{A representaion of the network effects present in the model, i.e. the constraints on the s
#'   parameters. See \code{\link{constrainedNBDAdata}}.}
#'   \item{baseline}{The baseline function used for each model.}
#'   \item{CONS.}{The constraint on each parameter, as taken from \code{constraintsVectMatrix}}
#'   \item{OFF}{The offset on each parameter, as taken from \code{offsetVectMatrix}}
#'   \item{convergence}{Was convergence reported by the optimization algorithm?}
#'   \item{loglik}{-log-likelihood for the model.}
#'   \item{s...}{Maximum likelihood estimates for s parameters.}
#'   \item{ASOCIAL...}{Maximum likelihood estimates for effects of ILVs on asocial learning.}
#'   \item{SOCIAL...}{Maximum likelihood estimates for effects of ILVs on social learning.}
#'   \item{ASOCIAL:SOCIAL...}{Maximum likelihood estimates for multiplicative effects of ILVs see \code{\link{tadaFit}}.}
#'   \item{SE...}{Standard errors for each parameter, set to 0 when a parameter was constrained. See \code{\link{tadaFit}}.}
#'   \item{aic}{AIC for the model.}
#'   \item{aicc}{AICc for the model.}
#'   \item{deltaAICc}{Difference in AICc or AIC from the best model.}
#'   \item{RelSupport}{Relative support for the model compared to the best model, calculated as exp(-0.5*deltaAICc).}
#'   \item{AkaikeWeight}{Akaike weight for the model. Can be interpretted as the probability that model has the highest
#'   predictive power (K-L information) out of the set of models considered.}
#'   }
#'@section tadaAICtable components:\describe{
#'   \item{@@nbdadata}{The unconstrained data the model is fitted to, as a list of nbdaData objects.}
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
#'   \item{@@MLEhaz}{Maximum likelihood estimates for baseline rate/hazard function (ordered by
#'   \code{constraintsVectMatrix}).}
#'   \item{@@SEhaz}{Standard errors for baseline rate/hazard function (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@typeVect}{Type of models (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@baselineVect}{The baseline function used for each model. (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@deltaAICc}{Difference in AICc or AIC from the best model. (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@RelSupport}{Relative support for the model compared to the best model, calculated as exp(-0.5*deltaAICc).
#'   (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@RelSupport}{Akaike weight for the model. (ordered by \code{constraintsVectMatrix}).}
#'   \item{@@printTable}{data.frame to be output by the print method for \code{tadaAICtable} (see above).}
#'}


#Function for implementing the initialization
tadaAICtable <-function(nbdadata,  constraintsVectMatrix,typeVect=NULL,baselineVect=NULL, offsetVectMatrix = NULL,
                        cores=1, modelsPerCorePerSet=NULL,writeProgressFile=F,statusBar=NULL,
                        noHazFunctParsCustom=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL),
                        startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,combineTables=F,
                        MLEs=NULL,SEs=NULL,MLEilv=NULL,SEilv=NULL,MLEint=NULL,SEint=NULL, MLEhaz=NULL,SEhaz=NULL,
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
    tadaAICtable_multiCore(nbdadata=nbdadata,constraintsVectMatrix=constraintsVectMatrix,cores=cores,typeVect=typeVect, offsetVectMatrix = offsetVectMatrix,baselineVect=baselineVect,
                           modelsPerCorePerSet=modelsPerCorePerSet,writeProgressFile=writeProgressFile,statusBar=statusBar,
                           startValue=startValue,method=method, gradient=gradient,iterations=iterations,
                           aicUse=aicUse,lowerList=lowerList,
                           noHazFunctParsCustom=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL))

  }else{
    if(is.null(statusBar))statusBar<-T
	  return(new("tadaAICtable",nbdadata= nbdadata, typeVect= typeVect, constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix = offsetVectMatrix, baselineVect=baselineVect,noHazFunctParsCustom=noHazFunctParsCustom,
	           hazFunct=hazFunct,cumHaz=cumHaz,startValue= startValue,method= method, gradient= gradient,
	           iterations= iterations,aicUse= aicUse,lowerList=lowerList,writeProgressFile=writeProgressFile,statusBar=statusBar, combineTables=combineTables,
	           MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint,MLEhaz=MLEhaz,SEhaz=SEhaz,
	           convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect))
  }
}

#Method for initializing addFit object- including model fitting
print.tadaAICtable<-function (tadaAICtable)
    {
		tadaAICtable@printTable
	}

#'Get the support for each baseline rate function from a tadaAICtable
#'
#'Calculates the support for each baseline rate function fitted in a set of models fitted using #'\code{\link{tadaAICtable}}.
#'
#'
#'@seealso \code{\link{oadaAICtable}}, \code{\link{tadaAICtable}}, \code{\link{typeSupport}},
#'\code{\link{networksSupport}}, \code{\link{modelAverageEstimates}}, \code{\link{variableSupport}},
#'\code{\link{unconditionalStdErr}}.
#'
#'@param nbdaAICtable an object of class \code{\link{tadaAICtable}}.
#'
#'@return matrix giving the support (total Akaike Weight) and number of models fitted for each combination.


baselineSupport<-function(tadaAICtable){
  #Calculate support for each combination of network constraints in the table
  support<-tapply(tadaAICtable@printTable$AkaikeWeight, tadaAICtable@printTable$baseline,sum)
  numbers<-tapply(tadaAICtable@printTable$AkaikeWeight, tadaAICtable@printTable$baseline,length)
  return(data.frame(support=support,numberOfModels=numbers))
}


#'Combine two or more tadaAICtables into a single tadaAICtable
#'
#'Takes two or more \code{\link{tadaAICtable}} objects, containing different models including the same number of networks and
#'the same ILVs, and combine them into a single table.
#'
#'This function can be used for a variety of practical reasons. If the \code{\link{tadaAICtable}} is interrupted and a partial
#'tadaAICtable recovered, the remaining models can be run as a separate set and combined with the first partial set. Alternatively,
#'the user might want to add more models to a set without re-running the entire set, or run different parts of a model set on
#'different computers to speed up computation time. It is also possible to run model sets with different networks and then combine
#'them using \code{combineTadaAICtables}, so long as the number of networks matches (though it is possible to have dummy networks
#'that are constrained to have s=0 in all models, in order to match network number). In such cases a modifier can be added to the
#'netcombo codes from each tadaAICtable so they can be distinguished by functions such as \code{\link{networksSupport}}. e.g.
#'if we have two tadaAICtable objects each with the netcombos "1:0","0:1","1:2". If we use
#'\code{netComboModifier= c("NetworksA:","NetworksB:")} we get netcombos "NetworksA:1:0","NetworksA:0:1","NetworksA:1:2",
#'"NetworksB:1:0","NetworksB:0:1","NetworksB:1:2" in the resulting tadaAICtable object.
#'
#'@seealso \code{\link{tadaAICtable}}
#'
#'@param tadaAICtableList a lift of \code{\link{tadaAICtable}} objects to be combined into a single table.
#'@param aicUse string specifying whether to use "aicc" or "aic".
#'@param netComboModifier optional character vector with length matching the length of \code{tadaAICtableList} to modify the
#'netcombo strings recorded from each individual \code{\link{tadaAICtable}} object (see below).
#'@return An object of class \code{tadaAICtable}.

combineTadaAICtables<-function(tadaAICtableList,aicUse="aicc",netComboModifier=rep("",length(tadaAICtableList))){

  i<-1
  tadaAICtableTemp<-tadaAICtableList[[i]]
  nbdadata<-tadaAICtableTemp@nbdadata

  writeProgressFile=F;
  startValue=NULL;
  method=NULL;
  gradient=NULL;
  iterations=NULL;
  lowerList=NULL;

  typeVect<-tadaAICtableTemp@typeVect;
  baselineVect<-tadaAICtableTemp@baselineVect
  constraintsVectMatrix<-tadaAICtableTemp@constraintsVectMatrix;
  offsetVectMatrix<-tadaAICtableTemp@offsetVectMatrix;
  netComboModifierVect<-rep(netComboModifier[i],dim(tadaAICtableTemp@constraintsVectMatrix)[1])
  maxNoHazFunctParsVect<-dim(tadaAICtableTemp@MLEhaz)[2]

  #Loop through the tadaAICtableList and create combined versions of typeVect, constraintsVectMatrix, and offsetVectMatrix
  for(i in 2:length(tadaAICtableList)){
    #Read in the next tadaAICtable
    tadaAICtableTemp<-tadaAICtableList[[i]]
    typeVect<-c(typeVect,tadaAICtableTemp@typeVect)
    baselineVect<-c(baselineVect,tadaAICtableTemp@baselineVect)
    constraintsVectMatrix<-rbind(constraintsVectMatrix,tadaAICtableTemp@constraintsVectMatrix)
    offsetVectMatrix<-rbind(offsetVectMatrix,tadaAICtableTemp@offsetVectMatrix)
    netComboModifierVect<-c(netComboModifierVect,rep(netComboModifier[i],dim(tadaAICtableTemp@constraintsVectMatrix)[1]))

    maxNoHazFunctParsVect<-max(maxNoHazFunctParsVect,dim(tadaAICtableTemp@MLEhaz)[2])


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
  MLEhaz<-matrix(NA,nrow=noModels,ncol=maxNoHazFunctParsVect,dimnames=list(1:noModels, paste("Baseline parameter",1:maxNoHazFunctParsVect,sep="")))
  SEhaz<-matrix(NA,nrow=noModels,ncol=maxNoHazFunctParsVect,dimnames=list(1:noModels, paste("SE Baseline parameter",1:maxNoHazFunctParsVect,sep="")))

  #Set up various vectors to record things about each model
  convergence<-loglik<-aic<-aicc<-seApprox<-rep(NA,noModels)

  #Loop through the tadaAICtableList and fill in the various matrices and vectors
  #Track the start point to fill in each time
  startPoint<-1
  for(i in 1:length(tadaAICtableList)){
    #Read in the next tadaAICtable
    tadaAICtableTemp<-tadaAICtableList[[i]]
    endPoint<-dim(tadaAICtableTemp@constraintsVectMatrix)[1]+startPoint-1

    MLEs[startPoint:endPoint,]<-tadaAICtableTemp@MLEs
    SEs[startPoint:endPoint,]<-tadaAICtableTemp@SEs
    MLEhaz[startPoint:endPoint,1:dim(tadaAICtableTemp@MLEhaz)[2]]<-tadaAICtableTemp@MLEhaz
    SEhaz[startPoint:endPoint,1:dim(tadaAICtableTemp@SEhaz)[2]]<-tadaAICtableTemp@SEhaz
    if(noILVasoc!=0){
      MLEilv[startPoint:endPoint,]<-tadaAICtableTemp@MLEilv
      SEilv[startPoint:endPoint,]<-tadaAICtableTemp@SEilv
      MLEint[startPoint:endPoint,]<-tadaAICtableTemp@MLEint
      SEint[startPoint:endPoint,]<-tadaAICtableTemp@SEint
    }
    convergence[startPoint:endPoint]<-tadaAICtableTemp@convergence
    loglik[startPoint:endPoint]<-tadaAICtableTemp@loglik
    aic[startPoint:endPoint]<-tadaAICtableTemp@aic
    aicc[startPoint:endPoint]<-tadaAICtableTemp@aicc
    startPoint<-endPoint+1
  }

  tableTemp<-tadaAICtable(nbdadata=nbdadata,  constraintsVectMatrix=constraintsVectMatrix,typeVect=typeVect, baselineVect=baselineVect,offsetVectMatrix = offsetVectMatrix,
                          startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,
                          combineTables=T,
                          MLEs=MLEs,SEs=SEs,MLEilv=MLEilv,SEilv=SEilv,MLEint=MLEint,SEint=SEint, MLEhaz=MLEhaz, SEhaz=SEhaz,
                          convergence=convergence,loglik=loglik,aic=aic,aicc=aicc,netComboModifierVect=netComboModifierVect)

  return(tableTemp)

}


#This function runs an tadaAICtable with multiple cores, but is called via tadaAICtable with the cores argument
#so does not need calling by the user
tadaAICtable_multiCore<-function(nbdadata,constraintsVectMatrix,cores,typeVect=NULL, offsetVectMatrix = NULL, baselineVect=NULL,modelsPerCorePerSet=NULL,writeProgressFile=F,
                                 statusBar=F,startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,
                                 noHazFunctParsCustom=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL)){
  noModels<-dim(constraintsVectMatrix)[1];
  #If cores is more than the number of models, reduce so two models are fit per core
  if(cores>noModels) cores<-(noModels*2)
  #How many models will be run on each core?
  modelsPerCore<-noModels%/%cores
  #How many will be left over to run at the end?
  remainderModels<-noModels%%cores
  #aicTable funcions do not work with a single model, so adjust to make sure a single model is not fitted


  #I will set the models to run in sets of 10 per core (default) then write the tadaAICtable so far to a
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
    if(is.null(baselineVect)){baselineVectTemp<-NULL}else{baselineVectTemp<-baselineVect[1:remainderModels]}
    if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[1:remainderModels,]}

    cumulativeAICtable<-tadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,baselineVect=baselineVectTemp,
                                     startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,statusBar = F,
                                     noHazFunctParsCustom=noHazFunctParsCustom,hazFunct=hazFunct,cumHaz=cumHaz)
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
      library(NBDA)
      #Identify which models need to be fitted in this core x set combination
      modelSet<-((set-1)*modelsPerCorePerSet*cores + (i-1)*modelsPerCorePerSet+remainderModels)+1:modelsPerCorePerSet
      #Cut down constraintsVectMatrix and, if necessary typeVect and offsetVectMatrix
      constraintsVectMatrixTemp<- rbind(constraintsVectMatrix[modelSet,])
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(baselineVect)){baselineVectTemp<-NULL}else{baselineVectTemp<-baselineVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-tadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,baselineVect=baselineVectTemp,
                                 startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,
                                 noHazFunctParsCustom=noHazFunctParsCustom,hazFunct=hazFunct,cumHaz=cumHaz)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative tadaAICtable so far
    cumulativeAICtable<-combineTadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
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
      constraintsVectMatrixTemp<- rbind(constraintsVectMatrix[modelSet,])
      if(is.null(typeVect)){typeVectTemp<-NULL}else{typeVectTemp<-typeVect[modelSet]}
      if(is.null(baselineVect)){baselineVectTemp<-NULL}else{baselineVectTemp<-baselineVect[modelSet]}
      if(is.null(offsetVectMatrix)){offsetVectMatrixTemp<-NULL}else{offsetVectMatrixTemp<-offsetVectMatrix[modelSet,]}
      #Fit the set of models required by this core x set
      tempAICtable<-tadaAICtable(nbdadata,constraintsVectMatrixTemp,typeVect=typeVectTemp,offsetVectMatrix=offsetVectMatrixTemp,baselineVect=baselineVectTemp,
                                 startValue=startValue,method=method, gradient=gradient,iterations=iterations,aicUse=aicUse,lowerList=lowerList,writeProgressFile=F,
                                 noHazFunctParsCustom=noHazFunctParsCustom,hazFunct=hazFunct,cumHaz=cumHaz)
      #Return the table, which is combined into a list with the tables from the other cores
      tempAICtable
    }
    #Stop the cluster
    stopImplicitCluster()
    stopCluster(cl)
    #Combine the list of models with the cumumative tadaAICtable so far
    cumulativeAICtable<-combineTadaAICtables(c(cumulativeAICtable,tablesFromSet),aicUse=aicUse)
  }
  if(writeProgressFile) save(cumulativeAICtable,file="finalAICtable_parallel.Rdata")
  if(statusBar)close(pb)

  return(cumulativeAICtable)
}
