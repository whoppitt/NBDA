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
          newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[1]))))
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
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata,convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 baselineVect=baselineVect, MLEhaz= MLEhaz,SEhaz= SEhaz,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)
		}else{
		  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), convergence= convergence, loglik= loglik,aic= aic,aicc= aicc,constraintsVectMatrix= constraintsVectMatrix, offsetVectMatrix= offsetVectMatrix,
		                 baselineVect=baselineVect, MLEhaz= MLEhaz,SEhaz= SEhaz,
		                 MLEs= MLEs,SEs= SEs,MLEilv= MLEilv,SEilv= SEilv,MLEint= MLEint,SEint= SEint,
		                 typeVect= newType,deltaAIC= deltaAIC,RelSupport= RelSupport,AkaikeWeight= AkaikeWeight,printTable=printTable)

		}
    }
)



#Function for implementing the initialization
tadaAICtable <-function(nbdadata,  constraintsVectMatrix,typeVect=NULL,baselineVect=NULL, offsetVectMatrix = NULL,
                        cores=1, modelsPerCorePerSet=NULL,writeProgressFile=F,statusBar=NULL,
                        noHazFunctParsCustom=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL),
                        startValue=NULL,method="nlminb", gradient=T,iterations=150,aicUse="aicc",lowerList=NULL,combineTables=F,
                        MLEs=NULL,SEs=NULL,MLEilv=NULL,SEilv=NULL,MLEint=NULL,SEint=NULL, MLEhaz=NULL,SEhaz=NULL,
                        convergence=NULL,loglik=NULL,aic=NULL,aicc=NULL,netComboModifierVect=""){
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


baselineSupport<-function(tadaAICtable){
  #Calculate support for each combination of network constraints in the table
  support<-tapply(tadaAICtable@printTable$AkaikeWeight, tadaAICtable@printTable$baseline,sum)
  numbers<-tapply(tadaAICtable@printTable$AkaikeWeight, tadaAICtable@printTable$baseline,length)
  return(data.frame(support=support,numberOfModels=numbers))
}


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
