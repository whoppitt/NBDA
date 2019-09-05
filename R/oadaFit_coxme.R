
#We need a fucntion to return the relevant coxme model for a given input:
return_coxme<-function(parVect, nbdadata,formula=NULL){

  if(is.character(nbdadata)){
    subdata <- nbdadatatemp<-eval(as.name(nbdadata[1]));
    coxmeData<-createCoxmeData(parVect,subdata)
    if (length(nbdadata)>1){
      for(i in 2:length(nbdadata)){
        subdata <- eval(as.name(nbdadata[i]));
        coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata))
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
    coxmeData<-createCoxmeData(parVect,nbdadata);
    nbdadatatemp<-nbdadata
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
  return(model)
}

return_asocial_coxme<-function(parVect, nbdadata, retainInt=TRUE,formula=NULL){

if(is.null(retainInt)){
  if(is.character(nbdadata)){
    retainInt<-FALSE
    for (i in 1:length(nbdadata)){
      nbdadatatemp2<-eval(as.name(nbdadata[i]));
      if(sum(nbdadatatemp2@offsetCorrection[,1])>0) retainInt<-TRUE
    }
  }else{
    retainInt<-sum(nbdadata@offsetCorrection[,1])>0
  }
}

if(is.character(nbdadata)){
  subdata <- nbdadatatemp<-eval(as.name(nbdadata[1]));
  noSParam<-dim(nbdadatatemp@stMetric)[2] #number of s parameters
  #Append 0s to parVect for the s parameters
  parVect<-c(rep(0,noSParam),parVect)
  coxmeData<-createCoxmeData(parVect,subdata,retainInt=retainInt)
  if (length(nbdadata)>1){
    for(i in 2:length(nbdadata)){
      subdata <- eval(as.name(nbdadata[i]));
      coxmeData<-rbind(coxmeData,createCoxmeData(parVect,subdata,retainInt=retainInt))
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
#  if(nbdadatatemp@random_effects[1]=="REabsent"){
#    formula<-formula
#  }else{
    formula<-paste(formula,paste("+(1|",nbdadatatemp@random_effects,")", collapse=""))
#  }
  formula<-as.formula(formula)
}
model<-coxme(formula=formula,data=coxmeData)
return(model)
}



#Define class of object for the fitted additive model
setClass("oadaFit_coxme",representation(nbdaMultiDiff="character",nbdadata="nbdaData",optimisation="list",loglik="numeric",aic="numeric",aicc="numeric",varNames="character",hessian="matrix",outputPar="numeric",se="numeric",type="character",REvar="list",fixedFormula="formula",randomFormula="list"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
          signature(.Object = "oadaFit_coxme"),
          function (.Object, nbdadata,type,startValue,lower,method,interval,gradient,iterations,standardErrors,formula,...)
          {

            #FOR NOW FIX standardErrors=F until I have them figured out for sure
            standardErrors<-F

            #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
            if(is.list(nbdadata)){
              nbdadatatemp<-nbdadata[[1]]
            }else{nbdadatatemp<-nbdadata}

            #calculate the number of each type of parameter
            noSParam <- dim(nbdadatatemp@stMetric)[2] #s parameters
            noILVasoc<- dim(nbdadatatemp@asocILVdata)[2] #ILV effects on asocial learning
            noILVint<- dim(nbdadatatemp@intILVdata)[2] #ILV effects on interaction (social learning)
            noILVmulti<- dim(nbdadatatemp@multiILVdata)[2] #ILV multiplicative model effects

            if(is.null(startValue)){
              #Fit a model without random effects to get the starting values for the coxme fit
              modelNoRE<-oadaFit(nbdadata= nbdadata,type=type,iterations=iterations,standardErrors=F,coxmeFit=F)
              if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){
                startValue<-modelNoRE@outputPar
              }else{
                startValue<-modelNoRE@outputPar[-((length(modelNoRE@outputPar)-noILVmulti+1):length(modelNoRE@outputPar))]
              }
            }

            if(type=="asocial"){

              #We need to know whether to remove the interaction variables. This depends on whether an offset is included for any of the s parameters in any of the diffusions.
              if(is.character(nbdadata)){
                retainInt<-FALSE
                for (i in 1:length(nbdadata)){
                  nbdadatatemp2<-eval(as.name(nbdadata[i]));
                  if(sum(nbdadatatemp2@offsetCorrection[,1])>0) retainInt<-TRUE
                }
              }else{
                retainInt<-sum(nbdadata@offsetCorrection[,1])>0
              }

              if(retainInt){
                if(nbdadatatemp@int_ilv[1]=="ILVabsent") noILVint<-0
              }else{
                noILVint<-0
              }

              if(nbdadatatemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
              if(nbdadatatemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0

              if((noILVasoc+noILVint)==0){
              #This now excludes the multi_ILVs since they are optimized in the coxme function

              noILVs<-0


              model<-return_coxme(parVect=rep(0,noSParam), nbdadata)
              loglik<--model$loglik[2]

              #Record MLEs
              outputPar<-model$coefficients;

              #Get sample size across all diffusions
              sampleSize<-sampSizeExtract(nbdadata);

              #Calculate aic and for model without social transmission
              aic<-2*noILVmulti+2*loglik;
              aicc<-2*noILVmulti*(sampleSize/(sampleSize-noILVmulti-1))+2*loglik;

              #To prevent a low AICc when there are more parameters than data!
              if(is.nan(aic)|is.nan(aicc)){}else{
Inf
              }

              #Extract names of variables
              if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){varNames<-"No variables"}else{varNames <- names(fixef(model))};

              #get standard errors from fitted coxme model

              if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){se<-NaN}else{
                ##extract vcov matrix
                vcov.mat <- as.matrix(vcov(model))
                se <- sqrt(diag(vcov.mat))
                names(se) <- varNames
              }
                hessianMat<-matrix(NA)

              if(is.null(outputPar)) outputPar<-NaN

              if(is.character(nbdadata)){
                callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadatatemp, optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)
              }else{
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)

              }

              }else{

                #Record asocialVar names
                if(nbdadatatemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadatatemp@asoc_ilv};
                if(nbdadatatemp@int_ilv[1]=="ILVabsent"|!retainInt){intVarNames<-NULL}else{intVarNames<-nbdadatatemp@int_ilv};
                if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadatatemp@multi_ilv};

                #Set staring values if not specified by the user
                if(is.null(startValue)) startValue<-rep(0,noILVasoc+noILVint);

                #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
                lower<-rep(-Inf,length(startValue));
                upper<-rep(Inf,length(startValue));
                #if(is.null(interval)) interval<-c(0,999);


                  fit1<-nlminb(start=startValue, objective= asocialLikelihood_coxme,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,control=list(iter.max=iterations));

                #Record MLEs
                outputPar<-fit1$par;

                #Update to include the parameters fitted in the coxme model
                model<-return_asocial_coxme(outputPar,nbdadata,retainInt=retainInt)
                outputPar<-c(outputPar,model$coefficients)

                #Get logLik
                loglik<-fit1$objective;

                #Get sample size across all diffusions
                sampleSize<-sampSizeExtract(nbdadata);


                #Calculate aic and for model without social transmission
                aic<-2*length(outputPar)+2*loglik;
                aicc<-2*(length(outputPar))*(sampleSize/(sampleSize-(length(outputPar))-1))+2*loglik;

                #To prevent a low AICc when there are more parameters than data!
                if(is.nan(aic)|is.nan(aicc)){}else{
  Inf
                }


                #		}
                #Extract names of variables
                varNames<-NULL
                if(!is.null(asocialVarNames))varNames<-c(varNames,paste("Asocial:",asocialVarNames))
                if(!is.null(intVarNames))varNames<-c(varNames,paste("Social:",intVarNames))
                if(!is.null(multiVarNames))varNames<-c(varNames,paste("Social= asocial:",multiVarNames))

                #Get hessian matrix and use it to get standard errors

                if(standardErrors){
                  if(multi_ilv[1]=="ILVabsent"){se<-NaN}else{
                    ##extract vcov matrix
                    vcov.mat <- as.matrix(vcov(model))
                    se <- c(rep(NA,length(fit1$par)),(sqrt(diag(vcov.mat))))
                    names(se) <- varNames
                    #For now just filling in the multi ones from the coxme model- figure out later if I can get the others
                  }
                  hessianMat<-matrix(NA)
                 }else{hessianMat<-NULL}

                if(is.null(hessianMat)){
                  se<-rep(NaN,length(outputPar))
                  hessianMat<-matrix(NA)
                }


                if(is.character(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadatatemp, optimisation=fit1,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)

                }
              }
            }
            else{


              if(nbdadatatemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
              if(nbdadatatemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
              if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

              #Record asocialVar names
              if(nbdadatatemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadatatemp@asoc_ilv};
              if(nbdadatatemp@int_ilv[1]=="ILVabsent"){intVarNames<-NULL}else{intVarNames<-nbdadatatemp@int_ilv};
              if(nbdadatatemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadatatemp@multi_ilv};

              #Set staring values if not specified by the user
              if(is.null(startValue)) startValue<-rep(0,noSParam+noILVasoc+noILVint);

              #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
              if(is.null(lower)){
                lower<-rep(-Inf,length(startValue));
                lower[1:noSParam]<-0;
              }
              upper<-rep(Inf,length(startValue));
              if(is.null(interval)) interval<-c(0,999);

              #Optimise for s
              #All being done without gradient at the moment, cannot get gradient to work
                 fit1<-nlminb(start=startValue, objective= oadalikelihood_coxme,lower=lower, upper=upper,nbdadata=nbdadata,control=list(iter.max=iterations));


                #Record MLEs
                outputPar<-fit1$par;

                #Update to include the parameters fitted in the coxme model
                model<-return_coxme(outputPar,nbdadata)
                outputPar<-c(outputPar,model$coefficients)

                #Get logLik
                loglik<-fit1$objective;

                #Get sample size across all diffusions
                sampleSize<-sampSizeExtract(nbdadata);


                #Calculate aic and for model without social transmission
                aic<-2*length(outputPar)+2*loglik;
                aicc<-2*(length(outputPar))*(sampleSize/(sampleSize-(length(outputPar))-1))+2*loglik;

                #To prevent a low AICc when there are more parameters than data!
                if(is.nan(aic)|is.nan(aicc)){}else{
                aicc<-Inf
                }


                #		}
                #Extract names of variables

                varNames<-paste("Social transmission",1:noSParam)
                if(!is.null(asocialVarNames))varNames<-c(varNames,paste("Asocial:",asocialVarNames))
                if(!is.null(intVarNames))varNames<-c(varNames,paste("Social:",intVarNames))
                if(!is.null(multiVarNames))varNames<-c(varNames,paste("Social= asocial:",multiVarNames))

                #Get hessian matrix and use it to get standard errors

                if(standardErrors){
                  if(multi_ilv[1]=="ILVabsent"){se<-NaN}else{
                    ##extract vcov matrix
                    vcov.mat <- as.matrix(vcov(model))
                    se <- c(rep(NA,length(fit1$par)),(sqrt(diag(vcov.mat))))
                    names(se) <- varNames
                    #For now just filling in the multi ones from the coxme model- figure out later if I can get the others
                  }
                  hessianMat<-matrix(NA)
                }else{hessianMat<-NULL}

                if(is.null(hessianMat)){
                  se<-rep(NaN,length(outputPar))
                  hessianMat<-matrix(NA)
                }


                if(is.character(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadatatemp, optimisation=fit1,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,REvar=model$vcoef,fixedFormula=model$formulaList$fixed,randomFormula=model$formulaList$random,...)

                }
            }
          }
            )


