#Editted for constrained model
# If type is specified as additive or multiplicative the nbdadata object is modified accordingly inside the likelihood and gradient functions
# This is not ideal for computation speed, but is just for backwards comaptibility with code written for a few analyses using the previous version 1.4


#Define class of object for the fitted additive model
setClass("tadaFit",representation(nbdaMultiDiff="character",nbdadata="list",optimisation="list",optim="list",loglik="numeric",aic="numeric",aicc="numeric",varNames="character",hessian="matrix",outputPar="numeric",se="numeric",type="character",baseline="character",noHazFunctPars="numeric",hazFunct="function",cumHaz="function"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
          signature(.Object = "tadaFit"),
          function (.Object, nbdadata,type,startValue,upper,lower,method,interval,gradient,iterations,standardErrors,baseline,noHazFunctPars,hazFunct,cumHaz,...)
          {

            if(baseline=="constant") noHazFunctPars<-1
            if(baseline=="gamma") noHazFunctPars<-2
            if(baseline=="weibull") noHazFunctPars<-2

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
              nbdadataTemp<-nbdadata[[1]]
            }else{nbdadataTemp<-nbdadata}

            if(is.na(nbdadataTemp@TADAtime1[1])){
              print("TADA times absent from nbdaData object")
              return(NULL)
            }

            #calculate the number of each type of parameter
            noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
            noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
            noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
            noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

            if(type=="asocial"){

              #We need to know whether to remove the interaction variables. This depends on whether an offset is included for any of the s parameters in any of the diffusions.
              if(is.list(nbdadata)){
                retainInt<-FALSE
                for (i in 1:length(nbdadata)){
                  nbdadataTemp2<-nbdadata[[i]];
                  if(sum(nbdadataTemp2@offsetCorrection[,1])>0) retainInt<-TRUE
                }
              }else{
                retainInt<-sum(nbdadata@offsetCorrection[,1])>0
              }

              if(retainInt){
                if(nbdadataTemp@int_ilv[1]=="ILVabsent") noILVint<-0
              }else{
                noILVint<-0
              }

              if(nbdadataTemp@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
              if(nbdadataTemp@multi_ilv[1]=="ILVabsent") noILVmulti<-0



                #Record asocialVar names
                if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
                if(nbdadataTemp@int_ilv[1]=="ILVabsent"|!retainInt){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
                if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};

                #Set staring values if not specified by the user
                if(is.null(startValue)) startValue<-c(mean(nbdadataTemp@TADAtime2),rep(1,noHazFunctPars-1),rep(0,noILVasoc+noILVint+noILVmulti));

                #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
                lower<-c(rep(0,noHazFunctPars),rep(-Inf,length(startValue)-noHazFunctPars));
                upper<-rep(Inf,length(startValue));
                #if(is.null(interval)) interval<-c(0,999);

                #Optimise for s
                fit1<-NULL
                if(gradient){
                  try(fit1<-nlminb(start=startValue, objective= asocialTadaLikelihood,gradient= asocialTadaGradient_fn, lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(iter.max=iterations)));
                }else{
                  try(fit1<-nlminb(start=startValue, objective= asocialTadaLikelihood,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(iter.max=iterations)));
                }
                if(is.null(fit1)){
                  print("Error in likeihood optimization");
                  return(NULL)
                }

                if(method=="both"){
                  if (is.null(fit1)){
                    try(fit2<-optim(par=fit1@par,fn=asocialTadaLikelihood,method="L-BFGS-B",gr=asocialTadaGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(maxit=iterations)))
                  }else{
                    try(fit2<-optim(par=startValue,fn=asocialTadaLikelihood,method="L-BFGS-B",gr=asocialTadaGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(maxit=iterations)))
                  }
                }else{fit2<-as.list(NA)}

                #Record MLEs
                outputPar<-fit1$par;

                loglik<-fit1$objective;

                #Get sample size across all diffusions
                sampleSize<-sampSizeExtract(nbdadata);


                #Calculate aic and for model without social transmission
                aic<-2*length(fit1$par)+2*loglik;
                aicc<-2*(length(fit1$par))*(sampleSize/(sampleSize-(length(fit1$par))-1))+2*loglik;

                #To prevent a low AICc when there are more parameters than data!
                if(is.nan(aic)|is.nan(aicc)){}else{
                  if(aicc<aic) aicc<-Inf;
                }


                #		}
                #Extract names of variables
                varNames<-NULL
                if(baseline=="constant") varNames<-"Scale (1/rate):"
                if(baseline=="weibull") varNames<-c("Scale (1/rate):", "Shape")
                if(baseline=="gamma") varNames<-c("Scale (1/rate):", "Shape")
                if(baseline=="custom") varNames<-paste("Baseline Parameter",1:noHazFunctPars)
                parCounter<-0
                if(!is.null(asocialVarNames)){
                  varNames<-c(varNames,paste(parCounter+(1:length(asocialVarNames)),"Asocial:",asocialVarNames))
                  parCounter<-parCounter+length(asocialVarNames)
                }
                if(!is.null(intVarNames)){
                  varNames<-c(varNames,paste(parCounter+(1:length(intVarNames)),"Social:",intVarNames))
                  parCounter<-parCounter+length(intVarNames)
                }
                if(!is.null(multiVarNames)){
                  varNames<-c(varNames,paste(parCounter+(1:length(multiVarNames)),"Social= asocial:",multiVarNames))
                }

                #Get hessian matrix and use it to get standard errors

                 if (standardErrors!=F){
                   hessianMat<-hessian(func=asocialTadaLikelihood,x=fit1$par,nbdadata=nbdadata,retainInt=retainInt,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz)
                 } else {hessianMat<-NULL}

                if(is.null(hessianMat)){
                  se<-rep(NaN,length(outputPar))
                  hessianMat<-matrix(NA)
                }else{
                  seTemp<-NULL
                  if(det(hessianMat)==0|is.na(det(hessianMat))|is.na(det(hessianMat))){
                    #If the determinant of the hessian matrix is zero it cannot be inverted
                    se<-rep(NaN,length(outputPar));
                  }else{
                    #Initialise varTemp to NaN so it is found if the hessianMat cannot be inverted
                    varTemp<-diag(hessianMat)
                    varTemp[]<-NaN
                    try(varTemp<-diag(solve(hessianMat)),silent=T);
                    varTemp[varTemp<0]<-NaN
                    seTemp<-sqrt(varTemp)
                    #Record SEs giving a zero for constrained values
                    if(is.null(seTemp)){
                      se<-rep(NaN,length(outputPar));
                    }else{
                      se<-seTemp;
                    }
                  }
                  if(!is.matrix(hessianMat))hessianMat<-matrix(hessianMat)
                }


                if(is.list(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,...)

                }
            }
            else{


              if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
              if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
              if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

              #Record asocialVar names
              if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
              if(nbdadataTemp@int_ilv[1]=="ILVabsent"){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
              if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};


              #Set staring values if not specified by the user
              if(is.null(startValue)) startValue<-c(mean(nbdadataTemp@TADAtime2),rep(1,noHazFunctPars-1),rep(0,noSParam+noILVasoc+noILVint+noILVmulti));

              #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
              if(is.null(lower)){
                lower<-rep(-Inf,length(startValue));
                lower[1:(noHazFunctPars+noSParam)]<-0;

              }
              if(is.null(upper)){
                  upper<-rep(Inf,length(startValue));
              }


              #Optimise for s
              fit1<-NULL
              if(gradient){
                try(fit1<-nlminb(start=startValue, objective= tadaLikelihood,gradient= tadaGradient_fn, lower=lower, upper=upper,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(iter.max=iterations)));
              }else{
                try(fit1<-nlminb(start=startValue, objective= tadaLikelihood,lower=lower, upper=upper,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(iter.max=iterations)));
              }
              if(is.null(fit1)){
                print("Error in likeihood optimization");
                return(NULL)
              }

              if(method=="both"){
                if (is.null(fit1)){
                  try(fit2<-optim(par=fit1$par,fn=tadaLikelihood,method="L-BFGS-B",gr=tadaGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(maxit=iterations)))
                }else{
                  try(fit2<-optim(par=startValue,fn=tadaLikelihood,method="L-BFGS-B",gr=tadaGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,control=list(maxit=iterations)))
                }
              }else{fit2<-as.list(NA)}



              #Record MLEs
              outputPar<-fit1$par;

              #Perform LRT for social transmission
              loglik<-fit1$objective;

              #Get sample size across all diffusions
              sampleSize<-sampSizeExtract(nbdadata);


              #Calculate aic and for model without social transmission
              aic<-2*length(fit1$par)+2*loglik;
              aicc<-2*(length(fit1$par))*(sampleSize/(sampleSize-(length(fit1$par))-1))+2*loglik;

              #To prevent a low AICc when there are more parameters than data!
              if(is.nan(aic)|is.nan(aicc)){}else{
                if(aicc<aic) aicc<-Inf;
              }


              #Extract names of variables
              varNames<-NULL
              if(baseline=="constant") varNames<-"Scale (1/rate):"
              if(baseline=="weibull") varNames<-c("Scale (1/rate):", "Shape")
              if(baseline=="gamma") varNames<-c("Scale (1/rate):", "Shape")
              if(baseline=="custom") varNames<-paste("Baseline Parameter",1:noHazFunctPars)
              parCounter<-0
              varNames<-c(varNames,paste(parCounter+(1:noSParam),"Social transmission",1:noSParam))
              parCounter<-parCounter+noSParam
              if(!is.null(asocialVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(asocialVarNames)),"Asocial:",asocialVarNames))
                parCounter<-parCounter+length(asocialVarNames)
              }
              if(!is.null(intVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(intVarNames)),"Social:",intVarNames))
                parCounter<-parCounter+length(intVarNames)
              }
              if(!is.null(multiVarNames)){
                varNames<-c(varNames,paste(parCounter+(1:length(multiVarNames)),"Social= asocial:",multiVarNames))
              }

              #Get hessian matrix and use it to get standard errors

              if (standardErrors!=F){
                hessianMat<-hessian(func=tadaLikelihood,x=fit1$par,nbdadata=nbdadata,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz)
              } else {hessianMat<-NULL}

              if(is.null(hessianMat)){
                se<-rep(NaN,length(outputPar))
                hessianMat<-matrix(NA)
              }else{
                seTemp<-NULL
                if(det(hessianMat)==0|is.na(det(hessianMat))){
                  #If the hessian matrix is not positive definite, SEs cannot be calculated
                  se<-rep(NaN,length(outputPar));
                }else{
                  #Initialise varTemp to NaN so it is found if the hessianMat cannot be inverted
                  varTemp<-diag(hessianMat)
                  varTemp[]<-NaN
                  try(varTemp<-diag(solve(hessianMat)),silent=T);
                  varTemp[varTemp<0]<-NaN
                  seTemp<-sqrt(varTemp)
                  #Record SEs giving a zero for constrained values
                  if(is.null(seTemp)){
                    se<-rep(NaN,length(outputPar));
                  }else{
                    se<-seTemp;
                  }
                }
                if(!is.matrix(hessianMat))hessianMat<-matrix(hessianMat)
              }

              if(is.list(nbdadata)){
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,...)
              }else{
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz,...)
              }


            }
          }
            )

#Function for implementing the initialization and choosing between normal and oada.coxme version
tadaFit<-function(nbdadata,type="social",startValue=NULL, upper=NULL,lower=NULL,interval=c(0,999), method="nlminb", gradient=T,iterations=150,standardErrors=T,baseline="constant",noHazFunctPars=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL)){
  if(type=="social"|type=="asocial"){
      return(new("tadaFit",nbdadata= nbdadata,type= type, startValue= startValue,upper=upper,lower=lower,interval= interval,method= method,gradient=gradient,iterations=iterations,standardErrors=standardErrors,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz))
  }else{
    print("Error: Invalid type of model")
    return()
  }
}



