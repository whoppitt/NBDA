#Editted for constrained model
# If type is specified as additive or multiplicative the nbdadata object is modified accordingly inside the likelihood and gradient functions
# This is not ideal for computation speed, but is just for backwards comaptibility with code written for a few analyses using the previous version 1.4


#Define class of object for the fitted additive model
setClass("oadaFit",representation(nbdaMultiDiff="character",nbdadata="nbdaData",optimisation="list",optim="list",loglik="numeric",aic="numeric",aicc="numeric",varNames="character",hessian="matrix",outputPar="numeric",se="numeric",type="character",SLdom="logical"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
          signature(.Object = "oadaFit"),
          function (.Object, nbdadata,type,startValue,lower,method,interval,gradient,iterations,standardErrors,SLdom,...)
          {

            #If there are multiple diffusions "borrow" the first diffusion to extract necessary parameters
            if(is.character(nbdadata)){
              nbdadataTemp<-eval(as.name(nbdadata[1]));
            }else{nbdadataTemp<-nbdadata}

            #calculate the number of each type of parameter
            noSParam <- dim(nbdadataTemp@stMetric)[2] #s parameters
            noILVasoc<- dim(nbdadataTemp@asocILVdata)[2] #ILV effects on asocial learning
            noILVint<- dim(nbdadataTemp@intILVdata)[2] #ILV effects on interaction (social learning)
            noILVmulti<- dim(nbdadataTemp@multiILVdata)[2] #ILV multiplicative model effects

            if(type=="asocial"){

              #We need to know whether to remove the interaction variables. This depends on whether an offset is included for any of the s parameters in any of the diffusions.
              if(is.character(nbdadata)){
                retainInt<-FALSE
                for (i in 1:length(nbdadata)){
                  nbdadataTemp2<-eval(as.name(nbdadata[i]));
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

              if((noILVasoc+noILVmulti+noILVint)==0){

              noILVs<-0
              #Record MLEs
              outputPar<-NaN;

              loglik<-asocialLikelihood(parVect=NULL, nbdadata)

              #Get sample size across all diffusions
              sampleSize<-sampSizeExtract(nbdadata);

              #Calculate aic and for model without social transmission
              aic<-2*loglik;
              aicc<-2*loglik;

              #To prevent a low AICc when there are more parameters than data!
              if(is.nan(aic)|is.nan(aicc)){}else{
                if(aicc<aic) aicc<-NaN;
              }

              #Extract names of variables
              varNames<-"No variables";

              #Get hessian matrix and use it to get standard errors
              se<-NaN
              hessianMat<-hessianNun<-matrix(NA)

              if(is.character(nbdadata)){
                callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadataTemp, optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat,se=se, type=type,SLdom=F,...)
              }else{
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)

              }

              }else{

                #Record asocialVar names
                if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
                if(nbdadataTemp@int_ilv[1]=="ILVabsent"|!retainInt){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
                if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};

                #Set staring values if not specified by the user
                if(is.null(startValue)) startValue<-rep(0,noILVasoc+noILVint+noILVmulti);

                #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
                lower<-rep(-Inf,length(startValue));
                upper<-rep(Inf,length(startValue));
                #if(is.null(interval)) interval<-c(0,999);

                #Optimise for s
                fit1<-NULL
                if(gradient){
                  try(fit1<-nlminb(start=startValue, objective= asocialLikelihood,gradient= asocialGradient_fn, lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,control=list(iter.max=iterations)));
                }else{
                  try(fit1<-nlminb(start=startValue, objective= asocialLikelihood,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,control=list(iter.max=iterations)));
                }
                if(is.null(fit1)){
                  print("Error in likeihood optimization");
                  return(NULL)
                }

                #if(standardErrors=="Numeric") method<-"both"
                # I am trying switching to using the numDeriv package to get these
                if(method=="both"){
                  if (is.null(fit1)){
                    try(fit2<-optim(par=fit1$par,fn=asocialLikelihood,method="L-BFGS-B",gr=asocialGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,control=list(maxit=iterations)))
                  }else{
                    try(fit2<-optim(par=startValue,fn=asocialLikelihood,method="L-BFGS-B",gr=asocialGradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,retainInt=retainInt,control=list(maxit=iterations)))
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
                  if(aicc<aic) aicc<-NaN;
                }


                #		}
                #Extract names of variables
               varNames<-NULL
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

                if(standardErrors=="Analytic"){
                  #Get hessian matrix and use it to get standard errors
                  hessianMat<-hessian_fn(fit1$par,nbdadata,type="asocial", retainInt=retainInt)
                }else{
                  if(standardErrors=="Numeric"){
                    #Get hessian matrix and use it to get standard errors
                    hessianMat<-hessian(func=asocialLikelihood,x=fit1$par,nbdadata=nbdadata,retainInt=retainInt)
                  }else{ hessianMat<-NULL}
                }

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


                if(is.character(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadataTemp, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)

                }
              }
            }
            else{
              if(!SLdom){

              if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
              if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
              if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

              #Record asocialVar names
              if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
              if(nbdadataTemp@int_ilv[1]=="ILVabsent"){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
              if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};

              #Set staring values if not specified by the user
              if(is.null(startValue)) startValue<-rep(0,noSParam+noILVasoc+noILVint+noILVmulti);

              #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
              if(is.null(lower)){
              lower<-rep(-Inf,length(startValue));
              lower[1:noSParam]<-0;
              }
              upper<-rep(Inf,length(startValue));
              if(is.null(interval)) interval<-c(0,999);

              #Optimise for s
              #All being done with nlminb at the moment
              # If type is specified as additive or multiplicative the nbdadata object is modified accordingly inside the likelihood and gradient functions
              # This is not ideal for computation speed, but is just for backwards comaptibility with a few bits of code written for the previous version 1.4
              fit1<-NULL
              if(gradient){
                try(fit1<-nlminb(start=startValue, objective= oadaLikelihood,gradient= gradient_fn, lower=lower, upper=upper,nbdadata=nbdadata,control=list(iter.max=iterations)));
              }else{
                try(fit1<-nlminb(start=startValue, objective= oadaLikelihood,lower=lower, upper=upper,nbdadata=nbdadata,type=type,control=list(iter.max=iterations)));
              }
              if(is.null(fit1)){
                print("Error in likeihood optimization");
                return(NULL)
              }

              #SEs currently set to numeric if int variables are included
              if(noILVint>0) standardErrors<-"Numeric"

              if(standardErrors=="Numeric") method<-"both"
              if(method=="both"){
                if (is.null(fit1)){
                  try(fit2<-optim(par=fit1$par,fn=oadaLikelihood,method="L-BFGS-B",gr=gradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,control=list(maxit=iterations)))
                }else{
                  try(fit2<-optim(par=startValue,fn=oadaLikelihood,method="L-BFGS-B",gr=gradient_fn,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,control=list(maxit=iterations)))
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
                if(aicc<aic) aicc<-NaN;
              }


              #		}

              #Extract names of variables
              parCounter<-0
              varNames<-paste(parCounter+(1:noSParam),"Social transmission",1:noSParam)
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



              if(standardErrors=="Analytic"){
                #Get hessian matrix and use it to get standard errors
                hessianMat<-hessian_fn(fit1$par,nbdadata,type="social", retainInt=retainInt)
              }else{
                if(standardErrors=="Numeric"){
                  #Get hessian matrix and use it to get standard errors
                  hessianMat<- hessian(func=oadaLikelihood,x=fit1$par,nbdadata=nbdaDataObject)

                }else{ hessianMat<-NULL}
              }

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

              if(is.character(nbdadata)){
                callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadataTemp, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)
              }else{
                callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)

              }
              }else
                {
                # If SLdom is T we fit a model in which social transmission is always dominant. i.e. asocial learning only works when
                # everyone has connection 0. To enable this the first s parameter is set to one and the other s parameters are estimated relative to it

                if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){noILVasoc<-0} #Ignore dummy ILV
                if(nbdadataTemp@int_ilv[1]=="ILVabsent"){noILVint<-0} #Ignore dummy ILV
                if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){noILVmulti<-0} #Ignore dummy ILV

                #Record asocialVar names
                if(nbdadataTemp@asoc_ilv[1]=="ILVabsent"){asocialVarNames<-NULL}else{asocialVarNames<-nbdadataTemp@asoc_ilv};
                if(nbdadataTemp@int_ilv[1]=="ILVabsent"){intVarNames<-NULL}else{intVarNames<-nbdadataTemp@int_ilv};
                if(nbdadataTemp@multi_ilv[1]=="ILVabsent"){multiVarNames<-NULL}else{multiVarNames<-nbdadataTemp@multi_ilv};

                #Set staring values if not specified by the user
                if(is.null(startValue)) startValue<-rep(0,(noSParam-1)+noILVasoc+noILVint+noILVmulti);

                #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
                if(is.null(lower)){
                  lower<-rep(-Inf,length(startValue));
                lower[1:noSParam]<-0;
                }
                upper<-rep(Inf,length(startValue));
                if(is.null(interval)) interval<-c(0,999);


                #Optimise for s
                #All being done with nlminb at the moment

                if(oadaLikelihood_SLdom(startValue,nbdadata)==Inf){
                  print("Starting values returned log-Lik=-Inf. This likely means that an individual with zero connection to informed individuals learned when others in the population had non-zero connections. This makes a social learning dominant model impossible for your data.")
                  return(NA)
                  }

                fit1<-NULL
                if(gradient){
                  try(fit1<-nlminb(start=startValue, objective= oadaLikelihood_SLdom,gradient= gradient_SLdom, lower=lower, upper=upper,nbdadata=nbdadata,control=list(iter.max=iterations)));
                }else{
                  try(fit1<-nlminb(start=startValue, objective= oadaLikelihood_SLdom,lower=lower, upper=upper,nbdadata=nbdadata,type=type,control=list(iter.max=iterations)));
                }
                if(is.null(fit1)){
                  print("Error in likeihood optimization");
                  return(NULL)
                }

               # if(standardErrors=="Numeric") method<-"both"
                if(method=="both"){
                  if (is.null(fit1)){
                    try(fit2<-optim(par=fit1$par,fn=oadaLikelihood_SLdom,method="L-BFGS-B",gr=gradient_SLdom,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,control=list(maxit=iterations)))
                  }else{
                    try(fit2<-optim(par=startValue,fn=oadaLikelihood_SLdom,method="L-BFGS-B",gr=gradient_SLdom,hessian=T,lower=lower, upper=upper,nbdadata=nbdadata,control=list(maxit=iterations)))
                  }
                }else{fit2<-as.list(NA)}


                #Record MLEs
                outputPar<-c(1,fit1$par);

                #Perform LRT for social transmission
                loglik<-fit1$objective;

                #Get sample size across all diffusions
                sampleSize<-sampSizeExtract(nbdadata);


                #Calculate aic and for model without social transmission
                aic<-2*(length(fit1$par)+1)+2*loglik;
                aicc<-2*((length(fit1$par))+1)*(sampleSize/(sampleSize-(length(fit1$par))-1))+2*loglik;

                # Here I add an extra parameter since I assume that a SLdom model is fitted when it seems like s
                # parameters are tending to Inf in a non-SLdom model- thus one is effectively fitting the dropped parameter to Inf
                # and thus it should count as a df.

                #To prevent a low AICc when there are more parameters than data!
                if(is.nan(aic)|is.nan(aicc)){}else{
                  if(aicc<aic) aicc<-NaN;
                }


                #		}

                #Extract names of variables
                parCounter<-0
                varNames<-paste(parCounter+(1:noSParam),"Social transmission",1:noSParam)
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



                if(standardErrors=="Analytic"){standardErrors<-"Numeric"}
                  #Hessian matrix not yet fixed for SLdom models


                if(standardErrors=="Analytic"){
                  #Get hessian matrix and use it to get standard errors
                  hessianMat<-hessian_fn(fit1$par,nbdadata)
                }else{
                  if(standardErrors=="Numeric"){
                    #Get hessian matrix and use it to get standard errors
                    hessianMat<-hessian(func=oadaLikelihood_SLdom,x=fit1$par,nbdadata=nbdaDataObject)
                  }else{ hessianMat<-NULL}
                }

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

                if(is.character(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff=nbdadata, nbdadata = nbdadataTemp, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)

                }
            }
            }
          }
            )

#Function for implementing the initialization and choosing between normal and oada.coxme version
oadaFit<-function(nbdadata,type="social",startValue=NULL, lower=NULL,interval=c(0,999), method="nlminb", gradient=T,iterations=150, standardErrors="Numeric",formula=NULL,coxmeFit=NULL,SLdom=F){
  if(type=="social"|type=="asocial"){
    if(is.null(coxmeFit)){
      #If a coxme model is not specified either way- fit a coxme model if random effects are specified and a normal model if not
      if(is.character(nbdadata)){
        nbdadataTemp<-eval(as.name(nbdadata[1]));
      }else{nbdadataTemp<-nbdadata}
      if(nbdadataTemp@random_effects[1]=="REabsent"){coxmeFit=F}else{coxmeFit=T}
    }
    if(coxmeFit){
      if(SLdom){
        print("Error: SLdom model cannot be fitted using coxme");
        return
      }else{
        return(new("oadaFit_coxme",nbdadata= nbdadata,type= type, startValue= startValue,lower=lower,interval= interval,method= method,gradient=gradient,iterations=iterations, standardErrors=standardErrors,formula=formula))
      }
    }else{
      return(new("oadaFit",nbdadata= nbdadata,type= type, startValue= startValue,lower=lower,interval= interval,method= method,gradient=gradient,iterations=iterations, standardErrors=standardErrors,SLdom=SLdom))
    }
  }else{
    print("Error: Invalid type of model")
    return()
  }
}



