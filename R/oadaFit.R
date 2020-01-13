
#Define class of object for the fitted additive model
setClass("oadaFit",representation(nbdaMultiDiff="character",nbdadata="list",optimisation="list",optim="list",loglik="numeric",aic="numeric",aicc="numeric",varNames="character",hessian="matrix",outputPar="numeric",se="numeric",type="character",SLdom="logical"));


#Method for initializing oadaFit object- including model fitting
setMethod("initialize",
          signature(.Object = "oadaFit"),
          function (.Object, nbdadata,type,startValue,lower,upper,method,interval,gradient,iterations,standardErrors,SLdom,...)
          {

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
              nbdadataTemp<-nbdadata[[1]]
            }else{nbdadataTemp<-nbdadata}


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
                  if(aicc<aic) aicc<-Inf;
                }

                #Extract names of variables
                varNames<-"No variables";

                #Get hessian matrix and use it to get standard errors
                se<-NaN
                hessianMat<-hessianNun<-matrix(NA)

                if(is.list(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat,se=se, type=type,SLdom=F,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=list(NULL),loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)

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
                  try(fit1<-nlminb(start=startValue, objective= asocialLikelihood,gradient= asocialGradient_fn, lower=lower,upper=upper, nbdadata=nbdadata,retainInt=retainInt,control=list(iter.max=iterations)));
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
                  if(aicc<aic) aicc<-Inf;
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


                if(is.list(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=F,...)

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
                if(is.null(upper)){
                  upper<-rep(Inf,length(startValue));
                }

                if(is.null(interval)) interval<-c(0,999);

                #Optimise for s
                #All being done with nlminb at the moment
                # If type is specified as additive or multiplicative the nbdadata object is modified accordingly inside the likelihood and gradient functions
                # This is not ideal for computation speed, but is just for backwards comaptibility with a few bits of code written for the previous version 1.4
                fit1<-NULL
                if(gradient){
                  try(fit1<-nlminb(start=startValue, objective= oadaLikelihood,gradient= gradient_fn, lower=lower, upper=upper,nbdadata=nbdadata,control=list(iter.max=iterations)));
                }else{
                  try(fit1<-nlminb(start=startValue, objective= oadaLikelihood,lower=lower, upper=upper,nbdadata=nbdadata,control=list(iter.max=iterations)));
                }
                if(is.null(fit1)){
                  print("Error in likeihood optimization");
                  return(NULL)
                }

                #SEs currently set to numeric if int variables are included
                if(noILVint>0) standardErrors<-"Numeric"

                #if(standardErrors=="Numeric") method<-"both"
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
                  if(aicc<aic) aicc<-Inf;
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
                    hessianMat<- hessian(func=oadaLikelihood,x=fit1$par,nbdadata=nbdadata)

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

                if(is.list(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata, optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)

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
                  if(aicc<aic) aicc<-Inf;
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
                    hessianMat<-hessian(func=oadaLikelihood_SLdom,x=fit1$par,nbdadata=nbdadata)
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

                if(is.list(nbdadata)){
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = nbdadata,  optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)
                }else{
                  callNextMethod(.Object, nbdaMultiDiff="NA", nbdadata = list(nbdadata), optimisation=fit1,optim=fit2,loglik=loglik, aic=aic,aicc=aicc,varNames=varNames, outputPar= outputPar, hessian = hessianMat ,se=se, type=type,SLdom=SLdom,...)

                }
              }
            }
          }
)



#'Fit an order of acquisition diffusion analysis (OADA) model
#'
#'\code{oadaFit} takes diffusion data in the form of an nbdaData object (\code{\link{nbdaData}}) or a list of nbdaData
#'objects (for multiple diffusions) and fits an order of acquisition diffusion analysis (OADA) model.
#'
#'The model is fitted using maximum likelihood methods, for TADA models use \code{\link{tadaFit}}.The ILVs and random
#'effects included in the model are determined by those present in the nbdaData object(s). All nbdaData objects
#'must contain the same social networks (assMatrix must match in the third dimension), the same individual level
#'variables (ILVs) in each of the asoc_ilv, int_ilv and multi_ilv slots and the same random effects in the random_effects
#'slot. If random effects are included, the model is fitted by calls to the \code{\link[coxme]{coxme}} function in the
#'coxme package. Random effects are assumed to operate multiplicatively, i.e. affect asocial and social learning
#'differences among individuals by the same amount. If more complex effects are required then a Bayesian approach is
#'recommended (not implemented in the NBDA package). trueTies specified in nbdaData object(s) are accounted for unless
#'random effects are included. This is done by adding the likelihood across all orders of acquisition consistent with the
#'trueTies.l This is highly computationally intensive and inadvisable for any but a few true ties.
#'
#'@seealso For TADA models use \code{\link{tadaFit}}. To obtain confidence intervals see \code{\link{profLikCI}}. For
#'further details about OADA see \url{https://www.sciencedirect.com/science/article/pii/S0022519310000081} and
#'\url{https://royalsocietypublishing.org/doi/full/10.1098/rstb.2016.0418}
#'
#'@param nbdadata an object of class nbdaData (\code{\link{nbdaData}}) to fit a model to a single diffusion or a list of
#'nbdaData objects to fit a model to multiple diffusions.
#'@param type a string specifying either "social" or "asocial" model. Usually asocial models have all s parameters
#'constrained =0 and all ILVs affecting only the rate of social learning are removed (i.e. those in the int_ilv slot(s)
#'of the nbdaData object(s)). However, if a non-zero offset is present on the social transmission component, e.g. when
#'constaining all s parameters to a specific value using \code{\link{constrainedNBDAdata}}, int_ilv variables are
#'retained. This situation occurs most commonly when the function is called internally by the \code{\link{profLikCI}}
#'function.
#'@param startValue optional numeric vector giving start values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model.
#'@param lower optional numeric vector giving lower values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be 0 for all s parameters and -Inf for coefficients of
#' ILVs.
#'@param upper optional numeric vector giving upper values for the maximum likelihood optimization. Length to match
#'the number of parameters fitted in the model. By default taken to be Inf for all parameters.
#'@param interval currently non-functioning argument: can be ignored.
#'@param method character string determining which optimization algorithm is used, defaulting to "nlminb" using the
#'\code{\link[stats]{nlminb}} function. If set to "both" the optim method \code{\link[stats]{optim}} is also used and the
#'results returned for both optimization procedures.
#'@param gradient logical indicating whether the gradient function should be used during optimization.
#'@param iterations numerical determining the maximum iterations to be used during optimization. Increasing this may solve
#'convergence issues.
#'@param standardErrors a string indicating how standard errors should be calculated. Defaults to "Numeric" which uses the
#'\code{\link[numDeriv]{hessian}} function. In some cases an analytical solution can be provided using "Analytic": this can
#'increase the speed of the fit for large models. In cases where the analytic solution is not available in the package, the
#'standard errors default back to numeric. There is no advantage in accuracy to using the analytical solution. If any other
#'string is provided, no standard errors are returned.
#'@param formula a formula can be provided to customise the model if being fitted using the coxme function. At the current
#'time users are advised to leave this =NULL and the appropriate formula, including all effect specified in the nbdaData
#'objects, will be built internally.
#'@param coxmeFit logical indicating whether the \code{\link[coxme]{coxme}} function should be used to fit the model. This
#'is set to NULL by default, meaning the \code{\link[coxme]{coxme}} function will be used if random effects are included.
#'If random effects are included and coxmeFit=F the random effects will be ignored.
#'@param SLdom logical determining whether a "social learning dominant" model should be fitted. This is useful in cases
#'where the user suspects that s parameters are all =Inf- which can occur if individuals with zero connections to informed
#'individuals only learn when all other naive individuals also have zero connections to informed individuals. This means
#'that the strength of social learning relative to asocial learning will be estimated as Inf. In such cases
#'an SLdom model enables the user to judge the size of s parameters for different networks relative to one another.
#'@return If coxmeFit=F or random effects are absent the function returns an object of class oadaFit. If coxmeFit=T or
#'random effects are included the function returns an object of class oadaFit_coxme.
#'@section oadaFit components:
#'The following components of the oadaFit object are of key importance for
#'interpreting the output: \describe{
#'   \item{@@outputPar}{The maximum likelihood estimates (MLEs) for the model parameters}
#'   \item{@@varNames}{The name of the variable corresponding to each of the parameter estimates. These are numbered so
#'   the user can easily identify parameters when obtaining confidence intervals using \code{\link{profLikCI}}. The s
#'   parameters are labelled "Social transmission N" with N giving the number of the network. ILV effects on asocial
#'   learning are preceded with "Asocial:". ILV effects on social learning are preceded with "Social:". "Multiplicative"
#'   ILV effects constrained to be equal on asocial and social learning are preceded with "Social=Asocial".}
#'   \item{@@se}{The standard error for each parameter. These can not always be derived so may be NaN. The user is advised
#'   to get confidence intervals for parameters using \code{\link{profLikCI}}.}
#'   \item{@@aic}{The AIC for the model.}
#'   \item{@@aicc}{The AICc for the model: AIC adjusted for sample size, with sample size taken to be the number of
#'   acquisition events.}
#'   \item{@@loglik}{The -log-likelihood for the model. Can be used to conduct likelihood ratio tests to test hypotheses.}
#'}The oadaData object also contains the following components:
#'\describe{
#'   \item{@@nbdadata}{The data the model is fitted to, as a list of nbdaData objects.}
#'   \item{@@optimisation}{The output of the \code{\link[stats]{nlminb}} optimization alogorithm, useful for assessing
#'   convergence of the model.}
#'   \item{@@optim}{The output of the \code{\link[stats]{optim}} optimization alogorithm, where used, useful for assessing
#'   convergence of the model.}
#'   \item{@@hessian}{The hessian matrix- giving the value of the second partial derivatives of the -log-likelihood with
#'   respect to the model parameters at the maximum likelihood estimators. Used to dervive the standard errors.}
#'   \item{@@type}{The model type: "asocial" or "social".}
#'   \item{@@SLdom}{Logical showing whether a "social learning dominant" model was fitted (see above).}
#'   }
#'@section Additional oadaFit_coxme components: In addtion to the components described for oadaFit objects, the following
#'slots are present in an oadaFit_coxme object:
#'\describe{
#'   \item{@@REvar}{The estimated variance of the random effects fitted by the \code{\link[coxme]{coxme}} function.}
#'   \item{@@fixedFormula}{The formula for the fixed effects input to the \code{\link[coxme]{coxme}} function. This will
#'   include any multiplicative ILVs fitted.}
#'   \item{@@randomFormula}{The formula for the random effects input to the \code{\link[coxme]{coxme}} function.}
#'   }
#'

#Function for implementing the initialization and choosing between normal and oada.coxme version
oadaFit<-function(nbdadata,type="social",startValue=NULL, lower=NULL,upper=NULL,interval=c(0,999), method="nlminb", gradient=T,iterations=150, standardErrors="Numeric",formula=NULL,coxmeFit=NULL,SLdom=F){

  #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
  if(is.character(nbdadata)){
    newNbdaData<-list()
    for(i in 1:length(nbdadata)){
      newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
    }
    nbdadata<-newNbdaData
  }

  #  if(class(nbdadata)!="nbdaData"){
  #    if(is.list(nbdadata)){
  #      if(class(nbdadata[[1]])!="nbdaData"){
  #        cat("OADA models can only be fitted to objects of class nbdaData.\nThe object provided is of class", class(nbdadata))
  #        return(NULL)
  #      }
  #    }else{
  #      cat("OADA models can only be fitted to objects of class nbdaData.\nThe object provided is of class", class(nbdadata))
  #      return(NULL)
  #    }
  #  }
  if(type=="social"|type=="asocial"){

    #If a multiple diffusion is specified as a character vector, convert to a list (for backwards compatibility)
    if(is.character(nbdadata)){
      newNbdaData<-list()
      for(i in 1:length(nbdadata)){
        newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
      }
      nbdadata<-newNbdaData
    }

    #Check if trueTies are present
    trueTiesPresent<-F
    if(is.list(nbdadata)){
      for(i in length(nbdadata)){
        nbdadataTemp<-nbdadata[[1]];
        if(!is.null(nbdadataTemp@trueTies[[1]])) trueTiesPresent<-T
      }

    }else{
      if(!is.null(nbdadata@trueTies[[1]])) trueTiesPresent<-T
    }
    #SLdom models cannot fit trueTies yet:
    if(trueTiesPresent&SLdom){
      print("Error: trueTies not currently supported for SLdom");
      return(NULL)
    }
    if(is.null(coxmeFit)){

      #If a coxme model is not specified either way- fit a coxme model if random effects are specified and a normal model if not
      if(is.list(nbdadata)){
        nbdadataTemp<-nbdadata[[1]];
      }else{nbdadataTemp<-nbdadata}
      if(nbdadataTemp@random_effects[1]=="REabsent"){coxmeFit=F}else{coxmeFit=T}
    }
    if(coxmeFit){
      if(trueTiesPresent){
        print("Error: trueTies not supported for coxmeFit");
        return(NULL)
      }
      if(SLdom){
        print("Error: SLdom model cannot be fitted using coxme");
        return (NULL)
      }else{
        return(new("oadaFit_coxme",nbdadata= nbdadata,type= type, startValue= startValue,lower=lower,upper=upper,interval= interval,method= method,gradient=gradient,iterations=iterations, standardErrors=standardErrors,formula=formula))
      }
    }else{
      return(new("oadaFit",nbdadata= nbdadata,type= type, startValue= startValue,lower=lower,upper=upper,interval= interval,method= method,gradient=gradient,iterations=iterations, standardErrors=standardErrors,SLdom=SLdom))
    }
  }else{
    print("Error: Invalid type of model")
    return(NULL)
  }
}

