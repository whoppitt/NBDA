
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
                newNbdaData<-c(newNbdaData,list(eval(as.name(nbdadata[i]))))
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



#'Fit a time of acquisition diffusion analysis (TADA) model
#'
#'A continuous TADA (cTADA) is fitted if an nbdaData object (\code{\link{nbdaData}}) or a list of nbdaData objects (for
#'multiple  diffusions) is provided. A discrete TADA (dTADA) is fitted if a dTADAData object (\code{\link{dTADAData}} or a
#'list of dTADAData objects (for multiple diffusions) is provided.
#'
#'The model is fitted using maximum likelihood methods, for OADA models use \code{\link{oadaFit}}.The ILVs included in the
#'model are determined by those present in the nbdaData or dTADAData object(s). All nbdaData/dTADAData objects must contain
#'the same social networks (assMatrix must match in the third dimension) and the same individual level variables (ILVs) in
#'each of the asoc_ilv, int_ilv and multi_ilv slots. Random effects are not included: if random effects are required then
#'an OADA \code{\link{oadaFit}} or a Bayesian TADA is 'recommended (the latter not implemented in the NBDA package).
#'
#'@seealso For OADA models use \code{\link{oadaFit}}. To obtain confidence intervals see \code{\link{profLikCI}}. For
#'further details about cTADA see \url{https://www.sciencedirect.com/science/article/pii/S0022519310000081} and
#'\url{https://royalsocietypublishing.org/doi/full/10.1098/rstb.2016.0418}. For further details about dTADA (the original
#'version of NBDA) see \url{https://royalsocietypublishing.org/doi/10.1098/rspb.2008.1824}. For further details about modelling
#'increasing or decreasing baseline rates see \url{https://link.springer.com/article/10.3758/LB.38.3.243}
#'
#'@param nbdadata for cTADA: an object of class nbdaData (\code{\link{nbdaData}}) to fit a model to a single diffusion or
#'a list of nbdaData objects to fit a model to multiple diffusions. For dTADA: an object of class dTADAData
#'(\code{\link{dTADAData}}) to fit a model to a single diffusion or a list of dTADAData objects to fit a model to multiple
#'diffusions.
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
#'@param standardErrors logical indicating whether standard errors should be calculated.
#'@param baseline string giving the baseline rate (hazard) function to be fitted. "constant" assumes that the baseline rate does
#'not change over time, fitting a single Scale parameter controlling the reference rate of asocial learning. "gamma" and "weibull"
#'both assume that the baseline rate of learning can increase or decrease over time, as determined by a second shape parameter.
#'Shape <1 indicates a decreasing baseline rate, and shape> 1 indicates an increasing baseline rate. "custom" allows the user to
#'provide their own baseline rate function (see below).
#'@param noHazFuncPars numercial giving the number of parameters in the baseline rate (hazard) function. Only necessary if
#'baseline="custom".
#'@param hazFunct a function returning the hazard function for the baseline rate of learning over time. This must return a
#'rate as a function of time,taking the form hazFunct(parameters,time). Only necessary if baseline="custom".
#'@param cumHaz a function giving the cumulative hazard function for the baseline rate of learning over time. This must return a
#'cumulative hazard as a function of time, taking the form cumHaz(parameters,time). Only necessary if baseline="custom".
#'
#'@return An object of class tadaFit is returned.
#'
#'@section tadaFit components:
#'The following components of the tadaFit object are of key importance for
#'interpreting the output: \describe{
#'   \item{@@outputPar}{The maximum likelihood estimates (MLEs) for the model parameters}
#'   \item{@@varNames}{The name of the variable corresponding to each of the parameter estimates. These are numbered so
#'   the user can easily identify parameters when obtaining confidence intervals using \code{\link{profLikCI}}. The s
#'   parameters are labelled "Social transmission N" with N giving the number of the network. ILV effects on asocial
#'   learning are preceded with "Asocial:". ILV effects on social learning are preceded with "Social:". "Multiplicative"
#'   ILV effects constrained to be equal on asocial and social learning are preceded with "Social=Asocial". "Scale" gives the
#'   parameter estimating the reference rate of asocial learning (scale= 1/rate). If gamma or weibull baseline functions are
#'   used a "Shape" parameter is also fitted. Shape <1 indicates a decreasing baseline rate, and shape> 1 indicates an #
#'   increasing baseline rate.}
#'   \item{@@se}{The standard error for each parameter. These can not always be derived so may be NaN. The user is advised
#'   to get confidence intervals for parameters using \code{\link{profLikCI}}.}
#'   \item{@@aic}{The AIC for the model.}
#'   \item{@@aicc}{The AICc for the model: AIC adjusted for sample size, with sample size taken to be the number of
#'   acquisition events.}
#'   \item{@@loglik}{The -log-likelihood for the model. Can be used to conduct likelihood ratio tests to test hypotheses.}
#'}The tadaData object also contains the following components:
#'\describe{
#'   \item{@@nbdadata}{The data the model is fitted to, as a list of nbdaData or dTADAData objects.}
#'   \item{@@optimisation}{The output of the \code{\link[stats]{nlminb}} optimization alogorithm, useful for assessing
#'   convergence of the model.}
#'   \item{@@optim}{The output of the \code{\link[stats]{optim}} optimization alogorithm, where used, useful for assessing
#'   convergence of the model.}
#'   \item{@@hessian}{The hessian matrix- giving the value of the second partial derivatives of the -log-likelihood with
#'   respect to the model parameters at the maximum likelihood estimators. Used to dervive the standard errors.}
#'   \item{@@type}{The model type: "asocial" or "social".}
#'   \item{@@baseline}{The baseline function used.}
#'   \item{@@noHazFunctPars}{The number of parameters fitted estimating the baseline function.}
#'   \item{@@hazFunct}{The custom hazard function used, if appropriate.}
#'   \item{@@cumHaz}{The custom cumulative hazard function used, if appropriate.}
#'
#'   }

#Function for implementing the initialization and choosing between normal and oada.coxme version
tadaFit<-function(nbdadata,type="social",startValue=NULL, upper=NULL,lower=NULL,interval=c(0,999), method="nlminb", gradient=T,iterations=150,standardErrors=T,baseline="constant",noHazFunctPars=NULL,hazFunct=function() return(NULL),cumHaz=function() return(NULL)){
  if(type=="social"|type=="asocial"){
      return(new("tadaFit",nbdadata= nbdadata,type= type, startValue= startValue,upper=upper,lower=lower,interval= interval,method= method,gradient=gradient,iterations=iterations,standardErrors=standardErrors,baseline=baseline,noHazFunctPars=noHazFunctPars,hazFunct=hazFunct,cumHaz=cumHaz))
  }else{
    print("Error: Invalid type of model")
    return()
  }
}



