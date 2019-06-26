

hessian_fn <- function(parVect, nbdadata,type="social",retainInt=NULL){
		if(type=="social"){
		  hessianMat<-hessian_fn.social(parVect,nbdadata)
		}else if (type=="asocial"){
		  hessianMat<-hessian_fn.asoc(parVect,nbdadata,retainInt=retainInt)
		}else{
			print("Invalid model type")
			hessianMat <-NULL
		}
		return(hessianMat)
}



hessian_fn.social <- function(parVect, nbdadata){


  if(is.character(nbdadata)){

    totalHessian <- matrix(rep(0, length(parVect)*length(parVect)),ncol=length(parVect));

    for(i in 1:length(nbdadata)){
      subdata <- eval(as.name(nbdadata[i]));
      totalHessian <- totalHessian + hessian_fn.social(parVect= parVect, nbdadata=subdata);
    }

    return(totalHessian);

  }else{

    #Define required function
    sumWithoutNA <- function(x) sum(na.omit(x))

    #calculate the number of each type of parameter
    noSParam <- dim(nbdadata@stMetric)[2] #s parameters
    noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
    noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
    noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

    if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
    if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
    if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

    includeInOADA<-nbdadata@event.id %in% nbdadata@event.id[nbdadata@status==1]
    #Exclude the lines of data corresponding to the final period to endtime, if necedssary
    datalength <- sum(includeInOADA)

    #Extract vector giving which naive individuals were present in the diffusion for each acqusition event
    presentInDiffusion<-nbdadata@ presentInDiffusion[includeInOADA]

    #assign different paramreter values to the right vectors
    sParam <- parVect[1:noSParam]
    asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
    intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]
    multiCoef<-parVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]

    if(nbdadata@asoc_ilv[1]=="ILVabsent") asocialCoef<-NULL
    if(nbdadata@int_ilv[1]=="ILVabsent") intCoef<-NULL
    if(nbdadata@multi_ilv[1]=="ILVabsent") multiCoef<-NULL

    # create a matrix of the coefficients to multiply by the observed data values, only if there are asocial variables
    if(nbdadata@asoc_ilv[1]=="ILVabsent"){
      asocialLP<-rep(0,datalength)
    }else{
      asocialCoef.mat <- matrix(data=rep(asocialCoef, datalength), nrow=datalength, byrow=T)
      asocial.sub <- nbdadata@asocILVdata[includeInOADA,]
      asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
    }
    asocialLP<-asocialLP+nbdadata@offsetCorrection[includeInOADA,2]

    # now do the same for the interaction variables
    if(nbdadata@int_ilv[1]=="ILVabsent"){
      socialLP<-rep(0,datalength)
    }else{
      intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
      int.sub <- nbdadata@intILVdata[includeInOADA,]
      socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
    }
    socialLP<-socialLP+nbdadata@offsetCorrection[includeInOADA,3]

    # now adjust both LPs for the variables specified to have a multiplicative effect (the same effect on asocial and social learning)
    if(nbdadata@multi_ilv[1]=="ILVabsent"){
      multiLP<-rep(0,datalength)
    }else{
      multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
      multi.sub <- nbdadata@multiILVdata[includeInOADA,]
      multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
    }
    multiLP<-multiLP+nbdadata@offsetCorrection[includeInOADA,4]
    asocialLP<-asocialLP+multiLP
    socialLP<-socialLP+multiLP

    sParam.mat <- matrix(data=rep(sParam, datalength), nrow=datalength, byrow=T) # create a matrix of sParams
    unscaled.st <- apply(sParam.mat*nbdadata@stMetric[includeInOADA,], MARGIN=1, FUN=sum)
    unscaled.st<-unscaled.st+nbdadata@offsetCorrection[includeInOADA,1]

    #The totalRate is set to zero for naive individuals not in the diffusion for a given event
    totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion

    socialRate<-(exp(socialLP)*unscaled.st)* presentInDiffusion

    # gradient for any s parameter - make sure you apply it to the relevant column of stMetric matrix

    #### S PARAMETERS

    s_hess <- matrix(0, nrow=noSParam, ncol=noSParam)
    if(nbdadata@asoc_ilv[1]!="ILVabsent"){s_ilv_asoc_hess<-matrix(0,nrow=noSParam, ncol= noILVasoc)}else{s_ilv_asoc_hess<-Ts_ilv_asoc_hess<-NULL}
    if(nbdadata@int_ilv[1]!="ILVabsent"){s_ilv_int_hess<-matrix(0,nrow=noSParam, ncol= noILVint)}else{s_ilv_int_hess<-Ts_ilv_int_hess<-NULL}
    if(nbdadata@multi_ilv[1]!="ILVabsent"){s_ilv_multi_hess<-matrix(0,nrow=noSParam, ncol= noILVmulti)}else{s_ilv_multi_hess<-Ts_ilv_multi_hess<-NULL}


    for (s1 in 1:noSParam){
      for(s2 in 1:s1){
        s_hess[s1,s2]<-s_hess[s2,s1]<-sum(tapply(exp(socialLP)*nbdadata@stMetric[,s1], INDEX=nbdadata@event.id, FUN=sum)*tapply(exp(socialLP)*nbdadata@stMetric[,s2], INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-(exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s1])*(exp(socialLP[nbdadata@status==1])*nbdadata@stMetric[nbdadata@status==1,s2])/(totalRate[nbdadata@status==1]^2))
      }
      if(nbdadata@asoc_ilv[1]!="ILVabsent"){
        for(ilv in 1: noILVasoc){
          s_ilv_asoc_hess[s1,ilv]<--sum((exp(asocialLP[nbdadata@status==1])*nbdadata@asocILVdata[nbdadata@status==1,ilv]*nbdadata@stMetric[nbdadata@status==1,s1])/(totalRate[nbdadata@status==1]^2)-tapply(exp(asocialLP)*nbdadata@asocILVdata[,ilv], INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@stMetric[,s1], INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2))
        }
        # closes s2 for loop
          Ts_ilv_asoc_hess<-t(s_ilv_asoc_hess)
      }
      if(nbdadata@int_ilv[1]!="ILVabsent"){
        for(ilv in 1: noILVint){
          s_ilv_int_hess[s1,ilv]<-sum(-
                                        (nbdadata@intILVdata[nbdadata@status==1,ilv]*nbdadata@stMetric[nbdadata@status==1,s1])/(totalRate[nbdadata@status==1]) +
                                        (socialRate[nbdadata@status==1]*nbdadata@intILVdata[nbdadata@status==1,ilv]*nbdadata@stMetric[nbdadata@status==1,s1])/(totalRate[nbdadata@status==1]^2)+
                                        (tapply(nbdadata@stMetric[,s1]*nbdadata@intILVdata[,ilv], INDEX=nbdadata@event.id, FUN=sum))/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))-
                                        (tapply(nbdadata@stMetric[,s1], INDEX=nbdadata@event.id, FUN=sum)*tapply(socialRate*nbdadata@intILVdata[,ilv], INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2))
          )
        }
        Ts_ilv_int_hess<-t(s_ilv_int_hess)
      }
      if(nbdadata@multi_ilv[1]!="ILVabsent"){
        for(ilv in 1: noILVmulti){
          s_ilv_multi_hess[s1,ilv]<-			sum(tapply(nbdadata@stMetric[,s1]*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(totalRate*nbdadata@multiILVdata[,ilv]*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-
                                             tapply(nbdadata@multiILVdata[,ilv]*nbdadata@stMetric[,s1]*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)))
        }
        Ts_ilv_multi_hess<-t(s_ilv_multi_hess)
      }

    }
    # closes s1 for loop


    #### ASOCIAL (ADDITIVE MODEL) PARAMETERS

    if(nbdadata@asoc_ilv[1]!="ILVabsent"){

      asocial_hess <-matrix(0,nrow= noILVasoc, ncol= noILVasoc)

        for (ilv1 in 1:length(nbdadata@asoc_ilv)){
        for(ilv2 in 1:ilv1){
          asocial_hess[ilv1,ilv2]<-asocial_hess[ilv2,ilv1]<-sum(tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@asocILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-
                                                                  tapply(nbdadata@asocILVdata[,ilv1]*nbdadata@asocILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)-
                                                                  (exp(asocialLP[nbdadata@status==1])^2*nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@asocILVdata[nbdadata@status==1,ilv2])/(totalRate[nbdadata@status==1]^2)+
                                                                  (exp(asocialLP[nbdadata@status==1])*nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@asocILVdata[nbdadata@status==1,ilv2])/(totalRate[nbdadata@status==1])
                                                                )
        }
      }

      if(nbdadata@int_ilv[1]!="ILVabsent"){
        asocial_int_hess <-matrix(0,nrow= noILVasoc, ncol= noILVint)
        for (ilv1 in 1:length(nbdadata@asoc_ilv)){
         for(ilv2 in 1:noILVint){
          asocial_int_hess[ilv1,ilv2]<-sum(   nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]*exp(asocialLP[nbdadata@status==1])*socialRate[nbdadata@status==1]/(totalRate[nbdadata@status==1]^2)+
                                                                             tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)
                                                                 )
        }
        }
        Tasocial_int_hess<-t(asocial_int_hess)
        }else{      asocial_int_hess <-Tasocial_int_hess<-NULL}

      if(nbdadata@multi_ilv[1]!="ILVabsent"){
        asocial_multi_hess <-matrix(0,nrow= noILVasoc, ncol= noILVmulti)
        for (ilv1 in 1:length(nbdadata@asoc_ilv)){
        for(ilv2 in 1:noILVmulti){
          asocial_multi_hess[ilv1,ilv2]<-sum(     tapply(nbdadata@asocILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))+
                                                    tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)
                                          )
        }

        }
        Tasocial_multi_hess<-t(asocial_multi_hess)
      }else{asocial_multi_hess <-Tasocial_multi_hess<-NULL}

    } else {
        asocial_hess <-asocial_multi_hess <-asocial_int_hess <-Tasocial_multi_hess <-Tasocial_int_hess <- NULL
    } # closes if(nbdadata@asoc[1]!="ILVabsent")



    #### INTERACTIVE PARAMETERS

    if(nbdadata@int_ilv[1]!="ILVabsent"){

      int_hess <-matrix(0,nrow= noILVint, ncol= noILVint)

      for (ilv1 in 1:length(nbdadata@int_ilv)){

        for(ilv2 in 1:ilv1){
          int_hess[ilv1,ilv2]<-int_hess[ilv2,ilv1]<-sum(
            (socialRate[nbdadata@status==1]^2)*nbdadata@intILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]/(totalRate[nbdadata@status==1]^2)-
              (socialRate[nbdadata@status==1])*nbdadata@intILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]/(totalRate[nbdadata@status==1])-
              tapply(nbdadata@intILVdata[,ilv1]*socialRate, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)+
              tapply(nbdadata@intILVdata[,ilv1]*nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
          )
        }
      } # closes loops through asocialVar

      if(nbdadata@multi_ilv[1]!="ILVabsent"){

        int_multi_hess <-matrix(0,nrow= noILVint, ncol= noILVmulti)

        for (ilv1 in 1:length(nbdadata@int_ilv)){

          for(ilv2 in 1:noILVmulti){
            int_multi_hess[ilv1,ilv2]<-sum(-
               tapply(nbdadata@intILVdata[,ilv1]*socialRate, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*totalRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)+
                tapply(nbdadata@intILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
            )
          }
        } # closes loops through asocialVar
        Tint_multi_hess<-t(int_multi_hess)
      }else{int_multi_hess<-Tint_multi_hess<-NULL}

    } else {int_hess <- int_multi_hess<-Tint_multi_hess<-NULL} # closes if(nbdadata@multi_ilv[1]!="ILVabsent")





    #### MULTIPLICATIVE PARAMETERS

    if(nbdadata@multi_ilv[1]!="ILVabsent"){

      multi_hess <-matrix(0,nrow= noILVmulti, ncol= noILVmulti)

      for (ilv1 in 1:length(nbdadata@multi_ilv)){
        for(ilv2 in 1:ilv1){
          multi_hess[ilv1,ilv2]<-multi_hess[ilv2,ilv1]<-sum(tapply(nbdadata@multiILVdata[,ilv1]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-tapply(nbdadata@multiILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
        }
      } # closes loops through asocialVar
    } else {multi_hess <- NULL} # closes if(nbdadata@multi_ilv[1]!="ILVabsent")

      hessian <- as.matrix(rbind(cbind(s_hess, s_ilv_asoc_hess,s_ilv_int_hess,s_ilv_multi_hess),
                       cbind(Ts_ilv_asoc_hess, asocial_hess, asocial_int_hess,asocial_multi_hess),
                       cbind(Ts_ilv_int_hess, Tasocial_int_hess,int_hess, int_multi_hess) ,
                       cbind(Ts_ilv_multi_hess, Tasocial_multi_hess,Tint_multi_hess,multi_hess)   ))


    return(-hessian)
  }
} # end function








hessian_fn.asoc <- function(parVect, nbdadata,retainInt=NULL){

  #We need to know whether to remove the interaction variables. This depends on whether an offset is included for any of the s parameters in any of the diffusions.
  #This will be passed on by the model fitting function, but if the function is called independently we need to calculate this here
  if(is.null(retainInt)){
    if(is.character(nbdadata)){
      retainInt<-FALSE
      for (i in 1:length(nbdadata)){
        nbdadataTemp2<-eval(as.name(nbdadata[i]));
        if(sum(nbdadataTemp2@offsetCorrection[,1])>0) retainInt<-TRUE
      }
    }else{
      retainInt<-sum(nbdadata@offsetCorrection[,1])>0
    }
  }


  if(is.character(nbdadata)){

    totalHessian <- matrix(rep(0, length(parVect)*length(parVect)),ncol=length(parVect));

    for(i in 1:length(nbdadata)){
      subdata <- eval(as.name(nbdadata[i]));
      totalHessian <- totalHessian + hessian_fn.asoc(parVect= parVect, nbdadata=subdata);
    }

    return(totalHessian);

  }else{



    #calculate the number of each type of parameter
    noSParam <- dim(nbdadata@stMetric)[2] #s parameters
    noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
    noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

    if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
    if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
    if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

    datalength <- length(nbdadata@id) #ILV effects on social transmission

    #Extract vector giving which naive individuals were present in the diffusion for each acqusition event
    presentInDiffusion<-nbdadata@ presentInDiffusion

    #Extend par to include 0 s parameters
    parVect<-c(rep(0,noSParam),parVect)

    #Allow for the fact that the user might provide offsets to the s parameters which might need to be accounted for
    if(retainInt){
      noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
      if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
      if(length(parVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
        cat("Error: parVect wrong length. \nNote a non-zero offset is provided for the s parameters. \nparVect must include values for the interaction effects\n")
        return(NA)
      }
    }else{
      noILVint<-0
      if(length(parVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
        cat("Error: parVect wrong length. \nNote a zero offset is provided for the s parameters. \nparVect must not include values for the interaction effects\n")
        return(NA)
      }
    }

    #assign different paramreter values to the right vectors
    sParam <- parVect[1:noSParam]
    asocialCoef <- parVect[(noSParam+1):(noSParam+ noILVasoc)]
    multiCoef<-parVect[(noSParam+noILVasoc+noILVint+1):(noSParam+ noILVasoc+noILVint+noILVmulti)]

    if(nbdadata@asoc_ilv[1]=="ILVabsent") asocialCoef<-NULL
    if(nbdadata@int_ilv[1]=="ILVabsent") intCoef<-NULL
    if(nbdadata@multi_ilv[1]=="ILVabsent") multiCoef<-NULL


    # create a matrix of the coefficients to multiply by the observed data values, only if there are asocial variables
    if(nbdadata@asoc_ilv[1]=="ILVabsent"){
      asocialLP<-rep(0,datalength)
    }else{
      asocialCoef.mat <- matrix(data=rep(asocialCoef, datalength), nrow=datalength, byrow=T)
      asocial.sub <- nbdadata@asocILVdata
      asocialLP <- apply(asocialCoef.mat*asocial.sub, MARGIN=1, FUN=sum)
    }
    asocialLP<-asocialLP+nbdadata@offsetCorrection[,2]

    # now calculate the multiplicative LP and add to the asocial LP
    if(nbdadata@multi_ilv[1]=="ILVabsent"){
      multiLP<-rep(0,datalength)
    }else{
      multiCoef.mat <- matrix(data=rep(multiCoef, datalength), nrow=datalength, byrow=T)
      multi.sub <- nbdadata@multiILVdata
      multiLP <- apply(multiCoef.mat*multi.sub, MARGIN=1, FUN=sum)
    }
    multiLP<-multiLP+nbdadata@offsetCorrection[,4]
    asocialLP<-asocialLP+multiLP

    unscaled.st<-nbdadata@offsetCorrection[,1]

    #Allow for the fact that the user might provide offsets to the s parameters which might need to be accounted for
    if(sum(unscaled.st)>0){

      #assign different paramreter values to the right vectors
      if(noILVint>0) intCoef<- parVect[(noSParam+noILVasoc+1):(noSParam+ noILVasoc+noILVint)]

      #interaction variables
      if(nbdadata@int_ilv[1]=="ILVabsent"){
        socialLP<-rep(0,datalength)
      }else{
        intCoef.mat <- matrix(data=rep(intCoef, datalength), nrow=datalength, byrow=T)
        int.sub <- nbdadata@intILVdata
        socialLP <- apply(intCoef.mat*int.sub, MARGIN=1, FUN=sum)
      }
      # calculate
      socialLP<-socialLP+nbdadata@offsetCorrection[,3]+multiLP
    }else{socialLP<-rep(0,datalength)}

    #The totalRate is set to zero for naive individuals not in the diffusion for a given event
    totalRate <- (exp(asocialLP) + exp(socialLP)*unscaled.st)* presentInDiffusion

    socialRate<-(exp(socialLP)*unscaled.st)* presentInDiffusion

    #### ASOCIAL PARAMETERS


    #### ASOCIAL (ADDITIVE MODEL) PARAMETERS

    if(nbdadata@asoc_ilv[1]!="ILVabsent"){

      asocial_hess <-matrix(0,nrow= noILVasoc, ncol= noILVasoc)

      for (ilv1 in 1:length(nbdadata@asoc_ilv)){
        for(ilv2 in 1:ilv1){
          asocial_hess[ilv1,ilv2]<-asocial_hess[ilv2,ilv1]<-sum(tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@asocILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-
                                                                  tapply(nbdadata@asocILVdata[,ilv1]*nbdadata@asocILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)-
                                                                  (exp(asocialLP[nbdadata@status==1])^2*nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@asocILVdata[nbdadata@status==1,ilv2])/(totalRate[nbdadata@status==1]^2)+
                                                                  (exp(asocialLP[nbdadata@status==1])*nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@asocILVdata[nbdadata@status==1,ilv2])/(totalRate[nbdadata@status==1])
          )
        }
      }

      if(nbdadata@int_ilv[1]!="ILVabsent"&retainInt){
        asocial_int_hess <-matrix(0,nrow= noILVasoc, ncol= noILVint)
        for (ilv1 in 1:length(nbdadata@asoc_ilv)){
          for(ilv2 in 1:noILVint){
            asocial_int_hess[ilv1,ilv2]<-sum(   nbdadata@asocILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]*exp(asocialLP[nbdadata@status==1])*socialRate[nbdadata@status==1]/(totalRate[nbdadata@status==1]^2)+
                                                  tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)
            )
          }
        }
        Tasocial_int_hess<-t(asocial_int_hess)
      }else{      asocial_int_hess <-Tasocial_int_hess<-NULL}

      if(nbdadata@multi_ilv[1]!="ILVabsent"){
        asocial_multi_hess <-matrix(0,nrow= noILVasoc, ncol= noILVmulti)
        for (ilv1 in 1:length(nbdadata@asoc_ilv)){
          for(ilv2 in 1:noILVmulti){
            asocial_multi_hess[ilv1,ilv2]<-sum(     tapply(nbdadata@asocILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))+
                                                      tapply(nbdadata@asocILVdata[,ilv1]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*exp(asocialLP)*presentInDiffusion, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)
            )
          }

        }
        Tasocial_multi_hess<-t(asocial_multi_hess)
      }else{asocial_multi_hess <-Tasocial_multi_hess<-NULL}

    } else {
      asocial_hess <-asocial_multi_hess <-asocial_int_hess <-Tasocial_multi_hess <-Tasocial_int_hess <- NULL
    } # closes if(nbdadata@asoc[1]!="ILVabsent")



    #### INTERACTIVE PARAMETERS

    if(nbdadata@int_ilv[1]!="ILVabsent"&retainInt){

      int_hess <-matrix(0,nrow= noILVint, ncol= noILVint)

      for (ilv1 in 1:length(nbdadata@int_ilv)){

        for(ilv2 in 1:ilv1){
          int_hess[ilv1,ilv2]<-int_hess[ilv2,ilv1]<-sum(
            (socialRate[nbdadata@status==1]^2)*nbdadata@intILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]/(totalRate[nbdadata@status==1]^2)-
              (socialRate[nbdadata@status==1])*nbdadata@intILVdata[nbdadata@status==1,ilv1]*nbdadata@intILVdata[nbdadata@status==1,ilv2]/(totalRate[nbdadata@status==1])-
              tapply(nbdadata@intILVdata[,ilv1]*socialRate, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)+
              tapply(nbdadata@intILVdata[,ilv1]*nbdadata@intILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
          )
        }
      } # closes loops through asocialVar

      if(nbdadata@multi_ilv[1]!="ILVabsent"){

        int_multi_hess <-matrix(0,nrow= noILVint, ncol= noILVmulti)

        for (ilv1 in 1:length(nbdadata@int_ilv)){

          for(ilv2 in 1:noILVmulti){
            int_multi_hess[ilv1,ilv2]<-sum(-
                                             tapply(nbdadata@intILVdata[,ilv1]*socialRate, INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*totalRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)+
                                             tapply(nbdadata@intILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*socialRate, INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
            )
          }
        } # closes loops through asocialVar
        Tint_multi_hess<-t(int_multi_hess)
      }else{int_multi_hess<-Tint_multi_hess<-NULL}

    } else {int_hess <- int_multi_hess<-Tint_multi_hess<-NULL} # closes if(nbdadata@multi_ilv[1]!="ILVabsent")





    #### MULTIPLICATIVE PARAMETERS

    if(nbdadata@multi_ilv[1]!="ILVabsent"){

      multi_hess <-matrix(0,nrow= noILVmulti, ncol= noILVmulti)

      for (ilv1 in 1:length(nbdadata@multi_ilv)){
        for(ilv2 in 1:ilv1){
          multi_hess[ilv1,ilv2]<-multi_hess[ilv2,ilv1]<-sum(tapply(nbdadata@multiILVdata[,ilv1]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)*tapply(nbdadata@multiILVdata[,ilv2]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)/(tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum)^2)-tapply(nbdadata@multiILVdata[,ilv1]*nbdadata@multiILVdata[,ilv2]*(totalRate), INDEX=nbdadata@event.id, FUN=sum)/tapply(totalRate, INDEX=nbdadata@event.id, FUN=sum))
        }
      } # closes loops through asocialVar
    } else {multi_hess <- NULL} # closes if(nbdadata@multi_ilv[1]!="ILVabsent")

    hessian <- as.matrix(rbind(cbind(asocial_hess, asocial_int_hess,asocial_multi_hess),
                               cbind(Tasocial_int_hess,int_hess, int_multi_hess) ,
                               cbind(Tasocial_multi_hess,Tint_multi_hess,multi_hess)   ))


    return(-hessian)

  }
} # end function


