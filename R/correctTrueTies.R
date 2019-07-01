#Functions to correct OADA likelihoods for trueTies
#We assume that time varying ILVs and dynamic networks do not change within the course of a tie

correctSingleTrueTie<-function(parVect, nbdadata,tiePosition){

  #Extract status at start of tie
  startStatus<-nbdadata@statusMatrix[,tiePosition[1]]

  if(sum(nbdadata@offsetCorrection)>0){
    print("true tie correction not supported where an offsetCorrection is included")
    return(NULL)
  }
  orderAcq<-nbdadata@orderAcq
  assMatrixIndex= rep(1,length(orderAcq))

  tempTieData<-nbdaData(label="tempTieData", assMatrix=nbdadata@assMatrix, asoc_ilv=nbdadata@asoc_ilv,int_ilv=nbdadata@int_ilv,multi_ilv=nbdadata@multi_ilv,
                        random_effects=nbdadata@random_effects, orderAcq=nbdadata@orderAcq[tiePosition],
                        ties=nbdadata@ties[tiePosition],updateTimes=updateTimes, demons=startStatus, presenceMatrix =nbdadata@presenceMatrix[,tiePosition],
                        assMatrixIndex= nbdadata@assMatrixIndex, weights=nbdadata@weights, asocialTreatment=nbdadata@asocialTreatment)

  givenOrderLogLik<-oadaLikelihood(parVect,tempTieData)



  #Now find the number of possible permutations within the tie

  perms<-permn(nbdadata@orderAcq[tiePosition]);


  #Now record the likelihood for each possible order within the tie
  likelihoodrecord<-vector(length=length(perms));
  for(i in 1:length(perms)){
    tempTieData<-nbdaData(label="tempTieData", assMatrix=nbdadata@assMatrix, asoc_ilv=nbdadata@asoc_ilv,int_ilv=nbdadata@int_ilv,multi_ilv=nbdadata@multi_ilv,
                          random_effects=nbdadata@random_effects, orderAcq=perms[[i]],
                          ties=nbdadata@ties[tiePosition],updateTimes=updateTimes, demons=startStatus, presenceMatrix =nbdadata@presenceMatrix[,tiePosition],
                          assMatrixIndex= nbdadata@assMatrixIndex, weights=nbdadata@weights, asocialTreatment=nbdadata@asocialTreatment)
    likelihoodrecord[i]<-oadaLikelihood(parVect,tempTieData)
  }

  #Calculate the total likelihood of the observed data within the tie
  logLikTie<--log(sum(exp(-likelihoodrecord)))

  #Calculate required adjustment by taking away the likelihood of the order given and adding the total likelihood of any order that results in the observed tie
  adjustment<--givenOrderLogLik+logLikTie;
  return(adjustment);
}

correctTrueTies<-function(parVect, nbdadata){

  adjustmentVector<-rep(NA,length(nbdadata@trueTies))

  for(i in 1:length(nbdadata@trueTies)){
    adjustmentVector[i]<-correctSingleTrueTie(parVect,nbdadata,tiePosition=nbdadata@trueTies[[i]])
  }
  return(sum(adjustmentVector))
}

asocialCorrectSingleTrueTie<-function(parVect, nbdadata,tiePosition){

  #Extract status at start of tie
  startStatus<-nbdadata@statusMatrix[,tiePosition[1]]

  if(sum(nbdadata@offsetCorrection)>0){
    print("true tie correction not supported where an offsetCorrection is included")
    return(NULL)
  }
  orderAcq<-nbdadata@orderAcq
  assMatrixIndex= rep(1,length(orderAcq))

  tempTieData<-nbdaData(label="tempTieData", assMatrix=nbdadata@assMatrix, asoc_ilv=nbdadata@asoc_ilv,int_ilv=nbdadata@int_ilv,multi_ilv=nbdadata@multi_ilv,
                        random_effects=nbdadata@random_effects, orderAcq=nbdadata@orderAcq[tiePosition],
                        ties=nbdadata@ties[tiePosition],updateTimes=updateTimes, demons=startStatus, presenceMatrix =nbdadata@presenceMatrix[,tiePosition],
                        assMatrixIndex= nbdadata@assMatrixIndex, weights=nbdadata@weights, asocialTreatment=nbdadata@asocialTreatment)

  givenOrderLogLik<-asocialLikelihood(parVect,tempTieData)

  #Now find the number of possible permutations within the tie

  perms<-permn(nbdadata@orderAcq[tiePosition]);


  #Now record the likelihood for each possible order within the tie
  likelihoodrecord<-vector(length=length(perms));
  for(i in 1:length(perms)){
    tempTieData<-nbdaData(label="tempTieData", assMatrix=nbdadata@assMatrix, asoc_ilv=nbdadata@asoc_ilv,int_ilv=nbdadata@int_ilv,multi_ilv=nbdadata@multi_ilv,
                          random_effects=nbdadata@random_effects, orderAcq=perms[[i]],
                          ties=nbdadata@ties[tiePosition],updateTimes=updateTimes, demons=startStatus, presenceMatrix =nbdadata@presenceMatrix[,tiePosition],
                          assMatrixIndex= nbdadata@assMatrixIndex, weights=nbdadata@weights, asocialTreatment=nbdadata@asocialTreatment)
    likelihoodrecord[i]<-asocialLikelihood(parVect,tempTieData)
  }

  #Calculate the total likelihood of the observed data within the tie
  logLikTie<--log(sum(exp(-likelihoodrecord)))

  #Calculate required adjustment by taking away the likelihood of the order given and adding the total likelihood of any order that results in the observed tie
  adjustment<--givenOrderLogLik+logLikTie;
  return(adjustment);
}

asocialCorrectTrueTies<-function(parVect, nbdadata){

  adjustmentVector<-rep(NA,length(nbdadata@trueTies))

  for(i in 1:length(nbdadata@trueTies)){
    adjustmentVector[i]<-asocialCorrectSingleTrueTie(parVect,nbdadata,tiePosition=nbdadata@trueTies[[i]])
  }
  return(sum(adjustmentVector))
}

