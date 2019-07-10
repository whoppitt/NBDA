# Editted for constrained model
# Wrapper functions also included

# A function that takes an oaData object and a set of constraints and creates a new oaData object, such that when a model is fitted the constraints are implemented
# Currently only works for additive and multiplicative not unconstrained
# For unconstrained I need another ILVdata part to the nbdaData object which I can put the constrained data for the interaction in

constrainedNBDAdata<-function(nbdadata,constraintsVect,offsetVect=NULL){

  #Calculate the number of parameters
  #calculate the number of each type of parameter
  noSParam <- dim(nbdadata@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interaction (social learning)
  noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

  if(nbdadata@asoc_ilv[1]=="ILVabsent") noILVasoc<-0
  if(nbdadata@int_ilv[1]=="ILVabsent") noILVint<-0
  if(nbdadata@multi_ilv[1]=="ILVabsent") noILVmulti<-0

  if(length(constraintsVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
    print("Error: constraintsVect must be of length equal to the number of variables in the input nbdadata object")
    return(NULL)
  }
  if(is.null(offsetVect)) offsetVect<-constraintsVect*0
  if(length(offsetVect)!=noSParam+noILVasoc+noILVint+noILVmulti){
      print("Error: offsetVect must be of length equal to the number of variables in the input nbdadata object")
      return(NULL)
    }

  sConstraintsVect<- constraintsVect[1:noSParam]
  sOffsetVect<- offsetVect[1:noSParam]

  if(nbdadata@asoc_ilv[1]=="ILVabsent"){
    asocILVConstraintsVect<-asocILVOffsetVect<-NULL
  }else{
    asocILVConstraintsVect<-constraintsVect[(noSParam+1):(noSParam+noILVasoc)]
    if(max(asocILVConstraintsVect)>0) asocILVConstraintsVect<-(asocILVConstraintsVect-min(asocILVConstraintsVect[asocILVConstraintsVect>0])+1)*(asocILVConstraintsVect>0)
    asocILVOffsetVect<-offsetVect[(noSParam+1):(noSParam+noILVasoc)]
  }

  if(nbdadata@int_ilv[1]=="ILVabsent"){
    intILVConstraintsVect<-intILVOffsetVect<-NULL
  }else{
    intILVConstraintsVect<-constraintsVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]
    if(max(intILVConstraintsVect)>0) intILVConstraintsVect<-(intILVConstraintsVect-min(intILVConstraintsVect[intILVConstraintsVect>0])+1)*(intILVConstraintsVect>0)
    intILVOffsetVect<-offsetVect[(noSParam+noILVasoc+1):(noSParam+noILVasoc+noILVint)]
  }

  if(nbdadata@multi_ilv[1]=="ILVabsent"){
    multiILVConstraintsVect<-multiILVOffsetVect<-NULL
  }else{
    multiILVConstraintsVect<-constraintsVect[(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]
    if(max(multiILVConstraintsVect)>0) multiILVConstraintsVect<-(multiILVConstraintsVect-min(multiILVConstraintsVect[multiILVConstraintsVect>0])+1)*(multiILVConstraintsVect>0)
    multiILVOffsetVect<-offsetVect[(noSParam+noILVasoc+noILVint+1):(noSParam+noILVasoc+noILVint+noILVmulti)]

  }

  newNBDAdata<-nbdadata

  tempstMetric<-NULL
  if(length(unique(sConstraintsVect)[unique(sConstraintsVect)>0])==0){
    tempstMetric<-matrix(NA)
  }else{
    for (i in unique(sConstraintsVect)[unique(sConstraintsVect)>0]){
    tempstMetric<-cbind(tempstMetric,apply(as.matrix(nbdadata@stMetric[,sConstraintsVect==i]),1,sum))
  }
  }
  newNBDAdata@stMetric<-tempstMetric

  tempasocILVdata<-newasoc<-NULL
  if(length(unique(asocILVConstraintsVect)[unique(asocILVConstraintsVect)>0])==0){
    tempasocILVdata<-as.matrix(nbdadata@time1*0)
    newasoc<-"ILVabsent"
  }else{
    for (i in unique(asocILVConstraintsVect)[unique(asocILVConstraintsVect)>0]){
      tempasocILVdata<-cbind(tempasocILVdata,apply(as.matrix(nbdadata@asocILVdata[,asocILVConstraintsVect==i]),1,sum))
      newasoc<-c(newasoc,paste(nbdadata@asoc_ilv[asocILVConstraintsVect==i],collapse="_"))
    }
  }
  newNBDAdata@asocILVdata<-tempasocILVdata
  newNBDAdata@asoc_ilv<-newasoc

  tempintILVdata<-newint<-NULL
  if(length(unique(intILVConstraintsVect)[unique(intILVConstraintsVect)>0])==0){
    tempintILVdata<-as.matrix(nbdadata@time1*0)
    newint<-"ILVabsent"
  }else{
    for (i in unique(intILVConstraintsVect)[unique(intILVConstraintsVect)>0]){
      tempintILVdata<-cbind(tempintILVdata,apply(as.matrix(nbdadata@intILVdata[,intILVConstraintsVect==i]),1,sum))
      newint<-c(newint,paste(nbdadata@int_ilv[intILVConstraintsVect==i],collapse="_"))
    }
  }
  newNBDAdata@intILVdata<-tempintILVdata
  newNBDAdata@int_ilv<-newint

  tempmultiILVdata<-newmulti<-NULL
  if(length(unique(multiILVConstraintsVect)[unique(multiILVConstraintsVect)>0])==0){
    tempmultiILVdata<-as.matrix(nbdadata@time1*0)
    newmulti<-"ILVabsent"
  }else{
    for (i in unique(multiILVConstraintsVect)[unique(multiILVConstraintsVect)>0]){
      tempmultiILVdata<-cbind(tempmultiILVdata,apply(as.matrix(nbdadata@multiILVdata[,multiILVConstraintsVect==i]),1,sum))
      newmulti<-c(newmulti,paste(nbdadata@multi_ilv[multiILVConstraintsVect==i],collapse="_"))
    }
  }
  newNBDAdata@multiILVdata<-tempmultiILVdata
  newNBDAdata@multi_ilv<-newmulti

  colnames(newNBDAdata@asocILVdata)<-newNBDAdata@asoc_ilv
  colnames(newNBDAdata@intILVdata)<-newNBDAdata@int_ilv
  colnames(newNBDAdata@multiILVdata)<-newNBDAdata@multi_ilv

  newNBDAdata@offsetCorrection[,1]<-nbdadata@offsetCorrection[,1]+apply(t(sOffsetVect*t(nbdadata@stMetric)),1,sum)
  newNBDAdata@offsetCorrection[,2]<-nbdadata@offsetCorrection[,2]+apply(t(asocILVOffsetVect*t(nbdadata@asocILVdata)),1,sum)
  newNBDAdata@offsetCorrection[,3]<-nbdadata@offsetCorrection[,3]+apply(t(intILVOffsetVect*t(nbdadata@intILVdata)),1,sum)
  newNBDAdata@offsetCorrection[,4]<-nbdadata@offsetCorrection[,4]+apply(t(multiILVOffsetVect*t(nbdadata@multiILVdata)),1,sum)

  return(newNBDAdata)
}


# These functions are wrappers for constrainNBDAdata that make particular commonly needed constraints quicker and easier

# This one constrains all non-multi ILV parameters to zero
# i.e. creates an nbdadata object for fitting a multiplicative model to
constrainToMultiOnly<-function(nbdadata){

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadata@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
  noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

  constraintsVect<-c(1:noSParam,rep(0,noILVasoc+noILVint),(noSParam+1):(noSParam+noILVmulti))

  return(constrainedNBDAdata(nbdadata,constraintsVect=constraintsVect))
}

# This one constrains all multi and int ILV parameters to zero
# i.e. creates an nbdadata object for fitting an additive model to
constrainToAddOnly<-function(nbdadata){

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadata@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
  noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

  constraintsVect<-c(1:(noSParam+noILVasoc),rep(0,noILVint+noILVmulti))

  return(constrainedNBDAdata(nbdadata,constraintsVect=constraintsVect))
}

# This one constrains all multi ILV parameters to zero
# i.e. creates an nbdadata object for fitting an unconstrained model to
removeMulti<-function(nbdadata){

  #calculate the number of each type of parameter
  noSParam <- dim(nbdadata@stMetric)[2] #s parameters
  noILVasoc<- dim(nbdadata@asocILVdata)[2] #ILV effects on asocial learning
  noILVint<- dim(nbdadata@intILVdata)[2] #ILV effects on interation (social learning)
  noILVmulti<- dim(nbdadata@multiILVdata)[2] #ILV multiplicative model effects

  constraintsVect<-c(1:(noSParam+noILVasoc+noILVint),rep(0,noILVmulti))

  return(constrainedNBDAdata(nbdadata,constraintsVect=constraintsVect))
}

# Takes an nbdadata object and constrains all asoc and int ILV parameters to zero
# But moves the asoc_ilv to the multi_ilv slot
# Thus an object with all the ILVs entered into the asoc_ilv slot can be converted
# into a multiplicative model object

convertAddToMulti<-function(nbdadata){

  nbdadata@multi_ilv<-nbdadata@asoc_ilv
  nbdadata@multiILVdata<-nbdadata@asocILVdata

  nbdadata@asoc_ilv<-"ILVabsent"
  nbdadata@int_ilv<-"ILVabsent"
  nbdadata@asocILVdata<-as.matrix(nbdadata@asocILVdata[,1]*0)
  nbdadata@intILVdata<-as.matrix(nbdadata@asocILVdata[,1]*0)

  return(nbdadata)
}



