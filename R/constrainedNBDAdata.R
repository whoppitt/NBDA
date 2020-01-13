
#'Put parameter constraints on an nbdaData object or dTADAData object
#'
#'Takes an object of class \code{\link{nbdaData}} or \code{\link{dTADAData}} and adds constraints to the parameters such that a
#'simpler model can be fitted. Offsets can also be added.
#'
#'Since the contents of the \code{\link{nbdaData}} or \code{\link{dTADAData}} object determine the form of an OADA or TADA fitted
#'using \code{\link{oadaFit}} or \code{\link{tadaFit}}, parameter constraints are made by applying them to the data object. This
#'is done by specifying a numerical vector \code{constraintsVect} and optionally, a  numerical vector
#'\code{offsetVect} that specifies the offset added to each parameter (see arguments above).
#'
#'@section Warning: \code{constrainedNBDAdata} cannot be used directly to fit a model with NO social effects, i.e. ALL s
#'parameters constrained. To do this, use \code{constrainedNBDAdata} to create an object with the required constraints among the
#'ILVs, then fit a model with  \code{type="asocial"}. To constrain all s paramter to have a specific values, add an
#'\code{offsetVect} using \code{constrainedNBDAdata} then fit the model with \code{type="asocial"}.
#'
#'
#'@seealso \code{\link{filteredNBDAdata}}, \code{\link{nbdaData}}, \code{\link{dTADAData}}.
#'
#'@param nbdadata an object of class \code{\link{nbdaData}} or \code{\link{dTADAData}}
#'@param constraintsVect a numerical vector specifying the constraints to be applied to the data object, of length matching the
#'number of parameters to be fitted in the correspondng OADA model (i.e. excluding baseline parameters for a TADA). Constraints
#'are specified for s parameters, then asoc_ilv, then int_ilv, then multi_ilv. If a paramter is assigned a value of zero, it is
#'constrained to have a value of zero (other values are then possible using  \code{offsetVect} below. Non-zero parameters are
#'then assigned a integer increasing from 1. Parameters assigned the same number are constrained to have the same value. e.g.
#'If we have an ndbaData object with 3 networks, 3 asoc_ilv, 3 int_ilv and 1 multi_ilv:
#'\code{constraintsVect=c(0,1,1,2,2,0,3,4,0,5)} constrains the first network to have no effect (s1=0), equivalent to removing it
#'from the model. The second two networks are constrained to have the same effect per unit connection (s2=s3). The first and
#'second ILVs in asoc_ilv are constrained to have the same effect, the third ILV in asoc_ilv is removed. The first and second
#'ILVs in int_ilv are unconstrained, and the third ILV is again removed. Finally the only multi_ilv remains unconstrained.
#'Note that parameters of different types cannot be constrained to have the same value, i.e. constraints must be within the s
#'parameters, asoc_ilv, int_ilv or multi_ilv categories. However, an ILV can be constrained to have the same effect on asocial and
#'social learning by creating a new object using \code{\link{nbdaData}} or \code{\link{dTADAData}}, removing it from asoc_ilv and
#'int_ilv and adding it to multi_ilv.
#'@param offsetVect an optional numerical vector specifying the offsets to be applied to the data object, of length matching
#'the number of parameters to be fitted in the correspondng OADA model (i.e. excluding baseline parameters for a TADA). An offset
#'is a coefficient for a predictor variable that is fixed to have a specific value, i.e. not fitted to the data. When combined
#'with an appropriate \code{constraintsVect}, \code{offsetVect} can be used to create a number of useful constrained models. e.g.
#'If we have an ndbaData object with 3 networks and 3 asoc_ilvs, we can use \code{constraintsVect=c(0,1,2,3,4,5)} and
#'\code{offsetVect=c(2,0,0,0,0,0)} and create a model in which s1 is constrained to have a value of s1=2. This works because
#'\code{constraintsVect} constrains s1=0 and then \code{offsetVect} adds an offset of 2 to s1.  This method is used internally
#'by \code{\link{profLikCI}} to calculate confidence intervals for a specific parameter.
#'Furthermore, we can create models
#'in which two parameters are constrained to differ by a specific amount: e.g. If we have an ndbaData object with 3 networks and
#'3 asoc_ilvs, we can use \code{constraintsVect=c(1,2,2,3,4,5)} and \code{offsetVect=c(0,0,1,0,0,0)} and create a model in
#'which s3 is constrained to have a value of s3=s2+1. This works because \code{constraintsVect} constrains s2=s3 and then
#'\code{offsetVect} adds an offset of 1 to s3. This method is used internally by \code{\link{profLikCI}} to calculate confidence
#'intervals for the difference bewteen two parameters.
#'@return An object of class \code{\link{nbdaData}} or \code{\link{dTADAData}} depending on the input data.


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

  newAssMatrixDims<-dim(nbdadata@assMatrix)

  newAssMatrixDims[3]<-length(unique(sConstraintsVect[sConstraintsVect>0]))
  newNBDAdata@assMatrix<-array(NA,dim=newAssMatrixDims)

  #Rearrange assMatrix (this is just used for trueTies)
  for(i in 1:length(unique(sConstraintsVect[sConstraintsVect>0]))){
  if(sum(sConstraintsVect==unique(sConstraintsVect[sConstraintsVect>0])[i])==1){
    newNBDAdata@assMatrix[,,i,]<-nbdadata@assMatrix[,,sConstraintsVect==unique(sConstraintsVect[sConstraintsVect>0])[i],]
  }else{
    newNBDAdata@assMatrix[,,i,]<-apply(array(nbdadata@assMatrix[,,sConstraintsVect==unique(sConstraintsVect[sConstraintsVect>0])[i],],dim=dim(nbdadata@assMatrix)),c(1,2,4),sum)
  }
  }

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

  if(length(asocILVOffsetVect)>0){
    tempAsocILVoffset<-t(asocILVOffsetVect*t(nbdadata@asocILVdata))
    tempAsocILVoffset[,apply(is.na(nbdadata@asocILVdata),2,sum)>0]<-0
    newNBDAdata@offsetCorrection[,2]<-nbdadata@offsetCorrection[,2]+apply(tempAsocILVoffset,1,sum)
  }else {
    newNBDAdata@offsetCorrection[,2]<-nbdadata@offsetCorrection[,2]
  }

  if(length(intILVOffsetVect)>0){
    tempIntILVoffset<-t(intILVOffsetVect*t(nbdadata@intILVdata))
    tempIntILVoffset[,apply(is.na(nbdadata@intILVdata),2,sum)>0]<-0
    newNBDAdata@offsetCorrection[,3]<-nbdadata@offsetCorrection[,3]+apply(tempIntILVoffset,1,sum)
  }else {
    newNBDAdata@offsetCorrection[,3]<-nbdadata@offsetCorrection[,3]
  }

  if(length(multiILVOffsetVect)>0){
    tempMultiILVoffset<-t(multiILVOffsetVect*t(nbdadata@multiILVdata))
    tempMultiILVoffset[,apply(is.na(nbdadata@multiILVdata),2,sum)>0]<-0
    newNBDAdata@offsetCorrection[,4]<-nbdadata@offsetCorrection[,4]+apply(tempMultiILVoffset,1,sum)
  }else {
    newNBDAdata@offsetCorrection[,4]<-nbdadata@offsetCorrection[,4]
  }


  return(newNBDAdata)
}


