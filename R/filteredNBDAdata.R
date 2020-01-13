#'Filter out individuals from an nbdaData object or dTADAData object.
#'
#'This function can be used to filter out individuals (or more generally lines of data). This removes them as learners only
#'so their influence as potential transmitters of information is not removed.
#'
#'The user may wish to remove individuals from the NBDA to examine their influence on the results, or because they have good
#'reason to consider them to be unrepresentative outliers. Removing these individuals from the diffusion when creating the
#'\code{\link{nbdaData}} or \code{\link{dTADAData}} object has the disadvantage that it removes those individuals as potential
#'transmitters of information thus breaking up pathways of social transmission. Using \code{filteredNBDAdata} to remove
#'individuals retains their potential influence as transmitters of information but removes the as learners when fitting the model.
#'
#'@seealso \code{\link{constrainedNBDAdata}}, \code{\link{nbdaData}}, \code{\link{dTADAData}}.
#'
#'@param nbdadata an object of class \code{\link{nbdaData}} or \code{\link{dTADAData}}
#'@param filter a string specifying the part of the nbdaData or dTADAData object to use as a filter variable. Setting
#'\code{filter="id"} is the best easiest way to filter specific individuals (see \code{exclude} below). Alternatively, the user
#'could filter out individuals based on a specific ILV, e.g. \code{filter="asocILVdata[,1]"} uses the first asoc_ilv,
#'\code{filter="intILVdata[,2]"} uses the second int_ilv and  \code{filter="multiILVdata[,1]"} uses the first multi_ilv.
#'@param exclude a vector specifying the cases to be excluded. If "id" is used as a filter the target individuals can be specified
#'by concatenating as a string the label of the diffusion, "_", and the row number of the individual in assMatrix. e.g. if the
#'diffusion label is "Diffusion1" and we wish to filter out individuals 5 and 11, we would specify
#'\code{exclude=c("Diffusion1_5","Diffusion1_11")}. If we have specified a specific ILV as a filter, we provide a vector giving the
#'values of the ILV to be excluded. e.g. if we have a sex ILV specifying males=-0.5 and females=0.5, we can use
#'\code{exclude=-0.5} to exclude males from the analysis.
#'
#'@return An object of class \code{\link{nbdaData}} or \code{\link{dTADAData}} depending on the input data.

filteredNBDAdata<-function(nbdadata, filter, exclude){

  if(!is.character(filter)) return("For filter please enter a string giving a vector within the nbdaData object to be used as a filter")

  filterVect <- eval(parse(text=paste("nbdadata@",filter,sep="")));

  if(length(filterVect)!=length(nbdadata@id)) return(paste("Please provide a filter of length",length(nbdadata@id)))

  nbdadataOld<-nbdadata

  for(i in 1:length(exclude)){
    nbdadata@label<-nbdadata@label[filterVect!=exclude[i]]
    nbdadata@event.id<-nbdadata@event.id[filterVect!=exclude[i]]
    nbdadata@id<-nbdadata@id[filterVect!=exclude[i]]
    nbdadata@time1<-nbdadata@time1[filterVect!=exclude[i]]
    nbdadata@time2<-nbdadata@time2[filterVect!=exclude[i]]
    nbdadata@status<-nbdadata@status[filterVect!=exclude[i]]
    nbdadata@presentInDiffusion<-nbdadata@presentInDiffusion[filterVect!=exclude[i]]
    nbdadata@stMetric<-as.matrix(nbdadata@stMetric[filterVect!=exclude[i],])
    nbdadata@asocILVdata<-as.matrix(nbdadata@asocILVdata[filterVect!=exclude[i],])
    nbdadata@intILVdata<-as.matrix(nbdadata@intILVdata[filterVect!=exclude[i],])
    nbdadata@multiILVdata<-as.matrix(nbdadata@multiILVdata[filterVect!=exclude[i],])
    nbdadata@randomEffectdata<-as.matrix(nbdadata@randomEffectdata[filterVect!=exclude[i],])
    nbdadata@offsetCorrection<-as.matrix(nbdadata@offsetCorrection[filterVect!=exclude[i],])
    filterVect<-filterVect[filterVect!=exclude[i]]
  }

  #Now we need to remove the non-learner data for any events that have been removed
  #Build an index saying if each line should be included
  includeIndex<-rep(NA,length(nbdadata@event.id))
  for(i in 1:length(nbdadata@event.id)){
    includeIndex[i]<-sum(nbdadata@event.id[i]==nbdadata@event.id[nbdadata@status==1])>0
  }
  #Now apply it to the data
  nbdadata@label<-nbdadata@label[includeIndex]
  nbdadata@event.id<-nbdadata@event.id[includeIndex]
  nbdadata@id<-nbdadata@id[includeIndex]
  nbdadata@time1<-nbdadata@time1[includeIndex]
  nbdadata@time2<-nbdadata@time2[includeIndex]
  nbdadata@status<-nbdadata@status[includeIndex]
  nbdadata@presentInDiffusion<-nbdadata@presentInDiffusion[includeIndex]
  nbdadata@stMetric<-as.matrix(nbdadata@stMetric[includeIndex,])
  nbdadata@asocILVdata<-as.matrix(nbdadata@asocILVdata[includeIndex,])
  nbdadata@intILVdata<-as.matrix(nbdadata@intILVdata[includeIndex,])
  nbdadata@multiILVdata<-as.matrix(nbdadata@multiILVdata[includeIndex,])
  nbdadata@randomEffectdata<-as.matrix(nbdadata@randomEffectdata[includeIndex,])
  nbdadata@offsetCorrection<-as.matrix(nbdadata@offsetCorrection[includeIndex,])

  dimnames(nbdadata@randomEffectdata)<-dimnames(nbdadataOld@randomEffectdata)
  dimnames(nbdadata@asocILVdata)<-dimnames(nbdadataOld@asocILVdata)
  dimnames(nbdadata@intILVdata)<-dimnames(nbdadataOld@intILVdata)
  dimnames(nbdadata@multiILVdata)<-dimnames(nbdadataOld@multiILVdata)



return(nbdadata)
}

