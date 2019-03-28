# This function can be used to filter out individuals (or more generally lines of data). This removes them as learners only
# so their influence as potential transmitters of information is not removed.
# This contrasts with removing an individual before using nbdaData() to create the data object.
# The user specifies which part of the nbdaData object is to be used as a filter as a string.
# e.g. to filter out by id, enter "id". To filter by the first variable in the asocILVmatrix, use "asocILVmatrix[,1]"
# Then provide a vector to exclude saying which cases whould be removed, e.g.
# newData<-filteredNBDAdata(nbdadata=Diffusion1, filter="id", exclude=c("Diffusion1_1","Diffusion_2"))
# Would remove individuals 1 and 2 from the data, assuming the diffusion label= Diffusion1

filteredNBDAdata<-function(nbdadata, filter, exclude){

  if(!is.character(filter)) return("For filter please enter a string giving a vector within the nbdaData object to be used as a filter")

  filterVect <- eval(parse(text=paste("nbdadata@",filter,sep="")));

  if(length(filterVect)!=length(nbdadata@id)) return(paste("Please provide a filter of length",length(nbdadata@id)))

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

return(nbdadata)
}
