#These functions create a data object specifically for the discrete TADA (dTADA).
#The form of the data is much the same except each line represents an individual that was naive for a single time period
#Time varying inputs to the function, e.g. networks, ILVs, presenceMatrix, also follow time structuring by time period rather than by event
#Random effects are not implemented but retained in the data structure for furture developments (e.g. Bayesian approach)

setClass("dTADAData", representation(label="vector", idname="vector", assMatrix="array", asoc_ilv="vector", int_ilv="vector",multi_ilv="vector",random_effects="vector",
                                    orderAcq="vector", timeAcq="vector",endTime="numeric",dTimePeriods="vector",
                                    updateTimes="vector",  demons="vector",
                                    weights="vector",statusMatrix="matrix", availabilityMatrix="matrix",presenceMatrix="matrix", event.id="vector", id="vector",
                                    time1="vector", time2="vector",TADAtime1="vector", TADAtime2="vector", status="vector", presentInDiffusion ="vector",
                                    assMatrixIndex="vector",asocialTreatment="character", stMetric="matrix", asocILVdata="matrix",intILVdata="matrix",
                                    multiILVdata="matrix",offsetCorrection="matrix",randomEffectdata="matrix"));

setMethod("initialize",
          signature(.Object = "dTADAData"),
          function (.Object, label, idname=NULL, assMatrix, asoc_ilv="ILVabsent", int_ilv="ILVabsent", multi_ilv="ILVabsent",random_effects="REabsent",orderAcq, timeAcq, endTime, dTimePeriods,
                    id=NA, event.id=NA, demons=NULL, updateTimes=NULL, presenceMatrix=NULL,
                    assMatrixIndex= rep(1,length(orderAcq)),weights=rep(1, dim(assMatrix)[1]), asocialTreatment="constant",offsetCorrection=NULL, ...)
          {

              #put the time varying association matrix into an object called assMatrixTV
              if(length(dim(assMatrix))==3){ assMatrixTV<- array(data=assMatrix,dim=c(dim(assMatrix),1))}else{assMatrixTV<-assMatrix}

              # create default numeric id vector for individuals if no names are provided, if it is, convert to a factor
              if(is.null(idname)) idname <- (1:dim(assMatrix)[1]);


              # each asoc_ vector should be a vector of character strings of matrices whose names correspond to the names of the individual-level variables (ILVs). the rows of each matrix should correspond to individuals and the columns to times at which the value of the ILV changes
              nAcq <- ifelse(any(is.na(orderAcq)), 0, length(orderAcq)) # number of acquisition events EXCLUDING demonstrators.

              #If no asoc variable is provided, set a single column of zeroes
              if(asoc_ilv[1]=="ILVabsent"|int_ilv[1]=="ILVabsent"|multi_ilv[1]=="ILVabsent"){
                ILVabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }
              if(random_effects[1]=="REabsent"){
                REabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }

              status <- presentInDiffusion <-vector(); # set up the status vector and presentInDiffusion vector

              # If there is just one asocial variable matrix for all events and times, then you will have a column matrix for each ILV, the length of the number of individuals
              # If there are more than one asocial variable matrices (i.e. time-varying covariates), then you will have a matrix for each ILV, with rows equal to the number of individuals, and columns equal to the number of acquisition events (because in OADA we are constraining this to be the case: ILVs can only change at the same time as acquisition events occur otherwise you can't obtain a marginal likelihood, Will says, only a partial likelihood)
              asoc_ilv.dim <- dim(eval(as.name(asoc_ilv[1])))[1] # specify the dimensions of assoc.array
              int_ilv.dim <- dim(eval(as.name(int_ilv[1])))[1] # specify the dimensions of assoc.array
              multi_ilv.dim <- dim(eval(as.name(multi_ilv[1])))[1] # specify the dimensions of assoc.array
              random_effects.dim <- dim(eval(as.name(random_effects[1])))[1] # specify the dimensions of assoc.array


              # create asoc.array to hold the asocial variables. depending on the treatment required: "timevarying" or "constant", create a one-matrix array or a multi-matrix array
              if (asocialTreatment=="constant"){
                asoc_ilv.array <- array(dim=c(asoc_ilv.dim, 1, length(asoc_ilv)))
                dimnames(asoc_ilv.array) <- list(NULL, NULL, asoc_ilv)
                int_ilv.array <- array(dim=c(int_ilv.dim, 1, length(int_ilv)))
                dimnames(int_ilv.array) <- list(NULL, NULL, int_ilv)
                multi_ilv.array <- array(dim=c(multi_ilv.dim, 1, length(multi_ilv)))
                dimnames(multi_ilv.array) <- list(NULL, NULL, multi_ilv)
                random_effects.array <- array(dim=c(random_effects.dim, 1, length(random_effects)))
                dimnames(random_effects.array) <- list(NULL, NULL, random_effects)

              } else {
                if (asocialTreatment=="timevarying"){
                  asoc_ilv.array <- array(dim=c(asoc_ilv.dim,endTime,length(asoc_ilv)))
                  dimnames(asoc_ilv.array) <- list(NULL, c(paste("time",c(1:endTime),sep="")), asoc_ilv)
                  int_ilv.array <- array(dim=c(int_ilv.dim,endTime,length(int_ilv)))
                  dimnames(int_ilv.array) <- list(NULL, c(paste("time",c(1:endTime),sep="")), int_ilv)
                  multi_ilv.array <- array(dim=c(multi_ilv.dim,endTime,length(multi_ilv)))
                  dimnames(multi_ilv.array) <- list(NULL, c(paste("time",c(1:endTime),sep="")), multi_ilv)
                  random_effects.array <- array(dim=c(random_effects.dim,endTime,length(random_effects)))
                  dimnames(random_effects.array) <- list(NULL, c(paste("time",c(1:endTime),sep="")), random_effects)
                }
              }

              # generate a matrix that contains the status of each individual at each acquisition event
              statusMatrix <- matrix(0, nrow=dim(assMatrix)[2], ncol=endTime+1)  # a matrix with as many rows as indivs and as many columns as acquisition events PLUS one for the demonstrators
              # create a list vector to hold the index of naive individuals after each acquisition event
              naive.id <- vector(mode="list", length=endTime)

              # if there are seeded demonstrators add the vector (which should have length dim(assMatrix)[1]) to the first column of the statusMatrix to show which individuals set out as skilled (status of 1)
              if(is.null(demons)){
                statusMatrix[,1] <- rep(0,dim(assMatrix)[1])
              } else {
                for(i in 1:(1+nAcq)){
                  statusMatrix[,i] <- demons
                }
              }

              availabilityMatrix <- statusMatrix # if there are ties the statusMatrix and the availabilityMatrix will differ (the latter gives who is available to be learned *from*). we create it as identical and modify it accordingly below

              #presenceMatrix gives who was present in the diffusion for each TIME PERIOD for a dTADA- set to 1s by default
              if(is.null(presenceMatrix)){
                presenceMatrix<-matrix(1,dim(statusMatrix)[1],dim(statusMatrix)[2]);
              }


              asocILVdata.naive <-intILVdata.naive<-multiILVdata.naive<-randomEffectdata.naive<-vector() # this will hold the individual level variables for the naive individuals

              # YOU WILL HAVE TO MAKE THIS WORK EVEN WHEN THERE ARE NO ASOCIAL VARIABLES... THINK ABOUT HOW YOU MIGHT DO THIS 20120822 (Theoni comment)
              # Will: I just put a dummy asoc variable in at the start with all 0s, when the model is fitted using oadaFit the constraints and offsets vector are
              # automatically modified to ignore the dummy variable (a zero is appended to the end of each, or 2 zeroes if type=unconstrained is specified)


              for(a in 1:length(asoc_ilv)){ # Loop through asocial variables - a loop
                  asoc_ilv.array[,,a] <- eval(as.name(asoc_ilv[a])) # evaluate each one in turn
                }
              for(a in 1:length(int_ilv)){ # Loop through asocial variables - a loop
                  int_ilv.array[,,a] <- eval(as.name(int_ilv[a])) # evaluate each one in turn
                }
              for(a in 1:length(multi_ilv)){ # Loop through asocial variables - a loop
                  multi_ilv.array[,,a] <- eval(as.name(multi_ilv[a])) # evaluate each one in turn
                }
              for(a in 1:length(random_effects)){ # Loop through asocial variables - a loop
                  random_effects.array[,,a] <- eval(as.name(random_effects)) # evaluate each one in turn
                }

              for (i in 1:endTime){ # Loop through time periods

                k <- ifelse(asocialTreatment=="constant", 1, i) # this makes sure you are treating ILVs correctly if they are constant and if they are time-varying
                statusMatrix[orderAcq[timeAcq==i],(i+1):(endTime+1)] <- 1 # give the individuals that acquired the trait a status of 1 and carry skilled status (1) through to all following acquisition events
                naive.id[[i]] <- which(statusMatrix[,i]==0) # index for naive individuals before the ith time Period event


              } # closes the i loop - nAcq (k is i or 1)
                            #Now correct the availabilityMatrix such that individuals who are not present for a time period cannot be learned from

              availabilityMatrix<-statusMatrix*presenceMatrix

              if(is.na(id[1])) {id <- paste(label,c(unlist(naive.id)), sep="_")} # id of naive individuals before each time period

              naive <- dim(assMatrix)[1]-apply(statusMatrix, 2, sum) # number of naive individuals remaining after each time period (last event will be the end of the diffusion for TADA data with incomplete diffusion)



              # work out the number of association matrices provided and set up stMetric matrix accordingly
              stMetric <- matrix(data=0, nrow=length(id), ncol=dim(assMatrix)[3])
              dimnames(stMetric) <- list(NULL, paste("stMetric",c(1:dim(assMatrix)[3]),sep=""))

              #############################################################
              # Loop through time periods - learner loop
              # learnMetric is the sum of network connections of individuals that learned at each acquisition events, to other informed individuals.
              # This will always be 0 for the first animal that learned unless there were demonstrators

              # time1 and time2 index the time period or "event period" corresponding to acquisitions
              #if(nAcq!=0){ I cut this since it should work without now I have changed nAcq to nAcq+1 for TADA
              time1 <-TADAtime1<- vector(); # time period index: time start
              time2 <-TADAtime2<- vector(); # time period index: time end
              event.id.temp <- vector(); # vector to indicate the event number (this is essentially indexing the naive individuals before each event)
              timeRunning<-0
              for(event in 1:endTime){
                time1 <- c(time1, rep(event-1, naive[event]))
                time2 <- c(time2, rep(event, naive[event]))
                TADAtime1 <- c(TADAtime1, rep(timeRunning, naive[event]))
                timeRunning<-timeRunning+dTimePeriods[event]
                TADAtime2 <- c(TADAtime2, rep(timeRunning, naive[event]))
                event.id.temp <- c(event.id.temp, rep(event, each=length(naive.id[[event]])))
              } # closes for loop through events


              if(is.na(event.id[1])) {event.id <- paste(label, event.id.temp, sep="_")}

              for(event in 1:endTime){ # event indexes the number of the time period

                assMatrix<-array(assMatrixTV[,,, assMatrixIndex[event]],dim=dim(assMatrixTV)[1:3])

                learner <- orderAcq[timeAcq==event] # learner is individual id of the animal that learned AT an event
                nonlearners <- naive.id[[event]] # nonlearners are individual id of the animals that were naive BEFORE an event

                status <- c(status, statusMatrix[unlist(naive.id[[event]]),event+1])
                presentInDiffusion<-c(presentInDiffusion,presenceMatrix[unlist(naive.id[[event]]),event])

                temp.stMetric <- vector() # reset this before the metrics for each event are calculated

                for (nonlearner in nonlearners){

                     # stMetric is the total assoc of the individuals that did NOT learn by that acquisition event, with all other already-informed individuals
                     m1 <- matrix(data=assMatrix[nonlearner,,], nrow=dim(assMatrix)[3], byrow=T) # matrix1
                     m2 <- (weights*statusMatrix[,event])*t(m1) # matrix2
                     v1 <- apply(X=m2, MARGIN=2, FUN=sum) # vector1 of rowsums
                     temp.stMetric <- rbind(temp.stMetric, v1)

                   } # closes nonlearner loop for MATRIX stMetric

                   stMetric[time2==event,] <- temp.stMetric


                   if(asoc_ilv[1]=="ILVabsent"){
                     ilv1 <-cbind("ILVabsent"=rep(0,length(nonlearners)))
                   }else{
                     if(asocialTreatment=="constant"){
                       ilv1 <- matrix(asoc_ilv.array[nonlearners, 1,],nrow=length(nonlearners))
                     }else{
                       ilv1 <- matrix(asoc_ilv.array[nonlearners, event,],nrow=length(nonlearners))
                     }
                   }# this makes sure the right column out of the asoc.array is used

                   if(int_ilv[1]=="ILVabsent"){
                     intilv1 <-cbind("ILVabsent"=rep(0,length(nonlearners)))
                   }else{
                     if(asocialTreatment=="constant"){
                       intilv1 <- matrix(int_ilv.array[nonlearners, 1,],nrow=length(nonlearners))
                     }else{
                       intilv1 <- matrix(int_ilv.array[nonlearners, event,],nrow=length(nonlearners))
                     }
                   }# this makes sure the right column out of the asoc.array is used

                   if(multi_ilv[1]=="ILVabsent"){
                     multiilv1 <-cbind("ILVabsent"=rep(0,length(nonlearners)))
                   }else{
                     if(asocialTreatment=="constant"){
                       multiilv1 <- matrix(multi_ilv.array[nonlearners, 1,],nrow=length(nonlearners))
                     }else{
                       multiilv1 <- matrix(multi_ilv.array[nonlearners, event,],nrow=length(nonlearners))
                     }
                   }# this makes sure the right column out of the asoc.array is used

                   if(random_effects[1]=="REabsent"){
                     randomeffect1 <-cbind("REabsent"=rep(0,length(nonlearners)))
                   }else{
                     if(asocialTreatment=="constant"){
                       randomeffect1 <- matrix(random_effects.array[nonlearners, 1,],nrow=length(nonlearners))
                     }else{
                       randomeffect1 <- matrix(random_effects.array[nonlearners, event,],nrow=length(nonlearners))
                     }
                   }# this makes sure the right column out of the asoc.array is used




                   asocILVdata.naive <- rbind(asocILVdata.naive, ilv1)
                   if(asoc_ilv[1]=="ILVabsent"){
                     attr(asocILVdata.naive, "dimnames") <- list(NULL,"ILVabsent")
                   }else{
                     attr(asocILVdata.naive, "dimnames") <- list(NULL,asoc_ilv)
                   }

                   intILVdata.naive <- rbind(intILVdata.naive, intilv1)
                   if(int_ilv[1]=="ILVabsent"){
                     attr(intILVdata.naive, "dimnames") <- list(NULL,"ILVabsent")
                   }else{
                     attr(intILVdata.naive, "dimnames") <- list(NULL,int_ilv)
                   }

                   multiILVdata.naive <- rbind(multiILVdata.naive, multiilv1)
                   if(multi_ilv[1]=="ILVabsent"){
                     attr(multiILVdata.naive, "dimnames") <- list(NULL,"ILVabsent")
                   }else{
                     attr(multiILVdata.naive, "dimnames") <- list(NULL,multi_ilv)
                   }
                   randomEffectdata.naive <- rbind(randomEffectdata.naive, randomeffect1)
                   if(random_effects[1]=="REabsent"){
                     attr(randomEffectdata.naive, "dimnames") <- list(NULL,"REabsent")
                   }else{
                     attr(randomEffectdata.naive, "dimnames") <- list(NULL,random_effects)
                   }


              } # closes event loop

              #############################################################
              label <- rep(label, length.out=length(id))

              if(is.null(demons)) demons <- NA;

              if(is.null(offsetCorrection)) offsetCorrection <- cbind(rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]));
              dimnames(offsetCorrection)[2]<-list(c("SocialOffsetCorr","AsocialILVOffsetCorr","InteractionOffsetCorr","MultiplicativeILVOffsetCorr"))

              callNextMethod(.Object, label=label, idname=idname, assMatrix=assMatrix, asoc_ilv=asoc_ilv, int_ilv=int_ilv,multi_ilv=multi_ilv,random_effects=random_effects, orderAcq=orderAcq, timeAcq=timeAcq, endTime=endTime,updateTimes=NA, demons=demons, weights=weights, statusMatrix=statusMatrix, availabilityMatrix=availabilityMatrix, event.id=event.id, id=id, time1=time1, time2=time2,TADAtime1=TADAtime1, TADAtime2=TADAtime2, status=status, presentInDiffusion= presentInDiffusion, presenceMatrix = presenceMatrix ,asocialTreatment=asocialTreatment, stMetric=stMetric, asocILVdata=asocILVdata.naive, intILVdata=intILVdata.naive, multiILVdata=multiILVdata.naive,randomEffectdata=randomEffectdata.naive,offsetCorrection=offsetCorrection,assMatrixIndex=assMatrixIndex)

          } # end function

) # end setMethod "initialize"

dTADAData <- function(label, idname=NULL, assMatrix, asoc=NULL, asoc_ilv="ILVabsent",int_ilv="ILVabsent",multi_ilv="ILVabsent", orderAcq, timeAcq,dTimePeriods=rep(1,endTime),
                      endTime, updateTimes=NULL, demons=NULL, presenceMatrix =NULL,assMatrixIndex= rep(1,endTime), weights=rep(1, dim(assMatrix)[1]), asocialTreatment="constant",offsetCorrection=NULL){

  # For backwards compatibility, allow the asoc argument and assume it refers to asoc_ilv

  if(!is.null(asoc)&asoc_ilv[1]=="ILVabsent") asoc_ilv<-asoc

  new("dTADAData",label=label,idname=idname,assMatrix=assMatrix,availabilityMatrix=availabilityMatrix,asoc_ilv=asoc_ilv,int_ilv=int_ilv,multi_ilv=multi_ilv,orderAcq=orderAcq,timeAcq=timeAcq,
      dTimePeriods=dTimePeriods,endTime=endTime,demons=demons, presenceMatrix = presenceMatrix, assMatrixIndex = assMatrixIndex,weights=weights,asocialTreatment=asocialTreatment);

}



