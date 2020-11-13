setClass("nbdaData", representation(label="vector", idname="vector", assMatrix="array", asoc_ilv="vector", int_ilv="vector",multi_ilv="vector",random_effects="vector",orderAcq="vector", timeAcq="vector",endTime="numeric",updateTimes="vector", ties="vector", trueTies="list", demons="vector", weights="vector",statusMatrix="matrix", availabilityMatrix="matrix",presenceMatrix="matrix", event.id="vector", id="vector", time1="vector", time2="vector",TADAtime1="vector", TADAtime2="vector", status="vector", presentInDiffusion ="vector", assMatrixIndex="vector",asocialTreatment="character", stMetric="matrix", asocILVdata="matrix",intILVdata="matrix",multiILVdata="matrix",offsetCorrection="matrix",randomEffectdata="matrix"));

setMethod("initialize",
          signature(.Object = "nbdaData"),
          function (.Object, label, idname=NULL, assMatrix, asoc_ilv="ILVabsent", int_ilv="ILVabsent", multi_ilv="ILVabsent",random_effects="REabsent",orderAcq, timeAcq, endTime, ties=NULL, trueTies=list(NULL), id=NA, event.id=NA, demons=NULL, updateTimes=NULL, presenceMatrix=NULL,
                    assMatrixIndex= rep(1,length(orderAcq)),weights=rep(1, dim(assMatrix)[1]), asocialTreatment="constant",offsetCorrection=NULL, ...)
          {

            if(is.na(timeAcq[1])){

              #put the time varying association matrix into an object called assMatrixTV
              if(length(dim(assMatrix))==3){ assMatrixTV<- array(data=assMatrix,dim=c(dim(assMatrix),1))}else{assMatrixTV<-assMatrix}

              # create default numeric id vector for individuals if no names are provided, if it is, convert to a factor
              if(is.null(idname)) idname <- (1:dim(assMatrix)[1]);
              # if there are no ties, make a vector of zeroes
              if(is.null(ties)) ties <- rep(0,length(orderAcq));



              # each asoc_ vector should be a vector of character strings of matrices whose names correspond to the names of the individual-level variables (ILVs). the rows of each matrix should correspond to individuals and the columns to times at which the value of the ILV changes
              nAcq <- ifelse(any(is.na(orderAcq)), 0, length(orderAcq)) # number of acquisition events EXCLUDING demonstrators.
              time1 <- vector(); # time period index: time start
              time2 <- vector(); # time period index: time end
              event.id.temp <- vector(); # vector to indicate the event number (this is essentially indexing the naive individuals before each event)

              #If no asoc variable is provided, set a single column of zeroes
              if(asoc_ilv[1]=="ILVabsent"|int_ilv[1]=="ILVabsent"|multi_ilv[1]=="ILVabsent"){
                ILVabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }
              if(random_effects[1]=="REabsent"){
                REabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }

              totalMetric <- vector() # total association of the individual that DID learn at an acquisition event, with all other individuals
              learnMetric <- vector(); # total associations of the individual that DID learn at an acquisition event, with all other individuals that have already learned

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
                  asoc_ilv.array <- array(dim=c(asoc_ilv.dim,nAcq,length(asoc_ilv)))
                  dimnames(asoc_ilv.array) <- list(NULL, c(paste("time",c(1:nAcq),sep="")), asoc_ilv)
                  int_ilv.array <- array(dim=c(int_ilv.dim,nAcq,length(int_ilv)))
                  dimnames(int_ilv.array) <- list(NULL, c(paste("time",c(1:nAcq),sep="")), int_ilv)
                  multi_ilv.array <- array(dim=c(multi_ilv.dim,nAcq,length(multi_ilv)))
                  dimnames(multi_ilv.array) <- list(NULL, c(paste("time",c(1:nAcq),sep="")), multi_ilv)
                  random_effects.array <- array(dim=c(random_effects.dim,nAcq,length(random_effects)))
                  dimnames(random_effects.array) <- list(NULL, c(paste("time",c(1:nAcq),sep="")), random_effects)
                }
              }

              # generate a matrix that contains the status of each individual at each acquisition event
              statusMatrix <- matrix(0, nrow=dim(assMatrix)[2], ncol=1+nAcq)  # a matrix with as many rows as indivs and as many columns as acquisition events PLUS one for the demonstrators
              # create a list vector to hold the index of naive individuals after each acquisition event
              naive.id <-naive.id.names<- vector(mode="list", length=nAcq)

              # if there are seeded demonstrators add the vector (which should have length dim(assMatrix)[1]) to the first column of the statusMatrix to show which individuals set out as skilled (status of 1)
              if(is.null(demons)){
                statusMatrix[,1] <- rep(0,dim(assMatrix)[1])
              } else {
                for(i in 1:(1+nAcq)){
                  statusMatrix[,i] <- demons
                }
              }

              availabilityMatrix <- statusMatrix # if there are ties the statusMatrix and the availabilityMatrix will differ (the latter gives who is available to be learned *from*). we create it as identical and modify it accordingly below

              #presenceMatrix gives who was present in the diffusion for each event- set to 1s by default
              if(is.null(presenceMatrix)){
                presenceMatrix<-statusMatrix;
                presenceMatrix[,]<-1;
              }else{
                #Add a column to the start of the presenceMatrix so the dimensions match statusMatrix
                presenceMatrix<-cbind(presenceMatrix[,1], presenceMatrix)
              }


              asocILVdata.naive <-intILVdata.naive<-multiILVdata.naive<-randomEffectdata.naive<-vector() # this will hold the individual level variables for the naive individuals

              # YOU WILL HAVE TO MAKE THIS WORK EVEN WHEN THERE ARE NO ASOCIAL VARIABLES... THINK ABOUT HOW YOU MIGHT DO THIS 20120822 (Theoni comment)
              # Will: I just put a dummy asoc variable in at the start with all 0s, when the model is fitted using oadaFit the constraints and offsets vector are
              # automatically modified to ignore the dummy variable (a zero is appended to the end of each, or 2 zeroes if type=unconstrained is specified)

              ############# to prevent errors when nAcq is zero
              if(nAcq==0){
                learnAsoc <-learnInt<-multiInt<- naive.id <- time1 <- time2 <- stMetric <- NA
                id <- id # this will be NA by default
              } else {
                #############

                # calculate the asocial learning variables for the learning individual at each step (at each acquisition event)
                learnAsoc <- matrix(nrow=length(asoc_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnAsoc) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnInt <- matrix(nrow=length(int_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnInt) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnMulti<- matrix(nrow=length(multi_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnMulti) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnRE<- matrix(nrow=length(random_effects), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnRE) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

              } # closes else

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
                random_effects.array[,,a] <- eval(as.name(random_effects[a])) # evaluate each one in turn
              }

              if(nAcq!=0){
                for (i in 1:nAcq){ # Loop through acquisition events - i loop

                  k <- ifelse(asocialTreatment=="constant", 1, i) # this makes sure you are treating ILVs correctly if they are constant and if they are time-varying

                  for(a in 1:length(asoc_ilv)){ # Loop through asocial variables - a loop
                    learnAsoc[a,i] <- asoc_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(int_ilv)){ # Loop through asocial variables - a loop
                    learnInt[a,i] <- int_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(multi_ilv)){ # Loop through asocial variables - a loop
                    learnMulti[a,i] <- multi_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(random_effects)){ # Loop through asocial variables - a loop
                    learnRE[a,i] <- random_effects.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }

                  statusMatrix[orderAcq[i],c((i+1):(nAcq+1))] <- 1 # give the individuals that acquired the trait a status of 1 and carry skilled status (1) through to all following acquisition events

  #WH this section was wrong- I  correct it below
                  # correct the status of the individuals that can be learned from if there are ties, because in that case they will not be the same as the skilled individuals
   #               if (ties[i]==0){
   #                 availabilityMatrix[orderAcq[i],] <- statusMatrix[orderAcq[i],]
   #               } else {
   #                 availabilityMatrix[orderAcq[i],] <- ifelse(length(orderAcq[i-1]), availabilityMatrix[orderAcq[i-1],], statusMatrix[orderAcq[i],])
   #               } # closes ties if statement

                  # if the event is recorded as tied with the previous event (ties[i]==1), it means that whoever learned in the previous event cannot be learned from for this event
                  # therefore if a tie is present for event i, we do not update the availabilityMatrix to match the statusMatrix until the ties is ended
                  if (ties[i]==0){
                    availabilityMatrix[,i] <- statusMatrix[,i]
                  } else {
                    availabilityMatrix[,i] <- availabilityMatrix[,i-1]
                  } # closes ties if statement


                  #Now correct the availabilityMatrix such that individuals who are not present for an event cannot be learned from
                  availabilityMatrix<-availabilityMatrix*presenceMatrix

                  naive.id[[i]] <- which(statusMatrix[,i]==0) # index for naive individuals before the ith acquisition event

                } # closes the i loop - nAcq (k is i or 1)
                availabilityMatrix[,nAcq+1] <- statusMatrix[,nAcq+1]
              } # closes the if statement for nAcq!=0



              if(is.na(id[1])) {id <- paste(label,c(unlist(naive.id)), sep="_")} # id of naive individuals before each acquisition event, including demonstrators

              naive <- dim(assMatrix)[1]-apply(statusMatrix, 2, sum) # number of naive individuals remaining after each acq event


              # work out the number of association matrices provided and set up stMetric matrix accordingly
              stMetric <- matrix(data=0, nrow=length(id), ncol=dim(assMatrix)[3])
              dimnames(stMetric) <- list(NULL, paste("stMetric",c(1:dim(assMatrix)[3]),sep=""))

              #############################################################
              # Loop through acquisition events - learner loop
              # learnMetric is the sum of network connections of individuals that learned at each acquisition events, to other informed individuals.
              # This will always be 0 for the first animal that learned unless there were demonstrators

              # time1 and time2 index the time period or "event period" corresponding to acquisitions
              if(nAcq!=0){
                for(event in 1:nAcq){
                  # it's a shame to have two identical loops but I need time1 and time2 to be ready for use when I come to calculate the social transmission metrics below
                  time1 <- c(time1, rep(event-1, naive[event]))
                  time2 <- c(time2, rep(event, naive[event]))
                  if(is.na(event.id[1])){
                    event.id.temp <- c(event.id.temp, rep(event, each=length(naive.id[[event]])))
                  }
                } # closes for loop through events

                if(is.na(event.id[1])) {event.id <- paste(label, event.id.temp, sep="_")}

                for(event in 1:nAcq){ # event indexes the number of the acquisition event

                  #Take the appropriate association matrix from the (weighted) time varying association matrix,
                  #as determined for that event by the assMatrixIndex vector
                  assMatrix<-array(assMatrixTV[,,, assMatrixIndex[event]],dim=dim(assMatrixTV)[1:3])

                  learner <- orderAcq[event] # learner is individual id of the animal that learned AT an event
                  nonlearners <- naive.id[[event]] # nonlearners are individual id of the animals that were naive BEFORE an event
                  status <- c(status, statusMatrix[unlist(naive.id[[event]]), event+1])
                  presentInDiffusion<-c(presentInDiffusion,presenceMatrix[unlist(naive.id[[event]]), event+1])

                  temp.stMetric <- vector() # reset this before the metrics for each event are calculated

                  for (nonlearner in nonlearners){

                    # stMetric is the total assoc of the individuals that had NOT learned prior to that acquisition event, with all other already-informed individuals
                    m1 <- matrix(data=assMatrix[nonlearner,,], nrow=dim(assMatrix)[3], byrow=T) # matrix1
                    m2 <- (weights*availabilityMatrix[,event])*t(m1) # matrix2
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
              } else { # closes if(nAcq!=0) statement
                asocILVdata.naive <-matrix(data=rep(0,length(asoc_ilv)), nrow=1, ncol=length(asoc_ilv))
                intILVdata.naive<-matrix(data=rep(0,length(int_ilv)), nrow=1, ncol=length(int_ilv))
                multiILVdata.naive<-matrix(data=rep(0,length(multi_ilv)), nrow=1, ncol=length(multi_ilv))
                randomEffectdata.naive<-matrix(data=rep(0,length(random_effects)), nrow=1, ncol=length(random_effects))

                attr(asocILVdata.naive, "dimnames") <- list(NULL,asoc_ilv)
                attr(multiILVdata.naive, "dimnames") <- list(NULL,int_ilv)
                attr(intILVdata.naive, "dimnames") <- list(NULL,multi_ilv)
                attr(randomEffectdata.naive, "dimnames") <- list(NULL,random_effects)
              } # closes else


              #############################################################
              label <- rep(label, length.out=length(id))

              if(is.null(demons)) demons <- NA;

              #Subtract the first column from presenceMatrix (added previously) so it again gives the presence of each individual for each event
              presenceMatrix<-presenceMatrix[,-1]

              if(is.null(offsetCorrection)) offsetCorrection <- cbind(rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]));
              dimnames(offsetCorrection)[2]<-list(c("SocialOffsetCorr","AsocialILVOffsetCorr","InteractionOffsetCorr","MultiplicativeILVOffsetCorr"))

              callNextMethod(.Object, label=label, idname=idname, assMatrix=assMatrixTV, asoc_ilv=asoc_ilv, int_ilv=int_ilv,multi_ilv=multi_ilv,random_effects=random_effects, orderAcq=orderAcq, timeAcq=timeAcq, endTime=endTime,updateTimes=NA, ties=ties, trueTies=trueTies, demons=demons, weights=weights, statusMatrix=statusMatrix, availabilityMatrix=availabilityMatrix, event.id=event.id, id=id, time1=time1, time2=time2, status=status, presentInDiffusion= presentInDiffusion, presenceMatrix = presenceMatrix ,asocialTreatment=asocialTreatment, stMetric=stMetric, asocILVdata=asocILVdata.naive, intILVdata=intILVdata.naive, multiILVdata=multiILVdata.naive,randomEffectdata=randomEffectdata.naive,offsetCorrection=offsetCorrection,assMatrixIndex=assMatrixIndex)
            }else
            {
              #TADA version to be inserted here

              #put the time varying association matrix into an object called assMatrixTV
              if(length(dim(assMatrix))==3){ assMatrixTV<- array(data=assMatrix,dim=c(dim(assMatrix),1))}else{assMatrixTV<-assMatrix}

              # create default numeric id vector for individuals if no names are provided, if it is, convert to a factor
              if(is.null(idname)) idname <- (1:dim(assMatrix)[1]);
              # if there are no ties, make a vector of zeroes
              if(is.null(ties)) ties <- rep(0,length(orderAcq));



              # each asoc_ vector should be a vector of character strings of matrices whose names correspond to the names of the individual-level variables (ILVs). the rows of each matrix should correspond to individuals and the columns to times at which the value of the ILV changes
              nAcq <- ifelse(any(is.na(orderAcq)), 0, length(orderAcq)) # number of acquisition events EXCLUDING demonstrators.
              time1 <- vector(); # time period index: time start
              time2 <- vector(); # time period index: time end
              event.id.temp <- vector(); # vector to indicate the event number (this is essentially indexing the naive individuals before each event)

              #If no asoc variable is provided, set a single column of zeroes
              if(asoc_ilv[1]=="ILVabsent"|int_ilv[1]=="ILVabsent"|multi_ilv[1]=="ILVabsent"){
                ILVabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }
              if(random_effects[1]=="REabsent"){
                REabsent <-matrix(data = rep(0, dim(assMatrix)[1]), nrow=dim(assMatrix)[1], byrow=F)
              }

              totalMetric <- vector() # total association of the individual that DID learn at an acquisition event, with all other individuals
              learnMetric <- vector(); # total associations of the individual that DID learn at an acquisition event, with all other individuals that have already learned

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
                  asoc_ilv.array <- array(dim=c(asoc_ilv.dim,(nAcq+1),length(asoc_ilv)))
                  dimnames(asoc_ilv.array) <- list(NULL, c(paste("time",c(1:(nAcq+1)),sep="")), asoc_ilv)
                  int_ilv.array <- array(dim=c(int_ilv.dim,(nAcq+1),length(int_ilv)))
                  dimnames(int_ilv.array) <- list(NULL, c(paste("time",c(1:(nAcq+1)),sep="")), int_ilv)
                  multi_ilv.array <- array(dim=c(multi_ilv.dim,(nAcq+1),length(multi_ilv)))
                  dimnames(multi_ilv.array) <- list(NULL, c(paste("time",c(1:(nAcq+1)),sep="")), multi_ilv)
                  random_effects.array <- array(dim=c(random_effects.dim,(nAcq+1),length(random_effects)))
                  dimnames(random_effects.array) <- list(NULL, c(paste("time",c(1:(nAcq+1)),sep="")), random_effects)
                }
              }

              # generate a matrix that contains the status of each individual at each acquisition event
              statusMatrix <- matrix(0, nrow=dim(assMatrix)[2], ncol=1+nAcq)  # a matrix with as many rows as indivs and as many columns as acquisition events PLUS one for the demonstrators
              # create a list vector to hold the index of naive individuals after each acquisition event
              naive.id <- vector(mode="list", length=nAcq+1)

              # if there are seeded demonstrators add the vector (which should have length dim(assMatrix)[1]) to the first column of the statusMatrix to show which individuals set out as skilled (status of 1)
              if(is.null(demons)){
                statusMatrix[,1] <- rep(0,dim(assMatrix)[1])
              } else {
                for(i in 1:(1+nAcq)){
                  statusMatrix[,i] <- demons
                }
              }

              availabilityMatrix <- statusMatrix # if there are ties the statusMatrix and the availabilityMatrix will differ (the latter gives who is available to be learned *from*). we create it as identical and modify it accordingly below

              #presenceMatrix gives who was present in the diffusion for each event- set to 1s by default
              if(is.null(presenceMatrix)){
                presenceMatrix<-statusMatrix;
                presenceMatrix[,]<-1;
              }else{
                #Add a column to the start of the presenceMatrix so the dimensions match statusMatrix
                presenceMatrix<-cbind(presenceMatrix[,1], presenceMatrix)
              }


              asocILVdata.naive <-intILVdata.naive<-multiILVdata.naive<-randomEffectdata.naive<-vector() # this will hold the individual level variables for the naive individuals

              # YOU WILL HAVE TO MAKE THIS WORK EVEN WHEN THERE ARE NO ASOCIAL VARIABLES... THINK ABOUT HOW YOU MIGHT DO THIS 20120822 (Theoni comment)
              # Will: I just put a dummy asoc variable in at the start with all 0s, when the model is fitted using oadaFit the constraints and offsets vector are
              # automatically modified to ignore the dummy variable (a zero is appended to the end of each, or 2 zeroes if type=unconstrained is specified)

              ############# to prevent errors when nAcq is zero
              if(nAcq==0){
                learnAsoc <-learnInt<-multiInt<- NA
                id <- id # this will be NA by default
              } else {
                #############

                # calculate the asocial learning variables for the learning individual at each step (at each acquisition event)
                learnAsoc <- matrix(nrow=length(asoc_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnAsoc) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnInt <- matrix(nrow=length(int_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnInt) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnMulti<- matrix(nrow=length(multi_ilv), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnMulti) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

                learnRE<- matrix(nrow=length(random_effects), ncol=nAcq) # a matrix with as many rows as individuals and as many columns as acquisition events
                dimnames(learnRE) <- list(NULL, c(paste("id",c(orderAcq),"time",c(1:nAcq),sep="")))

              } # closes else

              if(dim(eval(as.name(asoc_ilv[1])))[2]==(nAcq+1)|dim(eval(as.name(asoc_ilv[1])))[2]==1){
                for(a in 1:length(asoc_ilv)){ # Loop through asocial variables - a loop
                  asoc_ilv.array[,,a] <- eval(as.name(asoc_ilv[a])) # evaluate each one in turn
                }
              }else{
                # If no values are provided for the final period to endTime, the values are assumed to be the same as for the final event
                for(a in 1:length(asoc_ilv)){ # Loop through asocial variables - a loop
                  asoc_ilv.array[,1:nAcq,a] <- eval(as.name(asoc_ilv[a])) # evaluate each one in turn
                  asoc_ilv.array[,nAcq+1,a] <- asoc_ilv.array[,nAcq,a]
                }
              }


              if(dim(eval(as.name(int_ilv[1])))[2]==(nAcq+1)|dim(eval(as.name(int_ilv[1])))[2]==1){
                for(a in 1:length(int_ilv)){ # Loop through asocial variables - a loop
                  int_ilv.array[,,a] <- eval(as.name(int_ilv[a])) # evaluate each one in turn
                }
              }else{
                # If no values are provided for the final period to endTime, the values are assumed to be the same as for the final event
                for(a in 1:length(int_ilv)){ # Loop through asocial variables - a loop
                  int_ilv.array[,1:nAcq,a] <- eval(as.name(int_ilv[a])) # evaluate each one in turn
                  int_ilv.array[,nAcq+1,a] <- int_ilv.array[,nAcq,a]
                }
              }

              if(dim(eval(as.name(multi_ilv[1])))[2]==(nAcq+1)|dim(eval(as.name(multi_ilv[1])))[2]==1){
                for(a in 1:length(multi_ilv)){ # Loop through asocial variables - a loop
                  multi_ilv.array[,,a] <- eval(as.name(multi_ilv[a])) # evaluate each one in turn
                }
              }else{
                # If no values are provided for the final period to endTime, the values are assumed to be the same as for the final event
                for(a in 1:length(multi_ilv)){ # Loop through asocial variables - a loop
                  multi_ilv.array[,1:nAcq,a] <- eval(as.name(multi_ilv[a])) # evaluate each one in turn
                  multi_ilv.array[,nAcq+1,a] <- multi_ilv.array[,nAcq,a]
                }
              }

              if(dim(eval(as.name(random_effects[1])))[2]==(nAcq+1)|dim(eval(as.name(random_effects[1])))[2]==1){
                for(a in 1:length(random_effects)){ # Loop through asocial variables - a loop
                  random_effects.array[,,a] <- eval(as.name(random_effects)) # evaluate each one in turn
                }
              }else{
                # If no values are provided for the final period to endTime, the values are assumed to be the same as for the final event
                for(a in 1:length(random_effects)){ # Loop through asocial variables - a loop
                  random_effects.array[,1:nAcq,a] <- eval(as.name(random_effects[a])) # evaluate each one in turn
                  random_effects.array[,nAcq+1,a] <- random_effects.array[,nAcq,a]
                }
              }



              #  if(nAcq!=0){
              for (i in 1:(nAcq+1)){ # Loop through acquisition events - i loop

                if(i<=nAcq){
                  #exclude final period where no one learned
                  k <- ifelse(asocialTreatment=="constant", 1, i) # this makes sure you are treating ILVs correctly if they are constant and if they are time-varying

                  for(a in 1:length(asoc_ilv)){ # Loop through asocial variables - a loop
                    learnAsoc[a,i] <- asoc_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(int_ilv)){ # Loop through asocial variables - a loop
                    learnInt[a,i] <- int_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(multi_ilv)){ # Loop through asocial variables - a loop
                    learnMulti[a,i] <- multi_ilv.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }
                  for(a in 1:length(random_effects)){ # Loop through asocial variables - a loop
                    learnRE[a,i] <- random_effects.array[orderAcq[i],k,a] # fill the matrix with the rows of the asoc matrix that correspond to the individual that learned at each acquisition event. each column is an individual, each row is a type of variable
                  }

                  statusMatrix[orderAcq[i],c((i+1):(nAcq+1))] <- 1 # give the individuals that acquired the trait a status of 1 and carry skilled status (1) through to all following acquisition events


                  # if the event is recorded as tied with the previous event (ties[i]==1), it means that whoever learned in the previous event cannot be learned from for this event
                  # therefore if a tie is present for event i, we do not update the availabilityMatrix to match the statusMatrix until the ties is ended
                  if (ties[i]==0){
                     availabilityMatrix[,i] <- statusMatrix[,i]
                  } else {
                     availabilityMatrix[,i] <- availabilityMatrix[,i-1]
                  } # closes ties if statement
                }


                #Now correct the availabilityMatrix such that individuals who are not present for an event cannot be learned from
                availabilityMatrix<-availabilityMatrix*presenceMatrix

                naive.id[[i]] <- which(statusMatrix[,i]==0) # index for naive individuals before the ith acquisition event


              } # closes the i loop - nAcq (k is i or 1)
              availabilityMatrix[,nAcq+1] <- statusMatrix[,nAcq+1]
              #  } # closes the if statement for nAcq!=0


              if(is.na(id[1])) {id <- paste(label,c(unlist(naive.id)), sep="_")} # id of naive individuals before each acquisition event, including demonstrators

              naive <- dim(assMatrix)[1]-apply(statusMatrix, 2, sum) # number of naive individuals remaining after each acq event (last event will be the end of the diffusion for TADA data with incomplete diffusion)


              # work out the number of association matrices provided and set up stMetric matrix accordingly
              stMetric <- matrix(data=0, nrow=length(id), ncol=dim(assMatrix)[3])
              dimnames(stMetric) <- list(NULL, paste("stMetric",c(1:dim(assMatrix)[3]),sep=""))

              #############################################################
              # Loop through acquisition events - learner loop
              # learnMetric is the sum of network connections of individuals that learned at each acquisition events, to other informed individuals.
              # This will always be 0 for the first animal that learned unless there were demonstrators

              # time1 and time2 index the time period or "event period" corresponding to acquisitions
              #if(nAcq!=0){ I cut this since it should work without now I have changed nAcq to nAcq+1 for TADA
              for(event in 1:(nAcq+1)){
                # it's a shame to have two identical loops but I need time1 and time2 to be ready for use when I come to calculate the social transmission metrics below
                time1 <- c(time1, rep(event-1, naive[event]))
                time2 <- c(time2, rep(event, naive[event]))
                if(is.na(event.id[1])){
                  event.id.temp <- c(event.id.temp, rep(event, each=length(naive.id[[event]])))
                }
              } # closes for loop through events

              timeAcqNew<-c(0,timeAcq,endTime)
              TADAtime1<-timeAcqNew[time1+1]
              TADAtime2<-timeAcqNew[time2+1]


              if(is.na(event.id[1])) {event.id <- paste(label, event.id.temp, sep="_")}

              for(event in 1:(nAcq+1)){ # event indexes the number of the acquisition event- increased by 1 for TADA to allow for the endTime period

                #Take the appropriate association matrix from the (weighted) time varying association matrix,
                #as determined for that event by the assMatrixIndex vector
                if((length(assMatrixIndex)==nAcq)&(event==(nAcq+1))){
                  #If no separate assMatrix is specified for the end period it is assumed to be the same as for the final acquisition event
                  assMatrix<-array(assMatrixTV[,,, assMatrixIndex[nAcq]],dim=dim(assMatrixTV)[1:3])
                }else{
                  assMatrix<-array(assMatrixTV[,,, assMatrixIndex[event]],dim=dim(assMatrixTV)[1:3])
                }

                learner <- orderAcq[event] # learner is individual id of the animal that learned AT an event
                nonlearners <- naive.id[[event]] # nonlearners are individual id of the animals that were naive BEFORE an event

                if(length(nonlearners)>0){
                  #If everyone has learned by the final period of a TADA (i.e. up to endTime) the next section triggers errors- and we do not need a final period

                   status <- c(status, statusMatrix[unlist(naive.id[[event]]), min(nAcq+1,event+1)])
                   presentInDiffusion<-c(presentInDiffusion,presenceMatrix[unlist(naive.id[[event]]), min(nAcq+1,event+1)])

                   temp.stMetric <- vector() # reset this before the metrics for each event are calculated

                   for (nonlearner in nonlearners){

                     # stMetric is the total assoc of the individuals that did NOT learn by that acquisition event, with all other already-informed individuals
                     m1 <- matrix(data=assMatrix[nonlearner,,], nrow=dim(assMatrix)[3], byrow=T) # matrix1
                     m2 <- (weights*availabilityMatrix[,event])*t(m1) # matrix2
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
                }#closes if(length(nonlearners)>0) loop

              } # closes event loop

              #############################################################
              label <- rep(label, length.out=length(id))

              if(is.null(demons)) demons <- NA;

              #Subtract the first column from presenceMatrix (added previously) so it again gives the presence of each individual for each event
              presenceMatrix<-presenceMatrix[,-1]

              if(is.null(offsetCorrection)) offsetCorrection <- cbind(rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]),rep(0,dim(asocILVdata.naive)[1]));
              dimnames(offsetCorrection)[2]<-list(c("SocialOffsetCorr","AsocialILVOffsetCorr","InteractionOffsetCorr","MultiplicativeILVOffsetCorr"))

              callNextMethod(.Object, label=label, idname=idname, assMatrix=assMatrixTV, asoc_ilv=asoc_ilv, int_ilv=int_ilv,multi_ilv=multi_ilv,random_effects=random_effects, orderAcq=orderAcq, timeAcq=timeAcq, endTime=endTime,updateTimes=NA, ties=ties, trueTies=trueTies, demons=demons, weights=weights, statusMatrix=statusMatrix, availabilityMatrix=availabilityMatrix, event.id=event.id, id=id, time1=time1, time2=time2,TADAtime1=TADAtime1, TADAtime2=TADAtime2, status=status, presentInDiffusion= presentInDiffusion, presenceMatrix = presenceMatrix ,asocialTreatment=asocialTreatment, stMetric=stMetric, asocILVdata=asocILVdata.naive, intILVdata=intILVdata.naive, multiILVdata=multiILVdata.naive,randomEffectdata=randomEffectdata.naive,offsetCorrection=offsetCorrection,assMatrixIndex=assMatrixIndex)

            }
          } # end function

) # end setMethod "initialize"


#'Create an nbdaData object
#'
#'\code{nbdaData} creates an object of class nbdaData required for conducting a network-based diffusion analysis (NBDA) using
#'\code{\link{oadaFit}} or fitting a continuous TADA using \code{\link{tadaFit}}.
#'
#'An nbdaData object is required to contain the data for each diffusion to be used in the network based diffusion analysis
#'when using the order of acquisition diffusion analysis (OADA) method (\code{\link{oadaFit}} function) or continuous time of
#'acquistion diffusion analysis (cTADA) method (\code{\link{tadaFit}} function). When multiple diffusions are being modelled
#' a list of nbdaData objects is provided to the \code{\link{oadaFit}} or \code{\link{tadaFit}} function.
#'
#'@seealso To create a data object for discrete TADA models use \code{\link{dTADAData}}.
#'
#'@param label a string giving a name for the diffusion
#'@param idname an optional vector giving an id for each individual in the diffusion. Order to match that used in assMatrix.
#'@param assMatrix a three or four dimensional array specifying the association matrices or social network(s) to be included
#'in the analysis.
#'If assMatrix has three dimensions the networks are assumed to be static. The first dimension is the row of the social
#'network(s) and the second dimension is the column. Connections are taken to be from the column to the row, such that the
#'entry in row i and column j gives the connection from j to i. The third dimension contains the different networks to be
#'included in the analysis, i.e. assMatrix[,,1] gives the first network, assMatrix[,,2] the second and so on. If a fourth
#'dimension is included the network is taken to vary over time, with assMatrix[,,k,1] giving the values for network k in
#'time period 1, assMatrix[,,k,2] giving the values for network k in time period 2 and so on (see assMatrixIndex below).
#'@param asoc for backwards compatibility only, now replaced with asoc_ilv below.
#'@param asoc_ilv an optional character vector giving the names of individual-level variables (ILVs) having an effect on the
#'rate of asocial learning. Each entry in asoc_ilv must refer to an object containing the data for that variable. If
#'asocialTreatment= "constant", each variable must be a single column matrix with rows equal to the number of rows in
#'assMatrix, thus providing a single value for each individual in the diffusion. If asocialTreatment= "timevarying" then
#'each variable must be a matrix with columns equal to the number of acquisition event, and rows equal to the number of
#'rows in  assMatrix, thus providing a  value for each individual in the diffusion at the time of each acquisition event.
#'@param int_ilv a optional character vector giving the names of individual-level variables (ILVs) having an (interactive)
#'effect on the rate of social learning. These are specified as for asoc_ilv.
#'@param multi_ilv a optional character vector giving the names of individual-level variables (ILVs) having the same effect
#'on both asocial and social learning (referred to as the multiplicative model). These are specified as for asoc_ilv.
#'@param random_effects a optional character vector giving the names of random effects. Each variable must be a single
#'column matrix with rows equal to the number of rows in assMatrix, thus providing a single level for each individual in the
#'diffusion. These only have an effect when using \code{\link{oadaFit}}, and are assumed to operate equally on social and asocial
#'learning. If more complex effects are required, or random effects are required for a TADA then a Bayesian approach is
#'recommended (not implemented in the NBDA package).
#'@param orderAcq a numerical vector giving the order in which individuals acquired the target behaviour, with numbers
#'referring to rows of assMatrix.
#'@param timeAcq a numerical vector giving the time at which each individual acquired the target behaviour, given in the
#'order matching orderAcq. This is necessary for conducting a TADA with \code{\link{tadaFit}} but not if conducting an OADA with
#'\code{\link{oadaFit}}.
#'@param endTime numeric giving the time at which the diffusion ended. This is necessary for conducting a TADA with
#'\code{\link{tadaFit}} but not if conducting an OADA with \code{\link{oadaFit}}.
#'@param ties numeric binary vector specifying if each acquisition was "tied" with the acquisition before. e.g.
#'c(0,1,1,0,0,1,1) specifies that events 2 and 3 are tied, as are events 6 and 7. Events should be specified as ties if they
#'occurred too close in time for the individuals in question to have plausibly have learned from one another(in contrast to
#'trueTies below).
#'@param trueTies a list of numeric vectors specifying which events are trueTies, i.e. where we do not know the order in
#'which the events occurred. e.g. list(c(1,2),c(10,11,12)) specifies that we do not know the real order of events 1 and 2,
#'or of events 10, 11 and 12. This is only applicable to OADA using \code{\link{oadaFit}}. trueTies are accounted for by adding the
#'likelihood across all orders consistent with the trueTies, so can be highly computationally intensive for large trueTies or
#'large numbers of trueTies.
#'@param updateTimes non-functioning argument- can be ignored.
#'@param demons an optional binary numeric vector specifying which individuals are trained demonstrators or had otherwise already
#'acquired the target behaviour prior to the start of the diffusion. Length should match the number of rows of assMatrix.
#'e.g. c(0,0,1,0,0,0,1) specifies that individuals 3 and 7 are trained demonstrators.
#'@param presenceMatrix an optional binary matrix specifying who was present in the diffusion for each event- set to 1s by default.
#'Number of rows to match the number of individuals, and the number of columns to match the number of events (length of
#'orderAcq). 1 denotes that an individual was present in the diffusion for a given event, 0 denotes that an individual was
#'absent, and so could neither learn nor transmit the behaviour to others for that event. Set to 1s by default.
#'@param assMatrixIndex a numeric vector necessary if time-varying networks are used, i.e. if assMatrix is a four dimensional
#'array. This specifies which time period (4th dimension of assMatrix) is applicable to each event.
#'e.g. assMatrixIndex=c(1,1,1,2,2,2,3,3,3) specifies that assMatrix[,,,1] gives the networks for events 1-3, assMatrix[,,,2]
#'gives the networks for events 4-6 and assMatrix[,,,3] gives the networks for events 7-9.
#'@param weights an optional numeric vector giving the transmission weights for each individual, length to match the number of rows in
#'assMatrix. It is assumed that the rate at which information is transmitted FROM individual j is multiplied by entry j in
#'the weights vector.
#'@param asocialTreatment a string- "constant" if ILVs are assumed to be constant over time and "timevarying" if they are
#'assumed to be different for each acquistision event. See asoc_ilv, int_ilv and multi_ilv above.
#'@param offsetCorrection a four column matrix with rows equal to the number of individuals, specifying the offset to be
#'added to the social learning component and the linear predictors for the effect of ILVs on asocial learning, social
#'learning, or on both (multiplicative effect). This is used by other functions calling the \code{nbdaData} function and it
#'is not anticipated that users will need to use this argument directly.
#'@return Returns an object of class nbdaData, which can be used to conduct an NBDA using \code{\link{oadaFit}} or \code{\link{tadaFit}}
#'functions.


nbdaData <- function(label, idname=NULL, assMatrix, asoc=NULL, asoc_ilv="ILVabsent",int_ilv="ILVabsent",multi_ilv="ILVabsent",random_effects="REabsent", orderAcq, timeAcq=NA, endTime=max(timeAcq)+1,ties=NULL, trueTies=list(NULL), updateTimes=NULL, demons=NULL, presenceMatrix =NULL,assMatrixIndex= rep(1,length(orderAcq)), weights=rep(1, dim(assMatrix)[1]), asocialTreatment="constant",offsetCorrection=NULL){

  # For backwards compatibility, allow the asoc argument and assume it refers to asoc_ilv

  if(!is.null(asoc)&asoc_ilv[1]=="ILVabsent") asoc_ilv<-asoc

  new("nbdaData",label=label,idname=idname,assMatrix=assMatrix,availabilityMatrix=availabilityMatrix,asoc_ilv=asoc_ilv,int_ilv=int_ilv,multi_ilv=multi_ilv,random_effects=random_effects,orderAcq=orderAcq,timeAcq=timeAcq, endTime=endTime,ties=ties,trueTies=trueTies,demons=demons, presenceMatrix = presenceMatrix, assMatrixIndex = assMatrixIndex,weights=weights,asocialTreatment=asocialTreatment);

}

