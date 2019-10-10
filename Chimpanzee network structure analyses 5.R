#This script performs analyses of the structure of the observation network to test for preferences in who observed whom, based
#on kin, rank and age.

########################################################################
#LOAD IN THE DATA FILES
########################################################################

#Load in group 1 data
Grp1_Data<-read.csv("Grp1_Data.csv")
#Load in session time data for group 1
Grp1_Time<-read.csv("Grp1_Time_overview_WH.csv")
#Calculate the cumulative time at the start of each session
Grp1_Time$cumulativeStartTime<-c(0,cumsum(Grp1_Time$DurationInS)[-length(Grp1_Time$DurationInS)])
#Get the cumulative time across sessions for each event (in case we use TADA or stratified OADA)
Grp1_Data$cumulativeTime<-Grp1_Time$cumulativeStartTime[Grp1_Data$Session]+Grp1_Data$TimeSinceStartSession
#Load the association matrix
Grp1_association_matrix<-read.table("Grp1_association matrix.txt")
#Get a list of names from the association matrix
Grp1_indList<-row.names(Grp1_association_matrix)

#This returns the order of first solving
Grp1_oa_asText<-as.character(unique(Grp1_Data$solver))
#Convert to numerical vector
Grp1_oa<-rep(NA,length(Grp1_oa_asText))
for(i in 1:length(Grp1_oa))Grp1_oa[i]<-which(Grp1_oa_asText[i]==Grp1_indList)

#Get times of first solving
Grp1_ta<-sort(tapply(Grp1_Data$cumulativeTime,Grp1_Data$solver,min))
cbind(Grp1_oa,Grp1_ta)

#Load in group 2 data
Grp2_Data<-read.csv("Grp2_Data.csv")
#Load in session time data for group 2
Grp2_Time<-read.csv("Grp2_Time_overview_WH.csv")
#Calculate the cumulative time at the start of each session
Grp2_Time$cumulativeStartTime<-c(0,cumsum(Grp2_Time$DurationInS)[-length(Grp2_Time$DurationInS)])
#Get the cumulative time across sessions for each event (in case we use TADA or stratified OADA)
Grp2_Data$cumulativeTime<-Grp2_Time$cumulativeStartTime[Grp2_Data$Session]+Grp2_Data$TimeSinceStartSession
#Load the association matrix
Grp2_association_matrix<-read.table("Grp2_association matrix.txt")
#Get a list of names from the association matrix
Grp2_indList<-row.names(Grp2_association_matrix)


########################################################################
#GET ORDER OF ACQUISITION ACROSS GROUPS (This is just so we can use the same code to get the observation network below)
########################################################################

#This returns the order of first solving
Grp2_oa_asText<-as.character(unique(Grp2_Data$solver))
#Convert to numerical vector
Grp2_oa<-rep(NA,length(Grp2_oa_asText))
for(i in 1:length(Grp2_oa))Grp2_oa[i]<-which(Grp2_oa_asText[i]==Grp2_indList)
Grp2_oa
#Get times of first solving
Grp2_ta<-sort(tapply(Grp2_Data$cumulativeTime,Grp2_Data$solver,min))
cbind(Grp2_oa,Grp2_ta)

#Get order and time across both groups for the stratified OADA
both_oa<-c(Grp1_oa,Grp2_oa+length(Grp1_indList))
both_ta<-c(Grp1_ta,Grp2_ta)
both_oa<-both_oa[order(both_ta)]
both_ta<-both_ta[order(both_ta)]
cbind(both_oa,both_ta)

########################################################################
#CREATE OBSERVATION NETWORK
########################################################################

#Now create dynamic observation network for group 1
Grp1_observation_network<-array(0,dim=c(dim((Grp1_association_matrix)),2,length(both_oa)))
#Cycle through the data building up the number of observations prior between members of each dyad prior to each event
event<-1
for(i in 1:dim(Grp1_Data)[1]){
  if(Grp1_Data$cumulativeTime[i]>=both_ta[event]) event<-event+1
  Grp1_observation_network[as.character(Grp1_Data$observers[i])==Grp1_indList,as.character(Grp1_Data$solver[i])==Grp1_indList,1,event:length(both_oa)]<-Grp1_observation_network[as.character(Grp1_Data$observers[i])==Grp1_indList,as.character(Grp1_Data$solver[i])==Grp1_indList,1,event:length(both_oa)]+1
}
#Look at last event as a check
Grp1_observation_network[,,1,length(both_oa)]
#I notice individuals are often listed as observing their own solves? This doe not affect the analysis but I will set the diagonal to zero so we can see the connections more clearly
for (i in 1:length(both_oa))diag(Grp1_observation_network[,,1,i])<-0
Grp1_observation_network[,,1,length(both_oa)]

#Note that the identity of the observed individual does not matter to fit this kind of NBDA. But if we keep the information in network form we can test for biases in learning later if required.
#Also the NBDA package requires the network form for the data.

#Now create dynamic observation network for group 2
Grp2_observation_network<-array(0,dim=c(dim((Grp2_association_matrix)),2,length(both_oa)))
#Cycle through the data building up the number of observations prior between members of each dyad prior to each event
event<-1
for(i in 1:dim(Grp2_Data)[1]){
  if(event<=length(both_oa)){
    if(Grp2_Data$cumulativeTime[i]>=both_ta[event]) event<-event+1
  }
if(event<=length(both_oa)){
    Grp2_observation_network[as.character(Grp2_Data$observers[i])==Grp2_indList,as.character(Grp2_Data$solver[i])==Grp2_indList,2,event:length(both_oa)]<-Grp2_observation_network[as.character(Grp2_Data$observers[i])==Grp2_indList,as.character(Grp2_Data$solver[i])==Grp2_indList,2,event:length(both_oa)]+1
  }
}
#Look at last event as a check
Grp2_observation_network[,,2,length(both_oa)]
#I notice individuals are often listed as observing their own solves? This does not affect the analysis but I will set the diagonal to zero so we can see the connections more clearly
for (i in 1:length(both_oa))diag(Grp2_observation_network[,,2,i])<-0
Grp2_observation_network[,,2,length(both_oa)]

#Both networks contain some high numbers so I will divide them by 10 to aid model fitting
Grp1_observation_network<-Grp1_observation_network/10
Grp2_observation_network<-Grp2_observation_network/10

#combine the networks into a big one with 0 connections between groups

combined_observation_network<-array(0,dim=c(62,62,2,16))
combined_observation_network[1:21,1:21,1,]<-Grp1_observation_network[,,1,]
combined_observation_network[22:62,22:62,2,]<-Grp2_observation_network[,,2,]

########################################################################
#ORGANIZE ILVs
########################################################################

#Now extract the ILVs for each group, rank and age
Grp1_ILVs<-read.csv("Grp1_ILVs.csv")
Grp2_ILVs<-read.csv("Grp2_ILVs.csv")
names(Grp1_ILVs)[1]<-"name"
#I will standardize age, but then set so zero= oldest individual to aid model convergence
Grp1_ILVs$stAge<-(Grp1_ILVs$AgeAtStart-max(c(Grp1_ILVs$AgeAtStart,Grp2_ILVs$AgeAtStart)))/sd(c(Grp1_ILVs$AgeAtStart,Grp2_ILVs$AgeAtStart))
#I will transform rank so the highest ranked individual is 0 and the lowest ranked is 1
Grp1_ILVs$transRank<-(Grp1_ILVs$Rank-1)/(max(Grp1_ILVs$Rank)-1)



names(Grp2_ILVs)[1]<-"name"
Grp2_ILVs$stAge<-(Grp2_ILVs$AgeAtStart-max(c(Grp1_ILVs$AgeAtStart,Grp2_ILVs$AgeAtStart)))/sd(c(Grp1_ILVs$AgeAtStart,Grp2_ILVs$AgeAtStart))
Grp2_ILVs$transRank<-(Grp2_ILVs$Rank-1)/(max(Grp2_ILVs$Rank)-1)

stAge1<-cbind(Grp1_ILVs$stAge)
transRank1<-cbind(Grp1_ILVs$transRank)
male1<-cbind(Grp1_ILVs$Sex=="male")*1

stAge2<-cbind(Grp2_ILVs$stAge)
transRank2<-cbind(Grp2_ILVs$transRank)
male2<-cbind(Grp2_ILVs$Sex=="male")*1

stAge<-rbind(stAge1,stAge2)
transRank<-rbind(transRank1,transRank2)
male<-rbind(male1,male2)

ilvs<-c("stAge","transRank","male")



Grp1_maternal_kin<-as.matrix(read.table("Grp1_maternal kin.txt"))

#I notice we have some different individuals here to in the association network (I used the latter to generate the list of individuals in the
#NBDA and observation network). Here I will cut them down to indivuals in both
#> rownames(Grp1_association_matrix)
#[1] "bj"  "bob" "boo" "bre" "chr" "gen" "ger" "gir" "gon" "ils" "ing" "inn" "ire" "pal" "rac" "reg" "ren" "rit" "rus" "tar" "tob"
#> rownames(Grp1_maternal_kin)
#[1] "bj"  "bob" "bre" "chr" "gen" "gon" "ian" "ils" "ing" "inn" "ire" "pal" "rac" "reg" "ren" "rit" "rus" "tar" "tob"


Grp1_observationNet<-Grp1_observation_network[,,1,16]
dimnames(Grp1_observationNet)[[1]]<-Grp1_indList
dimnames(Grp1_observationNet)[[2]]<-Grp1_indList

Grp1_observationNet<-Grp1_observationNet[Grp1_indList%in%dimnames(Grp1_maternal_kin)[[1]],Grp1_indList%in%dimnames(Grp1_maternal_kin)[[1]]]
Grp1_maternal_kin<-Grp1_maternal_kin[rownames(Grp1_maternal_kin)%in%Grp1_indList,rownames(Grp1_maternal_kin)%in%Grp1_indList]
dimnames(Grp1_maternal_kin)[[1]]==dimnames(Grp1_observationNet)[[1]]
Grp1_newIndList<-dimnames(Grp1_observationNet)[[1]]

Grp2_maternal_kin<-as.matrix(read.table("Grp2_maternal kin.txt"))

Grp2_observationNet<-Grp2_observation_network[,,2,16]
dimnames(Grp2_observationNet)[[1]]<-Grp2_indList
dimnames(Grp2_observationNet)[[2]]<-Grp2_indList

Grp2_observationNet<-Grp2_observationNet[Grp2_indList%in%dimnames(Grp2_maternal_kin)[[1]],Grp2_indList%in%dimnames(Grp2_maternal_kin)[[1]]]
Grp2_maternal_kin<-Grp2_maternal_kin[rownames(Grp2_maternal_kin)%in%Grp2_indList,rownames(Grp2_maternal_kin)%in%Grp2_indList]
dimnames(Grp2_maternal_kin)[[1]]==dimnames(Grp2_observationNet)[[1]]
Grp2_newIndList<-dimnames(Grp2_observationNet)[[1]]

#I will convert the observation network to proportion observed, otherwise individuals who performed many times will have a disproportionate
#effect on the analysis

#Get total manipulations
Grp1_totals<-matrix(0,length(Grp1_newIndList),length(Grp1_newIndList))
Grp1_solverTable<-table(Grp1_Data$solver)
Grp1_solverTable
for(i in 1:length(Grp1_solverTable)) Grp1_totals[,Grp1_newIndList==names(Grp1_solverTable)[i]]<-Grp1_solverTable[i]
rbind(Grp1_newIndList,Grp1_totals)
Grp1_proportionObserved<-Grp1_observationNet*10/Grp1_totals
Grp1_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp1_proportionObserved[1,])
Grp1_maternal_kin<-Grp1_maternal_kin[,!is.nan(Grp1_proportionObserved[1,])]
Grp1_proportionObserved<-Grp1_proportionObserved[,!is.nan(Grp1_proportionObserved[1,])]



Grp2_totals<-matrix(0,length(Grp2_newIndList),length(Grp2_newIndList))
Grp2_solverTable<-table(Grp2_Data$solver)
Grp2_solverTable
for(i in 1:length(Grp2_solverTable)) Grp2_totals[,Grp2_newIndList==names(Grp2_solverTable)[i]]<-Grp2_solverTable[i]
rbind(Grp2_newIndList,Grp2_totals)
Grp2_proportionObserved<-Grp2_observationNet*10/Grp2_totals
Grp2_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp2_proportionObserved[1,])
Grp2_maternal_kin<-Grp2_maternal_kin[,!is.nan(Grp2_proportionObserved[1,])]
Grp2_proportionObserved<-Grp2_proportionObserved[,!is.nan(Grp2_proportionObserved[1,])]


#t-test confirms more observations between kin. But of course the p value it not valid
t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_maternal_kin),Grp2_maternal_kin))

#Welch Two Sample t-test
#
#data:  c(as.vector(Grp1_proportionObserved), as.vector(Grp2_proportionObserved)) by c(as.vector(Grp1_maternal_kin), Grp2_maternal_kin)
#t = -2.7218, df = 38.09, p-value = 0.009736
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.07103722 -0.01044105
#sample estimates:
#  mean in group 0 mean in group 1 
#0.01222363      0.05296276 

#Repeat for each group separately

t.test(c(as.vector(Grp1_proportionObserved))~c(as.vector(Grp1_maternal_kin)))
#Welch Two Sample t-test
#
#data:  c(as.vector(Grp1_proportionObserved)) by c(as.vector(Grp1_maternal_kin))
#t = -1.9282, df = 15.241, p-value = 0.07268
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0097461465  0.0004813747
#sample estimates:
#  mean in group 0 mean in group 1 
#     0.01012828      0.05645214 

#Very similar result from just group 1 (but with less power)

t.test(c(as.vector(Grp2_proportionObserved))~c(as.vector(Grp2_maternal_kin)))
#Welch Two Sample t-test
#
#data:  c(as.vector(Grp2_proportionObserved)) by c(as.vector(Grp2_maternal_kin))
#t = -1.9552, df = 20.97, p-value = 0.06402
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0077759142  0.0002401425
#sample estimates:
#  mean in group 0 mean in group 1 
#     0.01262533      0.05030419

#Very similar result from just group 2 (but with less power again)

#coefficients in each group are surprisingly similar!


testStat<-t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_maternal_kin),Grp2_maternal_kin))$statistic
testStatGrp1<-t.test(as.vector(Grp1_proportionObserved)~as.vector(Grp1_maternal_kin))$statistic
testStatGrp2<-t.test(as.vector(Grp2_proportionObserved)~as.vector(Grp2_maternal_kin))$statistic
                                          

#Get the observation net back to no.observation (it is divided by 10)
Grp1_observationNet_Perm<-Grp1_observationNet*10
Grp2_observationNet_Perm<-Grp2_observationNet*10


#The following procedure switches a random observation between two individuals for one solver and does the reverse switch for another solver
#This randomises the network whilst keeping the total number of solves and observations by each individual the same
#I then recalculate the proportion observed and re-run the t test, and repeat to get the null distribution

#This takes a long time, so if you want to run consider setting noPerm to 1000

number_of_switches<-10000
noPerm<-10000
nullDist<-nullDistG1<-nullDistG2<-rep(NA,noPerm)

for(perm in 1:noPerm){

for(switch in 1:number_of_switches){
#select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
#select one of their observations at random for switching
SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
#Find and individual to switch to.
SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
#Check the reverse switch is possible and reselect connection if not
while(sum(Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)!=randomSolver])==0){
  randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
  SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
  SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
}
#Now we need to find a viable reverse switch for other solvers
ReverseSwitch<-rmultinom(1,1,Grp1_observationNet_Perm[SwitchIn,]*(colnames(Grp1_observationNet_Perm)!=randomSolver))==1

#Now perform the switches
Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]-1
Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]+1
Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]+1
Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]-1


#select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
#select one of their observations at random for switching
SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
#Find and individual to switch to.
SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
#Check the reverse switch is possible and reselect connection if not
while(sum(Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)!=randomSolver])==0){
  randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
  SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
  SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
}
#Now we need to find a viable reverse switch for other solvers
ReverseSwitch<-rmultinom(1,1,Grp2_observationNet_Perm[SwitchIn,]*(colnames(Grp2_observationNet_Perm)!=randomSolver))==1

#Now perform the switches
Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]-1
Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]+1
Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]+1
Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]-1
}

#Recalculate the proportion observation
Grp1_proportionObserved_Perm<-Grp1_observationNet_Perm/Grp1_totals
#Drop the columns with NAs for observation matrix
Grp1_proportionObserved_Perm<-Grp1_proportionObserved_Perm[,!is.nan(Grp1_proportionObserved_Perm[1,])]

Grp2_proportionObserved_Perm<-Grp2_observationNet_Perm/Grp2_totals
#Drop the columns with NAs for observation matrix
Grp2_proportionObserved_Perm<-Grp2_proportionObserved_Perm[,!is.nan(Grp2_proportionObserved_Perm[1,])]

nullDist[perm]<-t.test(c(as.vector(Grp1_proportionObserved_Perm),as.vector(Grp2_proportionObserved_Perm))~c(as.vector(Grp1_maternal_kin),Grp2_maternal_kin))$statistic
nullDistG1[perm]<-t.test(as.vector(Grp1_proportionObserved_Perm)~as.vector(Grp1_maternal_kin))$statistic
nullDistG2[perm]<-t.test(as.vector(Grp2_proportionObserved_Perm)~as.vector(Grp2_maternal_kin))$statistic

}

nullDist<-c(nullDist,testStat)
nullDistG1<-c(nullDistG1,testStatGrp1)
nullDistG2<-c(nullDistG2,testStatGrp2)

save(nullDist,nullDistG1,nullDistG2,file="NullDistsKinTests.rData")

load(file="NullDistsKinTests.rData")

jpeg(file = "NullDistKin.jpg",
           width=1000, height=1000)
hist(nullDist)
abline(v=testStat,col=2)
dev.off()


mean(nullDist<=testStat)
#[1] 9.999e-05
mean(nullDistG1<=testStatGrp1)
#[1] 1
mean(nullDistG2<=testStatGrp2)
#[1] 9.999e-05

hist(nullDistG1)
abline(v=testStatGrp1,col=2)

hist(nullDist)
abline(v=testStat,col=2)






#Is there a bias to observe higher/lower ranK?

#Go back to full observation networks with no individuals cut out
Grp1_observationNet<-Grp1_observation_network[,,1,16]
Grp2_observationNet<-Grp2_observation_network[,,2,16]

Grp1_higherRankNet<-matrix(NA,nrow=length(Grp1_indList),ncol=length(Grp1_indList))

for(i in 1:length(Grp1_indList)){
  for(j in 1:length(Grp1_indList)){
    #Since higher ranks are lower numbers, we want a 1 if the row individual (observer) has a higher number than the column individual (performer)
    Grp1_higherRankNet[i,j]<-1*(transRank1[i]>transRank1[j])
  }}

Grp2_higherRankNet<-matrix(NA,nrow=length(Grp2_indList),ncol=length(Grp2_indList))

for(i in 1:length(Grp2_indList)){
  for(j in 1:length(Grp2_indList)){
    
    Grp2_higherRankNet[i,j]<-1*(transRank2[i]>transRank2[j])
  }}


#I will convert the observation network to proportion observed, otherwise individuals who performed many times will have a disproportionate
#effect on the analysis

#Get total manipulations 
Grp1_totals<-matrix(0,length(Grp1_indList),length(Grp1_indList))
Grp1_solverTable<-table(Grp1_Data$solver)
Grp1_solverTable
for(i in 1:length(Grp1_solverTable)) Grp1_totals[,Grp1_indList==names(Grp1_solverTable)[i]]<-Grp1_solverTable[i]
rbind(Grp1_indList,Grp1_totals)
Grp1_proportionObserved<-Grp1_observationNet/Grp1_totals
Grp1_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp1_proportionObserved[1,])
Grp1_higherRankNet<-Grp1_higherRankNet[,!is.nan(Grp1_proportionObserved[1,])]
Grp1_proportionObserved<-Grp1_proportionObserved[,!is.nan(Grp1_proportionObserved[1,])]



Grp2_totals<-matrix(0,length(Grp2_indList),length(Grp2_indList))
Grp2_solverTable<-table(Grp2_Data$solver)
Grp2_solverTable
for(i in 1:length(Grp2_solverTable)) Grp2_totals[,Grp2_indList==names(Grp2_solverTable)[i]]<-Grp2_solverTable[i]
rbind(Grp2_indList,Grp2_totals)
Grp2_proportionObserved<-Grp2_observationNet/Grp2_totals
Grp2_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp2_proportionObserved[1,])
Grp2_higherRankNet<-Grp2_higherRankNet[,!is.nan(Grp2_proportionObserved[1,])]
Grp2_proportionObserved<-Grp2_proportionObserved[,!is.nan(Grp2_proportionObserved[1,])]


t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_higherRankNet),Grp2_higherRankNet))

#Welch Two Sample t-test
#
#data:  c(as.vector(Grp1_proportionObserved), as.vector(Grp2_proportionObserved)) by c(as.vector(Grp1_higherRankNet), Grp2_higherRankNet)
#t = -2.8977, df = 552.14, p-value = 0.003909
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0020539046 -0.0003943089
#sample estimates:
#  mean in group 0 mean in group 1 
#0.0009234316    0.0021475384 

#Bias in the direction of observing higher ranks, but is it significant?

testStat<-t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_higherRankNet),Grp2_higherRankNet))$statistic


#Get the observation net back to no.observation (it is divided by 10)
Grp1_observationNet_Perm<-Grp1_observationNet*10
Grp2_observationNet_Perm<-Grp2_observationNet*10
colnames(Grp1_observationNet_Perm)<-Grp1_indList
colnames(Grp2_observationNet_Perm)<-Grp2_indList


#The following procedure switches a random observation between two individuals for one solver and does the reverse switch for another solver
#This randomises the network whilst keeping the total number of solves and observations by each individual the same
#I then recalculate the proportion observed and re-run the t test, and repeat to get the null distribution

number_of_switches<-10000
noPerm<-10000
nullDist<-rep(NA,noPerm)

for(perm in 1:noPerm){
  
  for(switch in 1:number_of_switches){
    #select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
    randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
    #select one of their observations at random for switching
    SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
    #Find and individual to switch to.
    SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    #Check the reverse switch is possible and reselect connection if not
    while(sum(Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)!=randomSolver])==0){
      randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
      SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
      SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    }
    #Now we need to find a viable reverse switch for other solvers
    ReverseSwitch<-rmultinom(1,1,Grp1_observationNet_Perm[SwitchIn,]*(colnames(Grp1_observationNet_Perm)!=randomSolver))==1
    
    #Now perform the switches
    Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]-1
    Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]+1
    Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]+1
    Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]-1
    
    
    #select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
    randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
    #select one of their observations at random for switching
    SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
    #Find and individual to switch to.
    SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    #Check the reverse switch is possible and reselect connection if not
    while(sum(Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)!=randomSolver])==0){
      randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
      SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
      SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    }
    #Now we need to find a viable reverse switch for other solvers
    ReverseSwitch<-rmultinom(1,1,Grp2_observationNet_Perm[SwitchIn,]*(colnames(Grp2_observationNet_Perm)!=randomSolver))==1
    
    #Now perform the switches
    Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]-1
    Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]+1
    Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]+1
    Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]-1
  }
  
  #Recalculate the proportion observation
  Grp1_proportionObserved_Perm<-Grp1_observationNet_Perm/Grp1_totals
  #Drop the columns with NAs for observation matrix
  Grp1_proportionObserved_Perm<-Grp1_proportionObserved_Perm[,!is.nan(Grp1_proportionObserved_Perm[1,])]
  
  Grp2_proportionObserved_Perm<-Grp2_observationNet_Perm/Grp2_totals
  #Drop the columns with NAs for observation matrix
  Grp2_proportionObserved_Perm<-Grp2_proportionObserved_Perm[,!is.nan(Grp2_proportionObserved_Perm[1,])]
  
  nullDist[perm]<-t.test(c(as.vector(Grp1_proportionObserved_Perm),as.vector(Grp2_proportionObserved_Perm))~c(as.vector(Grp1_higherRankNet),Grp2_higherRankNet))$statistic
}

nullDist<-c(nullDist,testStat)
hist(nullDist)
abline(v=testStat,col=2)

mean(nullDist<=testStat)
#[1] 0.05849415

#Weak evidence that higher ranks are preferentially watched





#Is there a bias to observe older/younger chimpanzees?

#Go back to full observation networks with no individuals cut out
Grp1_observationNet<-Grp1_observation_network[,,1,16]
Grp2_observationNet<-Grp2_observation_network[,,2,16]

Grp1_olderNet<-matrix(NA,nrow=length(Grp1_indList),ncol=length(Grp1_indList))

for(i in 1:length(Grp1_indList)){
  for(j in 1:length(Grp1_indList)){
    Grp1_olderNet[i,j]<-1*(stAge1[i]<stAge1[j])
  }}

Grp2_olderNet<-matrix(NA,nrow=length(Grp2_indList),ncol=length(Grp2_indList))

for(i in 1:length(Grp2_indList)){
  for(j in 1:length(Grp2_indList)){
    
    Grp2_olderNet[i,j]<-1*(stAge2[i]<stAge2[j])
  }}

#I will convert the observation network to proportion observed, otherwise individuals who performed many times will have a disproportionate
#effect on the analysis

#Get total manipulations 
Grp1_totals<-matrix(0,length(Grp1_indList),length(Grp1_indList))
Grp1_solverTable<-table(Grp1_Data$solver)
Grp1_solverTable
for(i in 1:length(Grp1_solverTable)) Grp1_totals[,Grp1_indList==names(Grp1_solverTable)[i]]<-Grp1_solverTable[i]
rbind(Grp1_indList,Grp1_totals)
Grp1_proportionObserved<-Grp1_observationNet/Grp1_totals
Grp1_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp1_proportionObserved[1,])
Grp1_olderNet<-Grp1_olderNet[,!is.nan(Grp1_proportionObserved[1,])]
Grp1_proportionObserved<-Grp1_proportionObserved[,!is.nan(Grp1_proportionObserved[1,])]



Grp2_totals<-matrix(0,length(Grp2_indList),length(Grp2_indList))
Grp2_solverTable<-table(Grp2_Data$solver)
Grp2_solverTable
for(i in 1:length(Grp2_solverTable)) Grp2_totals[,Grp2_indList==names(Grp2_solverTable)[i]]<-Grp2_solverTable[i]
rbind(Grp2_indList,Grp2_totals)
Grp2_proportionObserved<-Grp2_observationNet/Grp2_totals
Grp2_proportionObserved
#Drop the columns with NAs for both kin and observation matrices
!is.nan(Grp2_proportionObserved[1,])
Grp2_olderNet<-Grp2_olderNet[,!is.nan(Grp2_proportionObserved[1,])]
Grp2_proportionObserved<-Grp2_proportionObserved[,!is.nan(Grp2_proportionObserved[1,])]



t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_olderNet),Grp2_olderNet))

#Welch Two Sample t-test
#
#data:  c(as.vector(Grp1_proportionObserved), as.vector(Grp2_proportionObserved)) by c(as.vector(Grp1_olderNet), Grp2_olderNet)
#t = -2.4826, df = 552.65, p-value = 0.01334
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.0018689730 -0.0002178543
#sample estimates:
#  mean in group 0 mean in group 1 
#0.0009921824    0.0020355960

#Bias in the direction of observing older chimpanzees, but is it significant?

testStat<-t.test(c(as.vector(Grp1_proportionObserved),as.vector(Grp2_proportionObserved))~c(as.vector(Grp1_olderNet),Grp2_olderNet))$statistic


#The following procedure switches a random observation between two individuals for one solver and does the reverse switch for another solver
#This randomises the network whilst keeping the total number of solves and observations by each individual the same
#I then recalculate the proportion observed and re-run the t test, and repeat to get the null distribution

number_of_switches<-10000
noPerm<-10000
nullDist<-rep(NA,noPerm)

for(perm in 1:noPerm){
  
  for(switch in 1:number_of_switches){
    #select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
    randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
    #select one of their observations at random for switching
    SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
    #Find and individual to switch to.
    SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    #Check the reverse switch is possible and reselect connection if not
    while(sum(Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)!=randomSolver])==0){
      randomSolver<-names(Grp1_solverTable)[rmultinom(1,1,Grp1_solverTable*(Grp1_solverTable>1))]
      SwitchOut<-rmultinom(1,1,prob=Grp1_observationNet_Perm[,colnames(Grp1_observationNet_Perm)==randomSolver])==1
      SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    }
    #Now we need to find a viable reverse switch for other solvers
    ReverseSwitch<-rmultinom(1,1,Grp1_observationNet_Perm[SwitchIn,]*(colnames(Grp1_observationNet_Perm)!=randomSolver))==1
    
    #Now perform the switches
    Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchOut,colnames(Grp1_observationNet_Perm)==randomSolver]-1
    Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]<-Grp1_observationNet_Perm[SwitchIn,colnames(Grp1_observationNet_Perm)==randomSolver]+1
    Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchOut,ReverseSwitch]+1
    Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp1_observationNet_Perm[SwitchIn,ReverseSwitch]-1
    
    
    #select a solver at random in proprotion to number of observations, but excluding any observed only once since no switch is possible
    randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
    #select one of their observations at random for switching
    SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
    #Find and individual to switch to.
    SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    #Check the reverse switch is possible and reselect connection if not
    while(sum(Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)!=randomSolver])==0){
      randomSolver<-names(Grp2_solverTable)[rmultinom(1,1,Grp2_solverTable*(Grp2_solverTable>1))]
      SwitchOut<-rmultinom(1,1,prob=Grp2_observationNet_Perm[,colnames(Grp2_observationNet_Perm)==randomSolver])==1
      SwitchIn<-rmultinom(1,1,prob=1-SwitchOut)==1
    }
    #Now we need to find a viable reverse switch for other solvers
    ReverseSwitch<-rmultinom(1,1,Grp2_observationNet_Perm[SwitchIn,]*(colnames(Grp2_observationNet_Perm)!=randomSolver))==1
    
    #Now perform the switches
    Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchOut,colnames(Grp2_observationNet_Perm)==randomSolver]-1
    Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]<-Grp2_observationNet_Perm[SwitchIn,colnames(Grp2_observationNet_Perm)==randomSolver]+1
    Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchOut,ReverseSwitch]+1
    Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]<-Grp2_observationNet_Perm[SwitchIn,ReverseSwitch]-1
  }
  
  #Recalculate the proportion observation
  Grp1_proportionObserved_Perm<-Grp1_observationNet_Perm/Grp1_totals
  #Drop the columns with NAs for observation matrix
  Grp1_proportionObserved_Perm<-Grp1_proportionObserved_Perm[,!is.nan(Grp1_proportionObserved_Perm[1,])]
  
  Grp2_proportionObserved_Perm<-Grp2_observationNet_Perm/Grp2_totals
  #Drop the columns with NAs for observation matrix
  Grp2_proportionObserved_Perm<-Grp2_proportionObserved_Perm[,!is.nan(Grp2_proportionObserved_Perm[1,])]
  
  nullDist[perm]<-t.test(c(as.vector(Grp1_proportionObserved_Perm),as.vector(Grp2_proportionObserved_Perm))~c(as.vector(Grp1_olderNet),Grp2_olderNet))$statistic
}

nullDist<-c(nullDist,testStat)
hist(nullDist)
abline(v=testStat,col=2)

mean(nullDist<=testStat)
#[1] 0.1750825

#No real evidence that older individuals are preferentially copied
#You can see in the histogram (see attached png) that there is a bimodal null distribution. This is no a result of any dodgy autocorrelation in the permutations, we can see this if we plot
#the null distribution as a time series
plot(nullDist)
#See attached png

#So we have strong evidence that these chimpanzees preferred to watch maternal kin, and some evidence they preferred to watch higher ranks. However, no evidence that they
#preferred to watch older (or younger) individuals. This is using a data-stream permutation technique that controls for (keeps constant) the number of times each individual observed
#and was observed. In other words I switched around the individual observations in a way that maintained these things, rather than permuting the rows and columns of the matrix.

