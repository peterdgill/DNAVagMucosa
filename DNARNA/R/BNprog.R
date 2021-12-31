#choose from following:
#typeS<-"Penile swabs"
#typeS<-"Fingernail swabs"
#typeS<-"Boxershorts"
################
##input e.g. BNprog("Boxershorts","Results_DNA_rfu.csv")
BNprog<-function(typeS,datafile){
library(plyr)
library(poibin)
library(rstanarm)
library(ResourceSelection)
library(DescTools)
library(Rcpp)
library(fitdistrplus)
library(bootstrap)
####RUN PREAMBLE BELOW
##Prepare data
rep<-4000#no of bootstraps

#Set variables below
RFUtest<-2##log10(RFU) to test
TimeTest<-0##Time in h to test
###################
#declare dataframes
LRNoDNAPOIonly<-data.frame()## use this node to determine 1-AsT (calc Pr dir T)
LRVagOnly<-data.frame()##Vag mucosa node only - present
LRNoVagOnly<-data.frame()##Vag mucosa node only - not present
LRPOIVM<-data.frame()
LRPOIonly<-data.frame()
LRPOI<-data.frame()
LRVM<-data.frame()
LRnone<-data.frame()
test2<-data.frame()
mydat<-read.csv(datafile)
count=0
T=NULL
R=NULL
if (typeS=="Boxershorts"){
for (RFUtest in seq(from=2, to=6, by=1))
{
T=rbind(T,TimeTest)
R=rbind(R,RFUtest)
Results<-MainProgBox(TimeTest,RFUtest,rep,mydat,typeS)
#####Put results into arrays
LRPOIVM<-rbind.data.frame(LRPOIVM,Results$LRPOIVMr)
LRPOI<-rbind.data.frame(LRPOI,Results$LRPOIr)
LRVM<-rbind.data.frame(LRVM,Results$LRVMr)
LRnone<-rbind.data.frame(LRnone,Results$LRnoner)
LRVagOnly<-rbind.data.frame(LRVagOnly,Results$VagMucOnlyr)
LRNoVagOnly<-rbind.data.frame(LRNoVagOnly,Results$NoVagMucOnlyr)
LRPOIonly<-rbind.data.frame(LRPOIonly,Results$DNAPOIonlyr)
LRNoDNAPOIonly<-rbind.data.frame(LRNoDNAPOIonly,Results$NoDNAPOIonlyr)
}
}else{
for (TimeTest in seq(from=0, to=35,by=5))
{
for (RFUtest in seq(from=2, to=6, by=1))
{
T=rbind(T,TimeTest)
R=rbind(R,RFUtest)
Results<-MainProg(TimeTest,RFUtest,rep,mydat,typeS)
#####Put results into arrays
LRPOIVM<-rbind.data.frame(LRPOIVM,Results$LRPOIVMr)
LRPOI<-rbind.data.frame(LRPOI,Results$LRPOIr)
LRVM<-rbind.data.frame(LRVM,Results$LRVMr)
LRnone<-rbind.data.frame(LRnone,Results$LRnoner)
LRVagOnly<-rbind.data.frame(LRVagOnly,Results$VagMucOnlyr)
LRNoVagOnly<-rbind.data.frame(LRNoVagOnly,Results$NoVagMucOnlyr)
LRPOIonly<-rbind.data.frame(LRPOIonly,Results$DNAPOIonlyr)
LRNoDNAPOIonly<-rbind.data.frame(LRNoDNAPOIonly,Results$NoDNAPOIonlyr)
}}}#end for loops and if condition
#convert to log10
LRPOIVM<-log10(LRPOIVM)
LRPOI<-log10(LRPOI)
LRVM<-log10(LRVM)
LRnone<-log10(LRnone)
LRVagOnly<-log10(LRVagOnly)
LRNoVagOnly<-log10(LRNoVagOnly)
LRPOIonly<-log10(LRPOIonly)
LRNoDNAPOIonly<-log10(LRNoDNAPOIonly)
#####
#Prepare data for printing as .csv files
if (typeS=="Boxershorts"){
R=as.numeric(R)
LRPOIVM<-cbind(R,LRPOIVM)
LRPOI<-cbind(R,LRPOI)
LRVM<-cbind(R,LRVM)
LRnone<-cbind(R,LRnone)
LRVagOnly<-cbind(R,LRVagOnly)
LRNoVagOnly<-cbind(R,LRNoVagOnly)
LRPOIonly<-cbind(R,LRPOIonly)
LRNoDNAPOIonly<-cbind(R,LRNoDNAPOIonly)
colnames(LRNoDNAPOIonly)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOIonly)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRNoVagOnly)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRVagOnly)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOIVM)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOI)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRVM)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRnone)<-c("RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
}else{#run penile or fingernail swabs
T=as.numeric(T)
R=as.numeric(R)
LRPOIVM<-cbind(R,LRPOIVM)
LRPOIVM<-cbind(T,LRPOIVM)

LRPOI<-cbind(R,LRPOI)
LRPOI<-cbind(T,LRPOI)

LRVM<-cbind(R,LRVM)
LRVM<-cbind(T,LRVM)

LRnone<-cbind(R,LRnone)
LRnone<-cbind(T,LRnone)

LRVagOnly<-cbind(R,LRVagOnly)
LRVagOnly<-cbind(T,LRVagOnly)

LRNoVagOnly<-cbind(R,LRNoVagOnly)
LRNoVagOnly<-cbind(T,LRNoVagOnly)

LRPOIonly<-cbind(R,LRPOIonly)
LRPOIonly<-cbind(T,LRPOIonly)

LRNoDNAPOIonly<-cbind(R,LRNoDNAPOIonly)
LRNoDNAPOIonly<-cbind(T,LRNoDNAPOIonly)

colnames(LRNoDNAPOIonly)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOIonly)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRVagOnly)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRNoVagOnly)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOIVM)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRPOI)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRVM)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
colnames(LRnone)<-c("Time","RFU","0.025","0.05","0.25","0.5","0.75","0.95","0.975")
}#end if statement
write.csv(LRNoDNAPOIonly,"Dminus.csv")##DNA node only
write.csv(LRPOIonly,"Dplus.csv")##DNA node only
write.csv(LRNoVagOnly,"Vminus.csv")##Vag node only absence
write.csv(LRVagOnly,"Vplus.csv")## Vag node only presence
write.csv(LRPOIVM,"DplusVplus.csv")## POI +Vag node presence
write.csv(LRPOI,"DplusVminus")# POI + vag node vag- dna+
write.csv(LRVM,"DminusVplus.csv")# POI + Vag node vag+ dna-
write.csv(LRnone,"DminusVminus.csv")# vag- dna-
}#end function
##################################################################### PREAMBLE END

###########################MAIN FUNCTION
MainProg<-function(TimeTest,RFUtest,rep,mydat,typeS){
###Input mydat,rep,TimeTest,RFUtest
PGfile<-subset(mydat,Time_point != "Background" & Location ==typeS)
PGfile$Time<-as.numeric(PGfile$Time_point)
PGfile$POIcontrRfu<-as.numeric(PGfile$Ave_rfu_POI)
PGfile$POIcontrRfu[PGfile$POIcontrRfu==0]<-10##Note change 0RFU to 10RFU so that logs work, otherwise error occurs
####################

########
##Logistic program
#set plot space
Time=PGfile$Time
PGfile$binLR1m<-mapply(binx,log10(PGfile$POIcontrRfu),1)
log10bin=PGfile$binLR1m
#########################################
PGfile$binLR1m<-mapply(binx,log10(PGfile$POIcontrRfu),RFUtest)

log10bin=PGfile$binLR1m#set variables/change according to the input file
#########ADJUST time interval to correspond to direct T
Time=PGfile$Time
dat=as.data.frame(cbind(Time,log10bin))
#Priors can be specified as follows - not used
##priorI=normal(location = PriorIntercept, autoscale = TRUE)##If needed
##priorN=normal(location = PriorTime, autoscale = TRUE)##if needed
set.seed(101)
Model=stan_glm(log10bin~Time, family=binomial(link="logit"),dat)#logistic regression create model
#launch_shinystan(Model)
newdataX<-data.frame(Time=seq(0,36,by=0.2))#generate Time values to test
yweight<-predict(Model,newdata=newdataX,type="resp")#create model
####mydat<-as.data.frame(cbind(newdata=newdataX,yweight))
###SAVE COEFICIENTS IN SIM FILE
##This creates random parameters in var sims
sims <- as.matrix(Model)#generate 4000 random parameters Intercept and Time
colnames(sims)<-c("I","Tc")#change column names and make dataframe
sims<-as.data.frame(sims)#sims contains 2 coefs
##Calculate probs
direct<-mapply(Logistic,sims$I,sims$Tc,TimeTest)#Set to time zero time direct where Time= time since intercourse/offence - time samples collected
AsT<-direct###stored in file AsT

PGfile$binLR1m<-(mapply(binVag,PGfile$Vaginal_mucosa_result))
log10bin=PGfile$binLR1m##binary based upon detected/not detected VM
Time=PGfile$Time
RFU=log10(PGfile$POIcontrRfu)#log10(PGfile$Ave_rfu_POI)
#dat=as.data.frame(cbind(Time,log10bin))
dat=as.data.frame(cbind(RFU,log10bin))
set.seed(101)
Model=stan_glm(log10bin~RFU, family=binomial(link="logit"),dat)##model for quant
newdataX<-data.frame(RFU=seq(0,10,by=0.2))#generate Time values to test
yweight<-predict(Model,newdata=newdataX,type="resp")#create model

##########################
#Put the data into array simVagDir
##This creates random parameters in var sims
simVagDir <- as.matrix(Model)#generate 4000 random parameters Intercept and DNAquant
colnames(simVagDir)<-c("I","Tc")#change column names and make dataframe
simVagDir<-as.data.frame(simVagDir)#sims contains 2 coefs
##Find Probabilities and store in VagDir
VagDir<-mapply(Logistic,simVagDir$I,simVagDir$Tc,RFUtest)#Set to time zero time direct where Time= time since intercourse/offence - time samples collected
###########################
###########################
#Now analyse probabilities for secondary transfer using fitdistrplus and generate 4000 probabilties
###Analyse background DNA from cohabitee
#mydat<-read.csv("Results_DNA_rfu.csv")
PGfile<-subset(mydat,Time_point == "Background" & Location ==typeS)# select secndary T data only
data<-as.numeric(PGfile$Ave_rfu_POI)
set.seed(101)
Trans<-BstrpRFU(data,rep,10^RFUtest)##send to bootstrap turn RFUtest into ordinary number
SecT<-Trans$thetastar##Store the sec transfer probabilities for a given RFUtest
SecVag<-1/23 # SecVag is sec transfer for Vag mucosa - insufficient data to determine - need to plug in a prior here.
bac<-1/23 # assume background VM=0 as there are no data to help
###Prepare LR matrix
LR<-data.frame()
for (no in 1:rep){
Results<-(BNformulaHe(AsT[no],VagDir[no],SecT[no],SecVag,bac))
LR<-rbind.data.frame(LR,Results)
}
###Calculate quantiles for the results
##refers to pos results where labelled
VagMucOnly<-as.numeric(LR$VMPOILR)
VagMucOnly<-replace(VagMucOnly,is.infinite(VagMucOnly),NA)
VagMucOnlyr<-quantile(VagMucOnly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
##### No Vag Muc no DNA test
NoVagMucOnly<-as.numeric(LR$NoVMPOILR)
NoVagMucOnly<-replace(NoVagMucOnly,is.infinite(NoVagMucOnly),NA)
NoVagMucOnlyr<-quantile(NoVagMucOnly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
#####The POIDNA node only no vag muc test
DNAPOIonly<-as.numeric(LR$DNAPOILR)
DNAPOIonly<-replace(DNAPOIonly,is.infinite(DNAPOIonly),NA)
DNAPOIonlyr<-quantile(DNAPOIonly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
########POI DNA NODE ONLY 1-AST node
NoDNAPOIonly<-as.numeric(LR$NoDNAPOILR)
NoDNAPOIonly<-replace(NoDNAPOIonly,is.infinite(NoDNAPOIonly),NA)
NoDNAPOIonlyr<-quantile(NoDNAPOIonly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)


A<-as.numeric(LR$DNAVMLR)
A<-replace(A,is.infinite(A),NA)
B<-as.numeric(LR$DNANOVMLR)
B<-replace(B,is.infinite(B),NA)
C<-as.numeric(LR$VMNODNALR)
C<-replace(C,is.infinite(C),NA)
D<-as.numeric(LR$NOVMNODNALR)
D<-replace(D,is.infinite(D),NA)
LRPOIVMr<-quantile(A,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRPOIr<-quantile(B,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRVMr<-quantile(C,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRnoner<-quantile(D,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
Results<-(list(NoDNAPOIonlyr=NoDNAPOIonlyr,DNAPOIonlyr=DNAPOIonlyr,NoVagMucOnlyr=NoVagMucOnlyr,VagMucOnlyr=VagMucOnlyr,LRPOIVMr=LRPOIVMr,LRPOIr=LRPOIr,LRVMr=LRVMr,LRnoner=LRnoner))
}
###
#############################################
#############################################
###########################################
#FUNCTIONS#################################

BNformulaHe<-function(AsT,VagDir,SecT,SecVag,bac){
DNAPOILR<-(SecT*(1-AsT)+AsT)/SecT
VMPOILR<-((SecVag*(1-VagDir))+VagDir)/(SecVag+(bac*(1-SecVag)))
NoDNAPOILR<-1-AsT
NoVMPOILR<-1-VagDir
##The four possibilities are
#DNA +VM from POI
DNAVMLR<-DNAPOILR*VMPOILR
#DNA and no VM
DNANOVMLR<-DNAPOILR*NoVMPOILR
#VM and no DNA
VMNODNALR<-VMPOILR*NoDNAPOILR
#No Vm and no DNA
NOVMNODNALR<-NoDNAPOILR*NoVMPOILR
Res<-(list(VMPOILR=VMPOILR,DNAVMLR=DNAVMLR,DNANOVMLR=DNANOVMLR,VMNODNALR=VMNODNALR, NOVMNODNALR=NOVMNODNALR,DNAPOILR=DNAPOILR,VMPOILR=VMPOILR,NoDNAPOILR=NoDNAPOILR,NoVMPOILR=NoVMPOILR))
return(Res)
}



##Bin functions
binx <-function(x,Gtrx){if(x< Gtrx){y<-0} else {y<-1}
  return(y)}

binVag<-function(x){if(x=="Not detected"){y<-0} else {y<-1}
  return(y)}
#####################

######LOGISTIC FUNCTION Calculate logistic regression probability FUNCTION
Logistic<-function(I,Tc,H) {
PrH<-1/(1+exp(-(I+Tc*H)))
}
###BOOT FUNCTION
BstrpRFU<-function(data,rep,y){  ##contributor A and B packing only
set.seed(101)#set seed
Bstrp<-bootstrap(data,rep,kcalc3,y)#note that kcalc either sends to LN dist calculator
return(Bstrp)
}
#####################
kcalc3<-function(data,y){#test function with orig data
LNR<-LNCalc(data,y)##return the LN values for all the shedder types
return(LNR)
}
#########################LOGNORMAL CALCULATOR
LNCalc<-function(data,y){
fp=NULL#important to include this in order to initialise the variable
co<-length(data)
ndat<-subset(data,data>0.1) #extract data>0
###ndat<-as.numeric(ndat[,1])
k=co-length(ndat)#count no of data>0 # the number of zero observations
k<-(k/co)#calc proportion of k ie zero RFUs>0.1
##########
tryCatch({
fp<-fitdist(ndat,"lnorm")#fit log normal distr
},
error=function(e){
fp=NULL
})
if(is.null(fp)){
return(NA)
}else{
meanl=fp$estimate[1]#the meanlog coef
sdl=fp$estimate[2]#the sd coef
#set x
LNR<-LNdist(meanl,sdl,k,y)#calc the Pr
#####NOW FEED INTO THE BN
#LR<-MasterBN(0.75,0.75,expR)
return(LNR)
}
}
########################
######################
##the log normal formula###from expAnal prog
LNdist<-function(meanl,sdl,k,y){
invisible(LN<-plnorm(y,meanlog=meanl,sdlog=sdl))
expR<-1-(k+((1-k)*LN))
return(expR)

}


###BOXERSHORTS PROGRAM
MainProgBox<-function(TimeTest,RFUtest,rep,mydat,typeS){
###Input mydat,rep,TimeTest,RFUtest
PGfile<-subset(mydat,Time_point != "Background" & Location ==typeS)
PGfile$Time<-as.numeric(PGfile$Time_point)
PGfile$POIcontrRfu<-as.numeric(PGfile$Ave_rfu_POI)
PGfile$POIcontrRfu[PGfile$POIcontrRfu==0]<-10##Note change 0RFU to 10RFU so that logs work, otherwise error occurs
####################
if (typeS!="Boxershorts") {
########
##Logistic program
#set plot space
Time=PGfile$Time
PGfile$binLR1m<-mapply(binx,log10(PGfile$POIcontrRfu),1)
log10bin=PGfile$binLR1m
#########################################
PGfile$binLR1m<-mapply(binx,log10(PGfile$POIcontrRfu),RFUtest)

log10bin=PGfile$binLR1m#set variables/change according to the input file
#########ADJUST time interval to correspond to direct T
Time=PGfile$Time
dat=as.data.frame(cbind(Time,log10bin))
#Priors can be specified as follows - not used
##priorI=normal(location = PriorIntercept, autoscale = TRUE)##If needed
##priorN=normal(location = PriorTime, autoscale = TRUE)##if needed
set.seed(101)
Model=stan_glm(log10bin~Time, family=binomial(link="logit"),dat)#logistic regression create model
newdataX<-data.frame(Time=seq(0,36,by=0.2))#generate Time values to test
yweight<-predict(Model,newdata=newdataX,type="resp")#create model
####mydat<-as.data.frame(cbind(newdata=newdataX,yweight))
###SAVE COEFICIENTS IN SIM FILE
##This creates random parameters in var sims
sims <- as.matrix(Model)#generate 4000 random parameters Intercept and Time
colnames(sims)<-c("I","Tc")#change column names and make dataframe
sims<-as.data.frame(sims)#sims contains 2 coefs
##Calculate probs
direct<-mapply(Logistic,sims$I,sims$Tc,TimeTest)#Set to time zero time direct where Time= time since intercourse/offence - time samples collected
AsT<-direct###stored in file AsT
}else{###Boxershorts lognormal distribution
data<-as.numeric(PGfile$Ave_rfu_POI)
set.seed(101)
Trans<-BstrpRFU(data,rep,10^RFUtest)##send to bootstrap turn RFUtest into ordinary number
AsT<-Trans$thetastar##Store the direct transfer probabilities for a given RFUtest
}

PGfile$binLR1m<-(mapply(binVag,PGfile$Vaginal_mucosa_result))
log10bin=PGfile$binLR1m##binary based upon detected/not detected VM
Time=PGfile$Time
RFU=log10(PGfile$POIcontrRfu)#log10(PGfile$Ave_rfu_POI)
#dat=as.data.frame(cbind(Time,log10bin))
dat=as.data.frame(cbind(RFU,log10bin))
set.seed(101)
Model=stan_glm(log10bin~RFU, family=binomial(link="logit"),dat)##model for quant
newdataX<-data.frame(RFU=seq(0,10,by=0.2))#generate Time values to test
yweight<-predict(Model,newdata=newdataX,type="resp")#create model

##########################
#Put the data into array simVagDir
##This creates random parameters in var sims
simVagDir <- as.matrix(Model)#generate 4000 random parameters Intercept and DNAquant
colnames(simVagDir)<-c("I","Tc")#change column names and make dataframe
simVagDir<-as.data.frame(simVagDir)#sims contains 2 coefs
##Find Probabilities and store in VagDir
VagDir<-mapply(Logistic,simVagDir$I,simVagDir$Tc,RFUtest)#Set to time zero time direct where Time= time since intercourse/offence - time samples collected
###########################
###########################
#Now analyse probabilities for secondary transfer using fitdistrplus and generate 4000 probabilties
###Analyse background DNA from cohabitee
#mydat<-read.csv("Results_DNA_rfu.csv")
PGfile<-subset(mydat,Time_point == "Background" & Location ==typeS)# select secndary T data only
data<-as.numeric(PGfile$Ave_rfu_POI)
set.seed(101)
Trans<-BstrpRFU(data,rep,10^RFUtest)##send to bootstrap turn RFUtest into ordinary number
SecT<-Trans$thetastar##Store the sec transfer probabilities for a given RFUtest
SecVag<-1/23 # SecVag is sec transfer for Vag mucosa - insufficient data to determine - need to plug in a prior here.
bac<-1/23 # assume background VM=0 as there are no data to help
###Prepare LR matrix
LR<-data.frame()
for (no in 1:rep){
Results<-(BNformulaHe(AsT[no],VagDir[no],SecT[no],SecVag,bac))
LR<-rbind.data.frame(LR,Results)
}

###Calculate quantiles for the results
##refers to pos results where labelled
VagMucOnly<-as.numeric(LR$VMPOILR)
VagMucOnly<-replace(VagMucOnly,is.infinite(VagMucOnly),NA)
VagMucOnlyr<-quantile(VagMucOnly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)

##### No Vag Muc no DNA test
NoVagMucOnly<-as.numeric(LR$NoVMPOILR)
NoVagMucOnly<-replace(NoVagMucOnly,is.infinite(NoVagMucOnly),NA)
NoVagMucOnlyr<-quantile(NoVagMucOnly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
#####The POIDNA node only no vag muc test
DNAPOIonly<-as.numeric(LR$DNAPOILR)
DNAPOIonly<-replace(DNAPOIonly,is.infinite(DNAPOIonly),NA)
DNAPOIonlyr<-quantile(DNAPOIonly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)

########POI DNA NODE ONLY 1-AST node
NoDNAPOIonly<-as.numeric(LR$NoDNAPOILR)
NoDNAPOIonly<-replace(NoDNAPOIonly,is.infinite(NoDNAPOIonly),NA)
NoDNAPOIonlyr<-quantile(NoDNAPOIonly,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)

###Calculate quantiles for the results
##refers to pos results where labelled
A<-as.numeric(LR$DNAVMLR)
A<-replace(A,is.infinite(A),NA)
B<-as.numeric(LR$DNANOVMLR)
B<-replace(B,is.infinite(B),NA)
C<-as.numeric(LR$VMNODNALR)
C<-replace(C,is.infinite(C),NA)
D<-as.numeric(LR$NOVMNODNALR)
D<-replace(D,is.infinite(D),NA)
LRPOIVMr<-quantile(A,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRPOIr<-quantile(B,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRVMr<-quantile(C,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
LRnoner<-quantile(D,probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),na.rm = TRUE)
Results<-(list(NoDNAPOIonlyr=NoDNAPOIonlyr,DNAPOIonlyr=DNAPOIonlyr,NoVagMucOnlyr=NoVagMucOnlyr,VagMucOnlyr=VagMucOnlyr, LRPOIVMr=LRPOIVMr,LRPOIr=LRPOIr,LRVMr=LRVMr,LRnoner=LRnoner))
}
###########
#END
##########
##########
