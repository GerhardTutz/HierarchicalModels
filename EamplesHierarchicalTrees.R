

##     Examples ordinal hierarchical with Polytree

################################
### Fears data
########################################
########################################

#load("C:/Users/tutz/LRZ Sync+Share/ABookOrdinal/R/Families/GLES17angst.rda")
load("GLES17angst.rda")
summary(GLES)

dat<-GLES
k<-7
nameresp<-"Terrorism"
pred<-c("Age","Gender","Abitur")
preddis<-pred
preda<-pred

##########################################
## Ordinal hierarchical models
######################################

## with cml
fitc<-OrdinalHierTrees(dat,k,nameresp,pred,preda,preddis,fittype="clm",fam="cum")
fitc

fitv<-OrdinalHierTreesr(dat,k,nameresp,pred,preda,preddis,fittype="vglm",fam="cum")
fitv

### Fit directly with Polytrees


listcat<-list(list(c(1,2,3),4,c(5,6,7)),list(1,2,3),list(5,6,7)) 
listpred<-list(pred,preda,preddis)

fit<-PolyTreesGen(dat,k=7,nameresp,listcat,listpred,fittype="nclm",fam="cum")
fit

### adjacent categories model

## with vglm

fitvacat<-OrdinalHierTrees(dat,k,nameresp,pred,preda,preddis,fittype="vglm",fam="acat")
fitvacat


###################################
#### Neutral category separated
#######################################

predn<-pred

### cumulative model in non-neutral category
neuth<-OrdinalHierNeutTrees(dat,k,nameresp,pred,predn,fittype="vglm",hier="yes",fam="cum")
neuth

### adjacent categories model in non-neutral category
neuth2<-OrdinalHierNeutTrees(dat,k,nameresp,pred,predn,fittype="vglm",hier="yes",fam="acat")
neuth2

### hierarchical model in non-neutral category
neutc<-OrdinalHierNeutTrees(dat,k,nameresp,pred,predn,fittype,hier="no",fam="cum")
neutc



######################################
##### happiness even
###################################
###################################

happy <- readRDS("Happiness")

k<-10
nameresp<-"happicat"
pred<-c("age","Gender")
preddis<-pred
preda<-pred


fith<-OrdinalHier(happy,k,nameresp,pred,preda,preddis,fittype="clm")
fith$fit


#### generation happiness data with 10 categories variable happicat

library("CUB")
data(relgoods)
summary(relgoods)

### preparation data
myvars <- c("Gender", "Safety", "BirthYear", "EducationDegree", "WalkAlone","RelNeighbours", "Environment","Happiness")
newdata <- relgoods[myvars]

newdata$age <- 2014 -newdata$BirthYear

summary(newdata)
newdata <- newdata[newdata$age>=18,]
newdata <- newdata[newdata$age<=80,]
summary(newdata)

newdata <- na.omit(newdata) 
names(newdata)<-c("Gender", "Safety","Birthyear","Education","Walk","Neighbours","Environment","Happiness","age")

relgood2 <- newdata
relgood2$happicat<-floor((relgood2$Happiness/112)*10)+1 
relgood2$happicat<-as.factor(relgood2$happicat) 
#hcat<-round((relgood2$Happiness/110)*9,digits=0) 
summary(relgood2)

plot(relgood2$Happiness)
x <- relgood2$Happiness
h<-hist(x, breaks=9, col="blue", xlab="",
        main="")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="red", lwd=2)

summary(relgood2$happicat)

#saveRDS(relgood2, file = "C:\\Users\\tutz\\LRZ Sync+Share\\TuRegLikertHierarchical\\R\\Happiness")

###################  end generation









###########################################
##########  arthritis
######################################
######################################

arthritis <- readRDS("ArthritisUngrouped")

table(arthritis$agent,arthritis$resp)
k<-5
nameresp<-"resp"
pred<-c("agent")
preddis<-pred
preda<-pred

fita<-OrdinalHier(arthritis,k,nameresp,pred,preda,preddis,fittype="clm")
fita
###########################################



