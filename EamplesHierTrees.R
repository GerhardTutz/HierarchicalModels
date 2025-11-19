

##     Examples ordinal hierarchical with Polytree

################################
### Fears data
########################################
########################################


load("GLES17angst.rda")
summary(GLES)


source("ProgramsHierTrees.R")
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

fitv<-OrdinalHierTrees(dat,k,nameresp,pred,preda,preddis,fittype="vglm",fam="cum")
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
neutc<-OrdinalHierNeutTrees(dat,k,nameresp,pred,predn,fittype="vglm",hier="no",fam="cum")
neutc












