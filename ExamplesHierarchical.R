
######## Examples fitting of hierarchical models


################################################
#### fears data
###################################################
###################################################

#load("C:/Users/tutz/LRZ Sync+Share/ABookOrdinal/R/Families/GLES17angst.rda")
load("GLES17angst.rda")
summary(GLES)

dat<-GLES

##############################
#### fit hierarchical ordinal

k<-7
nameresp<-"Terrorism"

pred<-c("Age","Gender","Abitur")
#pred<-c("Age","Gender","EastWest","Abitur","Unemployment")

preddis<-pred
#preddis<-NULL
preda<-pred

hier<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "yes",preddisp=NULL)
hier
round(hier$`parm groupedresponse`[3:5,] ,digits=3)
round(hier$`parm nonneutralcategories`[5:10,] ,digits=3)
round(hier$testsequality ,digits=3)
round(hier$testsdispersion ,digits=3)

###################################
#### with dispersion

hierdisp<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "nyes",preddisp= "Age")
hierdisp
hierdisp$parmdisp


hierdisp2<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis="Age",tests = "nyes",preddisp= "Age")
hierdisp2

#### reduced predictors


##################################
### fit cumulative model
library("ordinal")

formula_str <- paste(nameresp, "~", paste(pred, collapse = " + "))
formula <- as.formula(formula_str)

dat[,nameresp] <- as.factor(dat[,nameresp])
fitcum <- clm(formula,  data=dat, link = "logit")
summary(fitcum)


################################ comparison with other programs

pr<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "nyes",preddisp=c("Age"))
pr
###  same results:
flex1<-PolyTreesFlex(dat ,k, nameresp,pred,preda,preddis,fittype="clm")
flex1 

round(pp$parm, digits = 3)
round(pp$testsequality, digits = 3)
round(pp$testsdispersion, digits = 3)


##################################################
#### separate neutral category

dat<-GLES

predn<-pred
prednonn<-pred
preda<-pred
preddis<-pred
sep<-  OrdinalHierNeutral(dat,k,nameresp,predn,prednonn,preda,preddis,tests="yes",preddisp=NULL,cum="nyes")  
sep
round(sep$parmsepneutral[2:4,], digits = 3)
round(sep$parmnonneutralcategories[5:10,], digits = 3)
round(sep$testsequality, digits = 3)
round(pp$testsdispersion, digits = 3)


## with cumulative model in non-neutral categories
sepcum<-  OrdinalHierNeutral(dat,k,nameresp,predn,prednonn,preda,preddis,tests="yes",preddisp=NULL,cum="yes")  
sepcum


### reduced versions, fewer predictors
predn<-c("Age", "Gender")
prednonn<-pred
preddis<-c("Age")

sepred<-  OrdinalHierNeutral(dat,k,nameresp,predn,prednonn,preda,preddis,tests="nyes",preddisp=NULL,cum="nyes")  
sepred


###################################################
##### happiness data: even number of categories from library("CUB"), data(relgoods)
###################################################
###################################################


dat <- readRDS("Happiness")
k<-10
nameresp<-"happicat"

pred<-c("age","Gender")

preddis<-pred
preda<-pred

hier<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "yes",preddisp=NULL)
hier
round(hier$`parm groupedresponse`[2:3,] ,digits=3)
round(hier$`parm nonneutralcategories`[9:12,] ,digits=3)
round(hier$testsequality ,digits=3)
round(hier$testsdispersion ,digits=3)


#### fitting with vgam

fitvglmh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=acat(parallel=TRUE),data=relgood2)
summary(fitvglmh)



