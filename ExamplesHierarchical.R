
######## Examples fitting of hierarchical models


################################################
#### fears data
###################################################
###################################################

#load("C:/Users/..../LRZ Sync+Share/ABookOrdinal/R/Families/GLES17angst.rda")
#setwd("C:\\Users\\tutz\\LRZ Sync+Share\\TuRegLikertHierarchical\\R\\ProgramaHier")

load("GLES17angst.rda")
summary(GLES)
dat<-GLES

source("ProgramsHierarchical.R")

##############################
#### fit basic hierarchical ordinal
####################################


k<-7
nameresp<-"Terrorism"

pred<-c("Age","Gender","Abitur")
preddis<-pred
preda<-pred

hier<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "LR",preddisp=NULL,equal=NULL)
hier

### fitting results:

round(hier$pargroupedresponse[3:5,] ,digits=3)
round(hier$parnonneutra[5:10,] ,digits=3)

### test results:
equality(hier,tests)
dispersion(hier,tests)

##loglikelihood and AIC
hier$loglikfullmodel
hier$AIC

### alternative with Wald tests

hierw<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "Wald",preddisp=NULL,equal=NULL)
hierw



###################################
#### with dispersion

hierdisp<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests= "no",preddisp= "Age",equal=NULL)
hierdisp
hierdisp$parmdisp
round(hierdisp$parmdisp,digits=3)

#### reduced predictors
hierdisp2<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis="Age",tests = "no",preddisp= "Age",equal=NULL)
hierdisp2




##################################
### fit cumulative model
library("ordinal")

formula_str <- paste(nameresp, "~", paste(pred, collapse = " + "))
formula <- as.formula(formula_str)

dat[,nameresp] <- as.factor(dat[,nameresp])
fitcum <- clm(formula,  data=dat, link = "logit")
summary(fitcum)

datcum<-dat
datcum$Terrorism<- as.ordered(dat$Terrorism)
fit <- vglm(formula, family = cumulative(parallel = TRUE), data = datcum)
fit 
AIC(fit)

fitc <- vglm(formula, family = cumulative(parallel = FALSE), data = datcum)
fitc 
AIC(fitc)

################################ check: comparison with other programs

pr<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "nyes",preddisp=c("Age"))
pr
###  same results:
source("ProgramsHierTrees.R")

flex1<-PolyTreesFlex(dat ,k, nameresp,pred,preda,preddis,fittype="clm")
flex1 

round(pp$parm, digits = 3)
round(pp$testsequality, digits = 3)
round(pp$testsdispersion, digits = 3)




###################################################
##################################################

#### separate neutral category
############################################

dat<-GLES

predn<-pred
prednonn<-pred
preda<-pred
preddis<-pred
sep<-  OrdinalHierNeutral(dat,k,nameresp,predn,prednonn,preda,preddis,tests="Wald",preddisp=NULL,cum="no")  
sep
# binary model that separates neutral:
round(sep$parmsepneutral[2:4,], digits = 3) 
# binary model disagreement versus agreement:
round(sep$fitnonneutral$pargroupedresponse[2:4,], digits = 3)
# model within disagreement and agreement:
round(sep$fitnonneutral$parnonneutral[5:10,], digits = 3)
round(sep$testsequality, digits = 3)
round(pp$testsdispersion, digits = 3)


## with cumulative model in non-neutral categories
sepcum<-  OrdinalHierNeutral(dat,k,nameresp,predn,prednonn,preda,preddis,tests="no",preddisp=NULL,cum="yes")  
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

hier<-OrdinalHierarchical(dat,k,nameresp,pred,preda,preddis,tests = "LR",preddisp=NULL,equal=NULL)
hier
round(hier$pargroupedresponse[2:3,],digits=3)
round(hier$parnonneutral[9:12,],digits=3)

equality(hier,tests)
dispersion(hier,tests)


### reduced model
preda=c("Gender") 
preddis=c("age","Gender")
preddisp=NULL
equal="Gender"

hierred<-OrdinalHierarchical(dat,k,nameresp,pred=c("age","Gender"),preda=c("Gender"),preddis=c("age","Gender"),tests = "no",preddisp=NULL,equal="Gender")
hierred


#### fitting with vgam

fitvglmh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=acat(parallel=TRUE),data=relgood2)
summary(fitvglmh)



