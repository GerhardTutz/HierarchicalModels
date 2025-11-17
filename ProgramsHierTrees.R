
#### Fits general hierarchical models by using cml and vglm


library("VGAM")
library("effectsize")
library("ordinal")






#######################################################
PolyTreesGen<- function(dat,k,nameresp,listcat,listpred,fittype, fam){

  ### general fitting procedure for cumulative trees with lists 
  #### 
  
  
  ## dat: data frame
  ## k: number of response categories
  ## namesresp: name response variable
  ## listcat: list of partitions of categories
  ## listpred: list of predictors for partitions in listcat
  ## fittype: "clm" fits with clm, otherwise vglm
  ## fam: "cum" for cumulative model or "acat" for adjacent categories model
  
  ### list names
  name_elems <- lapply(listcat, make_name)
  newlist <- as.list(name_elems)
  #####################
  
  dat$resp<-dat[,nameresp]  
  dat$resp<-as.numeric(dat$resp)
  
  ### grouped responses
  
  ## loop comp ### list component
  logliktotal<-0
  AICtotal<-0
  for (comp in 1:length(listcat)){
  dat$respgr<-dat$resp  ## grouped response
  
    lim<-dim(dat)[1]
    length(listcat[[1]])
    for (i in 1:lim) {
      for (ic in 1:length(listcat[[comp]])){
      if (dat$resp[i]%in%unlist(listcat[[comp]][[ic]]))dat$respgr[i]<-ic
    }}
    
    datred<-dat[dat$resp%in%unlist(listcat[[comp]]),]
      
      formula <- as.formula(paste("respgr", "~", paste(unlist(listpred[[comp]]), collapse = " + ")))
      formula
      
    if(fittype!="clm"){
      datred$respgr<-as.ordered(datred$respgr)
      if(fam=="cum")fitgr <- vglm(formula,family=cumulative(parallel=TRUE,reverse=TRUE),data=datred)
      if(fam=="acat")fitgr <- vglm(formula,family=acat(parallel=TRUE,reverse=TRUE),data=datred)
      #fitgra <- vgam(formula,family=cumulative(parallel=TRUE,reverse=TRUE),data=datred)  ### with vgam
      summary(fitgr)
      #summary(fitgra)
            #res<-list(listcat[[comp]],summary(fitgr))
            res<-list(newlist[[comp]],summary(fitgr))
            logliktotal<-logliktotal+logLik(fitgr)
            AICtotal<-AICtotal+AIC(fitgr)
                                        } #nonclm
    if(fittype=="clm"){
      datred$respgr <- as.factor(datred$respgr)
      if(fam=="cum")fitgr <- clm(formula,  data=datred, link = "logit")
      summary(fitgr)
      #res<-list(listcat[[comp]],summary(fitgr))
      res<-list(newlist[[comp]],summary(fitgr))
      logliktotal<-logliktotal+fitgr$logLik
      AICtotal<-AICtotal+AIC(fitgr)
          } #clm
    
    if (comp==1)listfit<-res
    if (comp>1)listfit<-append(listfit,res)
                 
    
      }  ## comp
 
  listfit
  logliktotal 
  AICtotal
  
  
  retlist<-list("loglikhiermodel"=logliktotal,"AIChiermodel"=AICtotal,"listfit"=listfit)  
  return(retlist)
}
########################################
    
#######################################
make_name <- function(sub) {
  if (!is.list(sub)) sub <- as.list(sub)
  parts <- vapply(sub, function(x) paste0(x, collapse = ""), character(1))
  paste(parts, collapse = "/")
}
######################################################


OrdinalHierTrees<- function(dat,k,nameresp,pred,preda,preddis,fittype,fam){
  
  ### generates and fits cumulative trees for k odd and even  
  #### uses clm and vglm
  #### 
  ## dat: data frame
  ## k: number of response categories
  ## namesresp: name response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories
  ##  preddis: names predictors disagreement categories
  ##  fittype: "clm" fits with clm, otherwise vglm
  ## fam: "cum" for cumulative model or "acat" for adjacent categories model
  
  listpred<-list(pred,preda,preddis)
  
  ### k odd 
  if(k/2 !=floor(k/2)){
  listcat<-list(list(seq(1,(k-1)/2,1),(k+1)/2,seq((k+1)/2+1,k,1)),as.list(seq(1,(k-1)/2,1)),as.list(seq((k+1)/2+1,k,1)))} 
  if(k/2 ==floor(k/2)){
    listcat<-list(list(seq(1,k/2,1),seq(k/2+1,k,1)),as.list(seq(1,k/2,1)),as.list(seq(k/2+1,k,1)))} 
  
  fit<-PolyTreesGen(dat,k,nameresp,listcat,listpred,fittype,fam)
  
  retlist<-list("fit"=fit)  
  return(retlist)
  
  
} 
#########################################

#PolyTreesNeut<- function(dat,k,nameresp,pred,predn,fittype,hier)
#OrdinalHierNeutTrees<- function(dat,k,nameresp,pred,predn,fittype,hier)

OrdinalHierNeutTrees<- function(dat,k,nameresp,pred,predn,fittype,hier,fam){
  
  ### fits cumulative trees for k odd and separated neutral category
  #### 
  #### 
  
  ## dat: data frame
  ## k: number of categories must be odd
  ## namesresp: response variable
  ## pred: names predictors within non-neutral categories
  ## predn: names predictors for neutral category or not
  ##  
  ##  fittype: "clm" fits with clm, otherwise vglm
  ##  hier: if "yes" for non-neutral categories a hierarchical model fitted , otherwise cumulative
  ## fam: "cum" for cumulative model or "acat" for adjacent categories model
  
  fit<-0

  if(hier=="yes"){
  #list(list(seq(1,(k-1)/2,1),(k+1)/2,seq((k+1)/2+1,k,1)),as.list(seq(1,(k-1)/2,1)),as.list(seq((k+1)/2+1,k,1)))
  listcat<-list(list(c(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1)), (k+1)/2 ), as.list(c(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1))) )
  listpred<-list(predn,pred)
 
  }

  if(hier!="yes"){
  #listcat<-list(list(c(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1)), (k+1)/2 ), as.list(c(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1))) )
    listcat<-list(list(c(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1)), (k+1)/2 ), list(seq(1,(k-1)/2,1),seq((k+1)/2+1,k,1))  ,as.list(seq(1,(k-1)/2,1)),as.list(seq((k+1)/2+1,k,1)))
    listpred<-list(predn,pred,pred,pred)
    }
  fit<-PolyTreesGen(dat,k,nameresp,listcat,listpred,fittype,fam)
  
  retlist<-list("fit"=fit)  
  return(retlist)
  
  
  } 
################################################
  
  