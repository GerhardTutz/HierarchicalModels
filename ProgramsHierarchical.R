




###################################################################
## Programs hierarchical model fit  
  
#####################

OrdinalHierarchical<- function(dat,k,nameresp,pred,preda,preddis,tests,preddisp, equal){
  
  ## dat: data frame
  ## k: number of response categories
  ## namesresp: response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories, can be NULL
  ##  preddis: names predictors disagreement categories, can be NULL
  ##  equal. names from preda and preddis that have common effects 
  ##  allows for tests 
  ##  predictors can be NULL (not in PolyTreesFlex)
  ##  tests="Wald" computes Wald Tests, only if pred=preda=preddis and  preddisp =NULL 
  
  ## preddisp  variables that are set as dispersion variables 
  ## if not NULL only full model and dispersion model fitted
  
  
  if (length(equal)==0) d<-MakeDesign(dat,k,nameresp,pred,preda,preddis)
  if (length(equal)>0) d<-MakeDesignEqual(dat,k,nameresp,pred,preda,preddis,equal)
  
  ###################################
  ###fits global/grouped
  
  incr<-.2
  knew<-dim(as.matrix(d$respmgr))[2]
  thr<- seq(1,knew-1,1)*incr
  beta<-c(thr,rep(.0,length(pred)))
  d$Xgr<-as.matrix(d$Xgr)
  
  #Loglikcum(beta,d$respmgr,d$ngrgr,d$Xgr)
  
  fitsgr <- optim(beta, Loglikcum, gr = scorecum,d$respmgr,d$ngrgr,d$Xgr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
  fitsgr
  errgr<- solve(fitsgr$hessian)
  stdgr<- sqrt(diag(errgr))
  probgr<-LoglikcumFit(fitsgr$par,d$respmgr,d$ngrgr,d$Xgr)
  
  
  ##output generation
  
  z<-fitsgr$par/stdgr
  pvalgr <-(1-pnorm(abs(z), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
  
  parmgr<- cbind(fitsgr$par,stdgr,z,pvalgr)
  parmgr<-as.data.frame(parmgr)
  
  names(parmgr)[1] <- "Estimates"
  names(parmgr)[2] <- "std err"
  names(parmgr)[3] <- "z-values"
  names(parmgr)[4] <- "p-values"
  if(length(pred)>0)row.names(parmgr)<-c(seq(1,(knew-1),1),paste(pred))
  
  
  ###################### end global
  
  #fits disagree,agree and tests if predisp not NULL 
  knew<-dim(as.matrix(d$respm))[2]
  
  if(length(preddisp) ==0){
    
    
    incr<-.1
    thr<- seq(1,knew-1,1)*incr
    beta<-c(thr,thr,rep(.0,length(preddis)),rep(.0,length(preda)))
    beta<-c(thr,thr,rep(.0,dim(d$X)[2]))
    #LoglikInt(beta,d$respm,d$ngr,d$X,d$indagr)
    #Loglikcum(beta,d$respm,d$ngr,d$X)
    
    fits <- optim(beta, LoglikInt, gr = scoreInt,d$respm,d$ngr,d$X,d$indagr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
    fits
    err<- solve(fits$hessian)
    std<- sqrt(diag(err))
    
    ##output generation
    
    z<-fits$par/std
    pval <-(1-pnorm(abs(z), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    parm<- cbind(fits$par,std,z,pval)
    parm<-as.data.frame(parm)
    
    names(parm)[1] <- "Estimates"
    names(parm)[2] <- "std err"
    names(parm)[3] <- "z-values"
    names(parm)[4] <- "p-values"
    
    
    if(length(equal)==0){
      if(length(preddis)*length(preda)>0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("disagree",preddis),paste("agree",preda))
      if(length(preddis)==0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("agree",preda))
      if(length(preda)==0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("disagree",preddis))
    }
    
    if(length(equal)>0){predwithout<-setdiff(pred, equal)
    predawithout<-setdiff(preda, equal)
    preddiswithout<-setdiff(preddis, equal)
    if(length(predawithout)*length(preddiswithout)>0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("equal",equal),paste("disagree",preddiswithout),paste("agree",predawithout))
    if(length(predawithout)==0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("equal",equal),paste("disagree",preddiswithout))
    if(length(preddiswithout)==0)row.names(parm)<-c(seq(1,2*(knew-1),1),paste("equal",equal),paste("agree",predawithout))
    }
    ### log-likelihood
    
    probnonneutr<-LoglikIntFit(fits$par,d$respm,d$ngr,d$X,d$indagr)
    sum <-0
    
    dat0<-d$datmod  ### extended data
    
    
    ### k odd 
    if(k/2 !=floor(k/2)){
      for (i in 1:dim(dat0)[1]) {respn<-dat0$resp0[i]
      obs<-dat0$obs[i]
      
      if(respn==(k+1)/2)prob<-probgr[i]
      probnow<-0
      if(respn!=(k+1)/2){
        for (s in 1:dim(d$X)[1]){  if(d$obsind[s]==obs) probnow<-probnonneutr[s] }
        prob<-probgr[i]*probnow
      }  
      
      #print (i)
      #print(probnow)
      sum<-sum+log(prob)
      
      }
      loglik<-sum
      dim(d$X)[2]
      #numpar<-length(pred)+length(preda)+length(preddis)+2+2*(k-3)/2
      numpar<-length(pred)+2+ dim(d$X)[2]+(k-3)
      AIC<--2*(loglik-numpar)
    } ### k odd
    
    
    ### k even
    if(k/2 ==floor(k/2)){
      for (i in 1:dim(dat0)[1]) {respn<-dat0$resp0[i]
      obs<-dat0$obs[i]
      
      #if(respn==(k+1)/2)prob<-probgr[i]
      probnow<-0
      #if(respn!=(k+1)/2){
      for (s in 1:dim(d$X)[1]){  if(d$obsind[s]==obs) probnow<-probnonneutr[s] }
      prob<-probgr[i]*probnow
      #}  
      
      #print (i)
      #print(probnow)
      sum<-sum+log(prob)
      
      }
      loglik<-sum
      #numpar<-length(pred)+length(preda)+length(preddis)+1+2*(k/2-1)
      numpar<-length(pred)+1+ dim(d$X)[2]+2*(k/2-1)
      AIC<--2*(loglik-numpar)
    } ### k even
    
    
    
    
    #### Tests
    if(length(equal) >0)tests<-"no"
    
    matrixtests<-0
    matrixdisp<-0
    
    if(tests=="Wald") {
      dbet<-length(beta)
      varnum<-length(pred)
      
      ### equality
      matrixtests<-matrix (0, varnum,2)
      for(v in 1:varnum){
        C<-matrix(0,1,dbet)
        var<-v
        C[1,2*(knew-1)+var]<-1
        C[1,2*(knew-1)+length(pred)+var]<--1
        #w<- C%*%fits$par*(C%*%err%*%t(C))^(-1)*C%*%fits$par
        w <- t(C %*% fits$par) %*% solve(C %*% err %*% t(C)) %*% (C %*% fits$par)
        pvalt<-1-pchisq(w, df=1) 
        matrixtests[v,1]<-w
        matrixtests[v,2]<-pvalt
      }
      matrixtests
      round(matrixtests, digits = 3)
      colnames(matrixtests)  <- c("z-values","p-values")
      row.names(matrixtests)<-c(paste(pred))
      
      ### dispersiontests
      matrixdisp<-matrix (0, varnum,2)
      for(v in 1:varnum){
        C<-matrix(0,1,dbet)
        var<-v
        C[1,2*(knew-1)+var]<-1
        C[1,2*(knew-1)+length(pred)+var]<-1
        w<- C%*%fits$par*(C%*%err%*%t(C))^(-1)*C%*%fits$par
        pvalt<-1-pchisq(w, df=1) 
        matrixdisp[v,1]<-w
        matrixdisp[v,2]<-pvalt
      }
      matrixdisp
      round(matrixdisp, digits = 3)
      colnames(matrixdisp)  <- c("z-values","p-values")
      if(length(preddis)*length(preda)>0)row.names(matrixdisp)<-c(paste(pred))
      if(length(preddis)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("agree",preda))
      if(length(preda)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("disagree",preddis))
    }   ### end Wald
    
    if(tests=="LR") {
      dbet<-length(beta)
      varnum<-length(pred)
      matrixtests<-matrix (0, varnum,2)
      ### equality
      
      for(v in 1:varnum){
        #C<-matrix(0,1,dbet)
        equal<-pred[v]
        d<-MakeDesignEqual(dat,k,nameresp,pred,preda,preddis,equal)
        
        incr<-.1
        thr<- seq(1,knew-1,1)*incr
        numpar<-2*length(thr)+dim(d$X)[2] 
        beta<-c(thr,thr,rep(.0,dim(d$X)[2]))
        fitst <- optim(beta, LoglikInt, gr = scoreInt,d$respm,d$ngr,d$X,d$indagr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = FALSE)
        
        lr<- -2*(-fitst$value+ fits$value)
        pvalt<-1-pchisq(lr, df=1) 
        matrixtests[v,1]<-lr
        matrixtests[v,2]<-pvalt
      }  ### v
      matrixtests
      round(matrixtests, digits = 3)
      colnames(matrixtests)  <- c("z-values","p-values")
      row.names(matrixtests)<-c(paste(pred))
      
      ### dispersiontests
      
      matrixdisp<-matrix (0, varnum,2)
      for(v in 1:varnum){
        #C<-matrix(0,1,dbet)
        disp<-pred[v]
        d<-MakeDesignDisp(dat,k,nameresp,pred,preda,preddis,disp)
        incr<-.1
        thr<- seq(1,knew-1,1)*incr
        numpar<-2*length(thr)+dim(d$X)[2] 
        beta<-c(thr,thr,rep(.0,dim(d$X)[2]))
        fitst <- optim(beta, LoglikInt, gr = scoreInt,d$respm,d$ngr,d$X,d$indagr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = FALSE)
        
        lr<- -2*(-fitst$value+ fits$value)
        pvalt<-1-pchisq(lr, df=1) 
        matrixdisp[v,1]<-lr
        matrixdisp[v,2]<-pvalt
      }  ### v
      matrixdisp
      round(matrixdisp, digits = 3)
      colnames(matrixdisp)  <- c("z-values","p-values")
      if(length(preddis)*length(preda)>0)row.names(matrixdisp)<-c(paste(pred))
      if(length(preddis)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("agree",preda))
      if(length(preda)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("disagree",preddis))
      
    }   ### end LR
    
    
    
    
    
  }  ## end no dispersion
  
  
  #### model with dispersion
  #### ignores equal
  
  if(length(preddisp) >0){
    incr<-.1
    thr<- seq(1,knew-1,1)*incr
    
    ##new
    d<-MakeDesignDisp(dat,k,nameresp,pred,preda,preddis,preddisp)
    d$X
    numpar<-2*length(thr)+dim(d$X)[2] 
    beta<-c(thr,thr,rep(.0,dim(d$X)[2]))
    length(beta)
    #LoglikInt(beta,d$respm,d$ngr,d$X,d$indagr)
    fitsp <- optim(beta, LoglikInt, gr = scoreInt,d$respm,d$ngr,d$X,d$indagr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
    fitsp
    
    err<- solve(fitsp$hessian)
    std<- sqrt(diag(err))
    
    z<-fitsp$par/std
    pval <-(1-pnorm(abs(z), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    ##output generation
    
    #parmdisp<- cbind(fitsp$par,std,z,pval)
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    
    parmdisp<- cbind(fitsp$par,std,z,pval)
    parmdisp<-as.data.frame(parmdisp)
    
    
    if(length(preddis)*length(preda)>0)row.names(parmdisp)<-c(seq(1,2*(knew-1),1),paste("dispersion",preddisp),paste("disagree",predwithout),paste("agree",predawithout))
    if(length(preddis)==0)row.names(parmdisp)<-c(seq(1,2*(knew-1),1),preddisp,paste("agree",preda))
    if(length(preda)==0)row.names(parmdisp)<-c(seq(1,2*(knew-1),1),preddisp,paste("disagree",preddis))
    parmdisp
    parmdisp<-as.data.frame(parmdisp)
    
    names(parmdisp)[1] <- "Estimates"
    names(parmdisp)[2] <- "std err"
    names(parmdisp)[3] <- "z-values"
    names(parmdisp)[4] <- "p-values"
    
    
    
    ### log-likelihood
    
    probnonneutr<-LoglikIntFit(fitsp$par,d$respm,d$ngr,d$X,d$indagr)
    sum <-0
    
    dat0<-d$datmod  ### extended data
    
    ### k odd 
    if(k/2 !=floor(k/2)){
      for (i in 1:dim(dat0)[1]) {respn<-dat0$resp0[i]
      obs<-dat0$obs[i]
      
      if(respn==(k+1)/2)prob<-probgr[i]
      probnow<-0
      if(respn!=(k+1)/2){
        for (s in 1:dim(d$X)[1]){  if(d$obsind[s]==obs) probnow<-probnonneutr[s] }
        prob<-probgr[i]*probnow
      }
      
      #print (i)
      #print(probnow)
      sum<-sum+log(prob)
      
      }  ### i
      
      loglikdisp<-sum
      #numpar<-length(pred)+length(preda)+length(preddis)+2+2*(k-3)/2
      #AICdisp<--2*(loglikdisp-(numpar-length(preddisp)))
      AICdisp<--2*(loglikdisp-numpar)
    }  ### k odd
    
    ### k even
    if(k/2 ==floor(k/2)){
      for (i in 1:dim(dat0)[1]) {respn<-dat0$resp0[i]
      obs<-dat0$obs[i]
      
      #if(respn==(k+1)/2)prob<-probgr[i]
      probnow<-0
      #if(respn!=(k+1)/2){
      for (s in 1:dim(d$X)[1]){  if(d$obsind[s]==obs) probnow<-probnonneutr[s] }
      prob<-probgr[i]*probnow
      #}
      
      #print (i)
      #print(probnow)
      sum<-sum+log(prob)
      
      }
      loglikdisp<-sum
      #numpar<-length(pred)+length(preda)+length(preddis)+1+2*(k/2-1)
      #AICdisp<--2*(loglikdisp-(numpar-length(preddisp)))
      AICdisp<--2*(loglikdisp-numpar)
    }  ### k even
    
    
    
    
    ##  output disp
    #newList <- list("loglikfullmodel"=loglik,"AIC"=AIC,"pargroupedresponse"=parmgr,"loglikgroupedmodel"=-fitsgr$value,"parnonneutral"=parm,"logliknonneutral"=-fits$value,"loglikdispmodel"=loglikdisp,"parmdisp"=parmdisp,
    #                "testsequality"=matrixtests,testsdispersion=matrixdisp,"loglikdispersionmodel"=loglikdisp,"AICdisp"=AICdisp,"parmdispersionmodel"=fitsp$par)#,"fits"=fits, "stderr"=std, "zval"=z )
    newList <- list("loglikdispersionmodel"=loglikdisp,"AICdisp"=AICdisp,"parmdispersionmodel"=parmdisp)#,"fits"=fits, "stderr"=std, "zval"=z )
    
    
  } ## disp
  ############################
  
  ## shortened output
  if(length(preddisp)==0)newList <- list("loglikfullmodel"=loglik,"AIC"=AIC, "pargroupedresponse"=parmgr,"loglikgroupedmodel"=-fitsgr$value,"parnonneutral"=parm,"logliknonneutral"=-fits$value,
                                         "testsequality"=matrixtests,testsdispersion=matrixdisp)#,"fits"=fits, "stderr"=std, "zval"=z )
  
  return(newList) 
  
  
}





#################################

MakeDesign <- function (dat,k,nameresp,pred,preda,preddis){  
  
  ## generates design and response fit disagree, agree categories, model <r and -x 
  
  ## dat: data frame
  ## k: number of categories
  ## namesresp: response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories
  ##  preddis: names predictors disagreement categories
  
  # output  
  ## $datmod is modified data, contains variables resp0 (response), respgr (grouped response, 3 categories if k odd, 2 if k even)
  
  ## respmgr: respmatrix for grouped response
  ## ngrgr: uninteresting, just ones
  ## Xgr:  explanatory variables for groped response
  
  
  ## respm: respmatrix for  response in non-neutral categories (3 categories if k=6 or 7)
  ## ngr: uninteresting
  ## X: explanatory variables for response in non-neutral categories
  ##    order: disagree, agree
  
  dat$resp0<-as.matrix(dat[,nameresp])  
  dat$resp0<-as.numeric(dat$resp0)  
  dat$respgr<-0
  dat$respneut<-1  ## if non-neutral: 1, if neutral: 2  
  
  ### k odd 
  if(k/2 !=floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k-1)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k+1)/2+1) dat$respgr[i]<-3
    if (dat$resp0[i]==(k+1)/2) dat$respgr[i]<-2
    dat$obs[i]=i
    if (dat$resp0[i]==(k+1)/2)dat$respneut[i]<-2}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],3)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k-1)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k+1)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k+1)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    Xdis<-cbind((datdis[,preddis]),0*(datdis[,preda]))
    Xa<-cbind((0*datagr[,preddis]),(datagr[,preda]))
    dim(Xdis)
    dim(Xa)
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    
    ### response matrix
    knew<-(k-1)/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
    } # end k
    
  
  #### k even
  if(k/2 ==floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k)/2+1) dat$respgr[i]<-2
    dat$obs[i]=i}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],2)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    Xdis<-cbind((datdis[,preddis]),0*(datdis[,preda]))
    Xa<-cbind((0*datagr[,preddis]),(datagr[,preda]))
    
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    
    ### response matrix
    knew<-k/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
  
    ### only for output needed
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    #for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
      } # end k
  
  
    newList <- list("respmgr"=respmgr, "ngrgr"= ngrgr,"Xgr"= Xgr,"respm"=respm, "ngr"= ngr,"X"= X,"indagr"=indagr, "datmod"=dat, "datdis"=datdis,"datagr"=datagr,
                    "obsind"=obsind, "respmsepneut"=respmsepneut,"ngsepneut"= ngsepneut,"Xsepneut"=Xsepneut)
    
    
    return(newList)}
###################################
  
  
  Loglikcum <- function (beta,respm,ngr,X){
    
    #### model Y <= r in predictor - x
    #### respm: responses as counts in n x k matrix 
    #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
    
    n<-dim(X)[1]
    k<- dim(respm)[2]
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k)
    ind<- matrix(0,n,k)
    for (i in 1:n){
      etam <- cbind(diag(k-1),-ones%*%X[i,])   ### predictor - x
      #prob[i,1]<- exp(etam[1,]%*%beta)/(1+exp(etam[1,]%*%beta))
      prob[i,1]<- plogis(etam[1,]%*%beta)
      if(k >=3){for (j in 2:(k-1)) {
        #prob[i,j]<-exp(etam[j,]%*%beta)/(1+exp(etam[j,]%*%beta)) -exp(etam[j-1,]%*%beta)/(1+exp(etam[j-1,]%*%beta))}
        prob[i,j]<-plogis(etam[j,]%*%beta) -plogis(etam[j-1,]%*%beta)}
        }
      prob[i,k]<- 1-sum(prob[i,])
      #ind[i,k]<-NegVal(prob[i,k])
    }
    
    ## loglik 
    loglik <-0
    for (i in 1:n)loglik <- loglik +respm[i,]%*%log(prob[i,])
    
    loglik <- -loglik ### negative log-likelihood
    return(loglik)
  }
  #################################
  

  
#########################################  
  LoglikInt <- function (beta,respm,ngr,X,indagr){
    
    #### model Y <= r in predictor - x
    #### respm: responses as counts in n x k matrix 
    #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
    #### X: design explanatory variables
    #### indagr
    
    n<-dim(X)[1]
    k<- dim(respm)[2]
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k)
    ind<- matrix(0,n,k)
    for (i in 1:n){
      etam <- cbind((1-indagr[i])*diag(k-1),indagr[i]*diag(k-1),-ones%*%X[i,])   ### predictor - x
      prob[i,1]<- plogis(etam[1,]%*%beta)
      if(k >=3){for (j in 2:(k-1)) {
        #prob[i,j]<-exp(etam[j,]%*%beta)/(1+exp(etam[j,]%*%beta)) -exp(etam[j-1,]%*%beta)/(1+exp(etam[j-1,]%*%beta))}
        prob[i,j]<-plogis(etam[j,]%*%beta) -plogis(etam[j-1,]%*%beta)}
      }
      prob[i,k]<- 1-sum(prob[i,])
      #ind[i,k]<-NegVal(prob[i,k])
    }
    
    ## loglik 
    loglik <-0
    for (i in 1:n)loglik <- loglik +respm[i,]%*%log(prob[i,])
    
    loglik <- -loglik ### negative log-likelihood
    return(loglik)
  }
#################################
  scorecum <- function (beta,respm,ngr,X){
    
    ### k-1 matrices
    n<-dim(X)[1]
    k<- dim(respm)[2]
    k1<-k-1
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k-1)
    cumprob<- matrix(0,n,k-1)
    dimdum<- length(beta)
    score <- matrix(0,dimdum,1)
    
    for (i in 1:n){
      #etam <- cbind(diag(k-1),-ones%*%X[i,])  ### design matrix X!
      etam <- cbind(diag(k-1),-ones%*%X[i,])   ### predictor - x
      ### probabilities
      prob[i,1]<- plogis(etam[1,]%*%beta)
      cumprob[i,1]<- prob[i,1]
      if(k  >=3){for (j in 2:(k-1)) {
        #prob[i,j]<-exp(etam[j,]%*%beta)/(1+exp(etam[j,]%*%beta)) -        exp(etam[j-1,]%*%beta)/(1+exp(etam[j-1,]%*%beta))
        prob[i,j]<-plogis(etam[j,]%*%beta) -  plogis(etam[j-1,]%*%beta)
        cumprob[i,j]<-plogis(etam[j,]%*%beta)}}
      
      #prob[i,k]<- 1-sum(prob[i,])
      
      ### covariance
      #vect<-matrix(prob[i,],k-1,1)
      if(k1 >=2)Sigma <- (diag(prob[i,])-prob[i,]%*%t(prob[i,]))/ngr[i]
      if(k1 <2) Sigma <- ((prob[i,])-prob[i,]*(prob[i,]))/ngr[i]
      ### Derivative 
      D<- matrix(0,k-1,k-1)
      for(ii in 1:k1) {
        for(r in 1:k1){if(ii <= r)D[ii,r]<- 1/((1-cumprob[i,r])*cumprob[i,r])
        }}
      D <- solve(D)
      #det(D)
      Sigmainv<-solve(Sigma)
      
      respshort<- respm[i,1:k1]/ngr[i]  ### relative frequencies
      score <- score +t(etam)%*%D%*%Sigmainv%*%(respshort-prob[i,])
    }## end i
    
    ## loglik 
    score <- -score ### negative log-likelihood
    return(score)
  }
  
  
  
  
#################################  
scoreInt <- function (beta,respm,ngr,X,indagr){
  
  ### k-1 matrices
  n<-dim(X)[1]
  k<- dim(respm)[2]
  k1<-k-1
  ones<- matrix(1,k-1,1)
  prob<- matrix(0,n,k-1)
  cumprob<- matrix(0,n,k-1)
  dimdum<- length(beta)
  score <- matrix(0,dimdum,1)
  
  for (i in 1:n){
    #etam <- cbind(diag(k-1),-ones%*%X[i,])  ### design matrix X!
    etam <- cbind((1-indagr[i])*diag(k-1),indagr[i]*diag(k-1),-ones%*%X[i,])   ### predictor - x
    ### probabilities
    prob[i,1]<- plogis(etam[1,]%*%beta)
    cumprob[i,1]<- prob[i,1]
    if(k  >=3){
      for (j in 2:(k-1)) {prob[i,j]<-plogis(etam[j,]%*%beta) -
      plogis(etam[j-1,]%*%beta)
      cumprob[i,j]<-plogis(etam[j,]%*%beta)}
              }
    
    #cumprob[i,]
    #prob[i,]
    #prob[i,]
    #prob[i,k]<- 1-sum(prob[i,])
    
    ### covariance
    #vect<-matrix(prob[i,],k-1,1)
    if(k1 >=2)Sigma <- (diag(prob[i,])-prob[i,]%*%t(prob[i,]))/ngr[i]
    if(k1 <2) Sigma <- ((prob[i,])-prob[i,]*(prob[i,]))/ngr[i]
    ### Derivative 
    D<- matrix(0,k-1,k-1)
    for(ii in 1:k1) {
      for(r in 1:k1){if(ii <= r)D[ii,r]<- 1/((1-cumprob[i,r])*cumprob[i,r])
      }}
    D <- solve(D)
    #det(D)
    Sigmainv<-solve(Sigma)
    
    respshort<- respm[i,1:k1]/ngr[i]  ### relative frequencies
    score <- score +t(etam)%*%D%*%Sigmainv%*%(respshort-prob[i,])
  }## end i
  
  ## loglik 
  score <- -score ### negative log-likelihood
  return(score)
}


##############################################

LoglikIntFit <- function (beta,respm,ngr,X,indagr){
  
  #### returns fitted probabilities
  #### model Y <= r in predictor - x
  #### respm: responses as counts in n x k matrix 
  #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
  
  #returns fitted probabilities
  n<-dim(X)[1]
  k<- dim(respm)[2]
  ones<- matrix(1,k-1,1)
  prob<- matrix(0,n,k)
  probobs<- matrix(0,n,1) #probabilities observed responses
  
  ind<- matrix(0,n,k)
  for (i in 1:n){
    etam <- cbind((1-indagr[i])*diag(k-1),indagr[i]*diag(k-1),-ones%*%X[i,])   ### predictor - x
    prob[i,1]<- plogis(etam[1,]%*%beta)
    if(k >=3){for (j in 2:(k-1)) {prob[i,j]<-plogis(etam[j,]%*%beta) -
      plogis(etam[j-1,]%*%beta)}
      #ind[i,j]<-NegVal(prob[i,j])
    }
    prob[i,k]<- 1-sum(prob[i,])
    for (ii in 1:k)if (respm[i,ii]==1)probobs[i,1]<-prob[i,ii]
    #ind[i,k]<-NegVal(prob[i,k])
  }
  
  ## loglik 
  #loglik <-0
  #for (i in 1:n)loglik <- loglik +respm[i,]%*%log(prob[i,])
  
  #loglik <- -loglik ### negative log-likelihood
  
  return(probobs) ## returns prob observed
}
#################################


  LoglikIntPen <- function (beta,respm,ngr,X,indagr,C, lambda){
    
    ### penalizes (C*beta)^2 with weight lamda
    #### model Y <= r in predictor - x
    #### respm: responses as counts in n x k matrix 
    #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
    
    n<-dim(X)[1]
    k<- dim(respm)[2]
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k)
    ind<- matrix(0,n,k)
    for (i in 1:n){
      etam <- cbind((1-indagr[i])*diag(k-1),indagr[i]*diag(k-1),-ones%*%X[i,])   ### predictor - x
      prob[i,1]<- plogis(etam[1,]%*%beta)
      if(k >=3){for (j in 2:(k-1)) {prob[i,j]<-plogis(etam[j,]%*%beta) -
        plogis(etam[j-1,]%*%beta)}
        #ind[i,j]<-NegVal(prob[i,j])
      }
      prob[i,k]<- 1-sum(prob[i,])
      #ind[i,k]<-NegVal(prob[i,k])
    }
    
    ## loglik 
    loglik <-0
    for (i in 1:n)loglik <- loglik +respm[i,]%*%log(prob[i,])
    loglik<-loglik-lambda*t(C%*%as.matrix(beta))%*%C%*%as.matrix(beta)
    
    loglik <- -loglik ### negative log-likelihood
    return(loglik)
  }

  ######################################
  scoreIntPen <- function (beta,respm,ngr,X,indagr,C, lambda){
    
    ### k-1 matrices
    n<-dim(X)[1]
    k<- dim(respm)[2]
    k1<-k-1
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k-1)
    cumprob<- matrix(0,n,k-1)
    dimdum<- length(beta)
    score <- matrix(0,dimdum,1)
    
    for (i in 1:n){
      #etam <- cbind(diag(k-1),-ones%*%X[i,])  ### design matrix X!
      etam <- cbind((1-indagr[i])*diag(k-1),indagr[i]*diag(k-1),-ones%*%X[i,])   ### predictor - x
      ### probabilities
      prob[i,1]<- plogis(etam[1,]%*%beta)
      cumprob[i,1]<- prob[i,1]
      if(k >=3){for (j in 2:(k-1)) {prob[i,j]<-plogis(etam[j,]%*%beta) -
        plogis(etam[j-1,]%*%beta)
        cumprob[i,j]<-plogis(etam[j,]%*%beta)}}
      
      #prob[i,k]<- 1-sum(prob[i,])
      
      ### covariance
      #vect<-matrix(prob[i,],k-1,1)
      if(k1 >=2)Sigma <- (diag(prob[i,])-prob[i,]%*%t(prob[i,]))/ngr[i]
      if(k1 <2) Sigma <- ((prob[i,])-prob[i,]*(prob[i,]))/ngr[i]
      ### Derivative 
      D<- matrix(0,k-1,k-1)
      for(ii in 1:k1) {
        for(r in 1:k1){if(ii <= r)D[ii,r]<- 1/((1-cumprob[i,r])*cumprob[i,r])
        }}
      D <- solve(D)
      #det(D)
      Sigmainv<-solve(Sigma)
      
      respshort<- respm[i,1:k1]/ngr[i]  ### relative frequencies
      score <- score +t(etam)%*%D%*%Sigmainv%*%(respshort-prob[i,])
    }## end i
    
    ## loglik 
    #score<-score-2*lambda*as.numeric(C%*%as.matrix(beta))*t(C)
    score<-score- lambda* 2* t(C)%*%C%*%beta  #as.numeric(C%*%as.matrix(beta))*t(C)
    score <- -score ### negative log-likelihood
    return(score)
  }
  
#######################################
  LoglikcumFit <- function (beta,respm,ngr,X){
    
    #### returns fitted probabilities
    #### model Y <= r in predictor - x
    #### respm: responses as counts in n x k matrix 
    #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
    
    #returns fitted probabilities
    n<-dim(X)[1]
    k<- dim(respm)[2]
    ones<- matrix(1,k-1,1)
    prob<- matrix(0,n,k)
    probobs<- matrix(0,n,1) #probabilities observed responses
    
    ind<- matrix(0,n,k)
    for (i in 1:n){
      etam <- cbind(diag(k-1),-ones%*%X[i,])   ### predictor - x
      prob[i,1]<- plogis(etam[1,]%*%beta)
      if(k >=3){for (j in 2:(k-1)) {prob[i,j]<-plogis(etam[j,]%*%beta) -
        plogis(etam[j-1,]%*%beta)}
        #ind[i,j]<-NegVal(prob[i,j])
      }
      prob[i,k]<- 1-sum(prob[i,])
      for (ii in 1:k)if (respm[i,ii]==1)probobs[i,1]<-prob[i,ii]
      #ind[i,k]<-NegVal(prob[i,k])
    }
    
    ## loglik 
    #loglik <-0
    #for (i in 1:n)loglik <- loglik +respm[i,]%*%log(prob[i,])
    
    #loglik <- -loglik ### negative log-likelihood
    
    return(probobs) ## returns prob observed
  } 
  ###################################################



##################################

OrdinalHierNeutral<- function(dat,k,nameresp,predn,prednonn,preda,preddis,tests,preddisp, cum){
  
  ## dat: data frame
  ## k: number of categories
  ## namesresp: response variable
  ##  predn: names predictors first node, neutral versus non-neutral
  ##  prednonn: predictors agreemaent non-agreement
  ##  preda: names predictors agreement categories, can be NULL
  ##  preddis: names predictors disagreement categories, can be NULL
  
  ##  allows for tests  
  ##  predictors can be NULL (not in PolyTreesFlex)
  ##  tests="Wald" only if pred=preda=preddis, only Wald test available
  
  ## preddisp  variables that are set as dispersion variables 
  ## if not NULL only full model and dispersion model fitted
  
  ## cum: if cum ="yes" in second step cumulative model fitted
  ##      with predn and prednonn
  
  d<-MakeDesign(dat,k,nameresp,predn,preda,preddis)
  d  
  
  ## fit neutral
  
  ###################################
  ### fits neutral versus non-neutral
  
  incr<-.2
  knew<-dim(as.matrix(d$respmsepneut))[2]
  thr<- seq(1,knew-1,1)*incr
  beta<-c(thr,rep(.0,length(predn)))
  d$Xsepneut<-as.matrix(d$Xsepneut)
  
  #Loglikcum(beta,d$respmgr,d$ngrgr,d$Xgr)
  fitsgr <- optim(beta, Loglikcum, gr = scorecum,d$respmsepneut,d$ngsepneut,d$Xsepneut,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
  fitsgr
  errgr<- solve(fitsgr$hessian)
  stdgr<- sqrt(diag(errgr))
  probgr<-LoglikcumFit(fitsgr$par,d$respmsepneut,d$Xsepneut,d$Xsepneut)
  #plot(d$datmod$respneut,probgr)
  #plot(probgr,d$datmod$respneut)
  #table(d$datmod$respneut)
  
  ##output generation
  
  z<-fitsgr$par/stdgr
  pvalgr <-(1-pnorm(abs(z), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
  
  parmgr<- cbind(fitsgr$par,stdgr,z,pvalgr)
  parmgr<-as.data.frame(parmgr)
  
  names(parmgr)[1] <- "Estimates"
  names(parmgr)[2] <- "std err"
  names(parmgr)[3] <- "z-values"
  names(parmgr)[4] <- "p-values"
  if(length(predn)>0)row.names(parmgr)<-c(seq(1,(knew-1),1),paste(predn))
  
  
  ###################### end neutral
  
  ################fits disagree,agree hierarchical model
  
  ## non neutral data set: datred
  dat$resp0<-as.matrix(dat[,nameresp])  
  dat$resp0<-as.numeric(dat$resp0)
  
  dat$red<-dat$resp0
  for(i in 1:dim(dat)[1]){if(dat$resp[i] ==((k+1)/2)) {dat$red[i] <- NA}}
  for(i in 1:dim(dat)[1]){if(dat$resp[i] >=((k+3)/2)) {dat$red[i] <- dat$red[i]-1}}
  
  datred <- na.omit(dat) 
  
  if(cum !="yes"){
  
  
  hierred<-OrdinalHierarchical(datred,k-1,nameresp="red",prednonn,preda,preddis,tests ="no",preddisp=NULL,equal=NULL)
  if(all(prednonn==preda) & all(prednonn==preddis) )hierred<-OrdinalHierarchical(datred,k-1,nameresp="red",prednonn,preda,preddis,tests = tests,preddisp=NULL,equal=NULL)
  
  
  knew<-dim(as.matrix(d$respm))[2]
  incr<-.1
  thr<- seq(1,knew-1,1)*incr
  beta<-c(thr,thr,rep(.0,length(preddis)),rep(.0,length(preda)))
  #LoglikInt(beta,d$respm,d$ngr,d$X,d$indagr)
  
  fits <- optim(beta, LoglikInt, gr = scoreInt,d$respm,d$ngr,d$X,d$indagr,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
  fits
  err <- tryCatch(solve(fits$hessian), error = function(e) MASS::ginv(fits$hessian))
  std<- sqrt(diag(err))
  
  
  ### log-likelihood  
  
  loglik<--fitsgr$value+hierred$loglikfullmodel
  numpar<-1+length(predn)+1+length(prednonn)+length(preda)+length(preddis)+2*(k-3)/2
  AIC<--2*(loglik-numpar)
  
  
  #### Tests
  
  matrixtests<-0
  matrixdisp<-0
  if(tests=="Wald") {
    dbet<-length(beta)
    varnum<-length(pred)
    
    ### equality
    matrixtests<-matrix (0, varnum,2)
    for(v in 1:varnum){
      C<-matrix(0,1,dbet)
      var<-v
      C[1,2*(knew-1)+var]<-1
      C[1,2*(knew-1)+length(pred)+var]<--1
      w<- C%*%fits$par*(C%*%err%*%t(C))^(-1)*C%*%fits$par
      pvalt<-1-pchisq(w, df=1) 
      matrixtests[v,1]<-w
      matrixtests[v,2]<-pvalt
    }
    matrixtests
    round(matrixtests, digits = 3)
    colnames(matrixtests)  <- c("z-values","p-values")
    row.names(matrixtests)<-c(paste(pred))
    
    ### dispersiontests
    matrixdisp<-matrix (0, varnum,2)
    for(v in 1:varnum){
      C<-matrix(0,1,dbet)
      var<-v
      C[1,2*(knew-1)+var]<-1
      C[1,2*(knew-1)+length(pred)+var]<-1
      w<- C%*%fits$par*(C%*%err%*%t(C))^(-1)*C%*%fits$par
      pvalt<-1-pchisq(w, df=1) 
      matrixdisp[v,1]<-w
      matrixdisp[v,2]<-pvalt
    }
    matrixdisp
    round(matrixdisp, digits = 3)
    colnames(matrixdisp)  <- c("z-values","p-values")
    if(length(preddis)*length(preda)>0)row.names(matrixdisp)<-c(paste(pred))
    if(length(preddis)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("agree",preda))
    if(length(preda)==0)row.names(matrixdisp)<-c(seq(1,2*(knew-1),1),paste("disagree",preddis))
  }  ## tests
  
  if(length(preddisp)==0)newList <- list("loglikfullmodel"=loglik,"AIC"=AIC, "parmsepneutral"=parmgr,"loglikgroupedmodel"=-fitsgr$value,
                                         "fitnonneutral"=hierred,"numpar"=numpar)
  
  
  
  }  ## cum nonyes 
  ################fits disagree,agree cumulative model
  
  if(cum =="yes"){
    
    Xnonn<-as.matrix(datred[,prednonn])  ### !!!!!!
    respmnonn<- matrix(0,dim(Xnonn)[1],k-1)
    for (i in 1:dim(Xnonn)[1]) respmnonn[i,datred$red[i]]<-1
    ngrnonn <- rep(1,dim(Xnonn)[1])
    
    
    incr<-.2
    knew<-dim(as.matrix(respmnonn))[2]
    thr<- seq(1,knew-1,1)*incr
    beta<-c(thr,rep(.0,length(prednonn)))
    Xnonn<-as.matrix(Xnonn)
    
    #Loglikcum(beta,respmnonn,ngrnonn,Xnonn)
    fitsnonn <- optim(beta, Loglikcum, gr = scorecum,respmnonn,ngrnonn,Xnonn,method = "BFGS",lower = -Inf, upper = Inf,control = list(), hessian = TRUE)
    fitsnonn
    errnonn<- solve(fitsnonn$hessian)
    stdnonn<- sqrt(diag(errnonn))
    #probgr<-LoglikcumFit(fitsgr$par,d$respmgr,d$ngrgr,d$Xgr)
    
    
    ##output generation
    
    z<-fitsnonn$par/stdnonn
    pvalgr <-(1-pnorm(abs(z), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    parmnonn<- cbind(fitsnonn$par,stdgr,z,pvalgr)
    parmnonn<-as.data.frame(parmnonn)
    
    names(parmnonn)[1] <- "Estimates"
    names(parmnonn)[2] <- "std err"
    names(parmnonn)[3] <- "z-values"
    names(parmnonn)[4] <- "p-values"
    if(length(prednonn)>0)row.names(parmnonn)<-c(seq(1,(knew-1),1),paste(prednonn))
    
    ### log-likelihood  
    
    logliknonn<--fitsgr$value-fitsnonn$value
    numpar<-1+length(predn)+knew-1+length(prednonn)
    AICnonn<--2*(logliknonn-numpar)
    
    newList <- list("loglikfullmodel"=logliknonn,"AIC"=AICnonn, "parmsepneutral"=parmgr,"loglikgroupedmodel"=-fitsgr$value,"fitnonneutral"=parmnonn,"numpar"=numpar)
    
    
  } # end cum yes
  
    ############################
  
  return(newList) 
  
  
}
##################################


MakeDesignDisp <- function(dat,k,nameresp,pred,preda,preddis,preddisp){  
  
  ## generates design and response fit disagree, agree categories, model <r and -x 
  
  ## dat: data frame
  ## k: number of categories
  ## namesresp: response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories
  ##  preddis: names predictors disagreement categories
  
  # output  
  ## $datmod is modified data, contains variables resp0 (response), respgr (grouped response, 3 categories if k odd, 2 if k even)
  
  ## respmgr: respmatrix for grouped response
  ## ngrgr: uninteresting, just ones
  ## Xgr:  explanatory variables for groped response
  
  
  ## respm: respmatrix for  response in non-neutral categories (3 categories if k=6 or 7)
  ## ngr: uninteresting
  ## X: explanatory variables for response in non-neutral categories
  #     with dispersion effects
  ##    order: disagree, agree
  
  dat$resp0<-as.matrix(dat[,nameresp])  
  dat$resp0<-as.numeric(dat$resp0)  
  dat$respgr<-0
  dat$respneut<-1  ## if non-neutral: 1, if neutral: 2  
  
  ### k odd 
  if(k/2 !=floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k-1)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k+1)/2+1) dat$respgr[i]<-3
    if (dat$resp0[i]==(k+1)/2) dat$respgr[i]<-2
    dat$obs[i]=i
    if (dat$resp0[i]==(k+1)/2)dat$respneut[i]<-2}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],3)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k-1)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k+1)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k+1)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    
    ## design
    
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    
    Xdis<-cbind((datdis[,preddisp]),datdis[,predwithout],0*(datdis[,predawithout]))
    Xa<-cbind(-datagr[,preddisp],0*datagr[,predwithout],(datagr[,predawithout]))
    dim(Xdis)
    dim(Xa)
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    names(X)
    ### response matrix
    knew<-(k-1)/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
  } # end k
  
  
  #### k even
  if(k/2 ==floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k)/2+1) dat$respgr[i]<-2
    dat$obs[i]=i}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],2)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    ## design
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    
    Xdis<-cbind(datdis[,preddisp],datdis[,predwithout],0*(datdis[,predawithout]))
    Xa<-cbind(-datdis[,preddisp],0*datagr[,predwithout],(datagr[,predawithout]))
    
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    
    ### response matrix
    knew<-k/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
    
    ### only for output needed
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    #for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
  } # end k
  
  
  newList <- list("respmgr"=respmgr, "ngrgr"= ngrgr,"Xgr"= Xgr,"respm"=respm, "ngr"= ngr,"X"= X,"indagr"=indagr, "datmod"=dat, "datdis"=datdis,"datagr"=datagr,
                  "obsind"=obsind, "respmsepneut"=respmsepneut,"ngsepneut"= ngsepneut,"Xsepneut"=Xsepneut)
  
  
  return(newList)}
###################################

MakeDesignEqualOld <- function(dat,k,nameresp,pred,preda,equal){  
  
  ## generates design and response fit disagree, agree categories, model <r and -x 
  
  ## dat: data frame
  ## k: number of categories
  ## namesresp: response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories
  ##  preddis: names predictors disagreement categories
  ##  equal: variable that has same effect in disagreement agreement categories
  # output  
  ## $datmod is modified data, contains variables resp0 (response), respgr (grouped response, 3 categories if k odd, 2 if k even)
  
  ## respmgr: respmatrix for grouped response
  ## ngrgr: uninteresting, just ones
  ## Xgr:  explanatory variables for groped response
  
  
  ## respm: respmatrix for  response in non-neutral categories (3 categories if k=6 or 7)
  ## ngr: uninteresting
  ## X: explanatory variables for response in non-neutral categories
  #     with dispersion effects
  ##    order: disagree, agree
  
  preddisp<-equal
  dat$resp0<-as.matrix(dat[,nameresp])  
  dat$resp0<-as.numeric(dat$resp0)  
  dat$respgr<-0
  dat$respneut<-1  ## if non-neutral: 1, if neutral: 2  
  
  ### k odd 
  if(k/2 !=floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k-1)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k+1)/2+1) dat$respgr[i]<-3
    if (dat$resp0[i]==(k+1)/2) dat$respgr[i]<-2
    dat$obs[i]=i
    if (dat$resp0[i]==(k+1)/2)dat$respneut[i]<-2}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],3)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k-1)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k+1)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k+1)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    
    ## design
    
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    
    Xdis<-cbind((datdis[,preddisp]),datdis[,predwithout],0*(datdis[,predawithout]))
    Xa<-cbind(datagr[,preddisp],0*datagr[,predwithout],(datagr[,predawithout]))
    dim(Xdis)
    dim(Xa)
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    names(X)
    ### response matrix
    knew<-(k-1)/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
  } # end k
  
  
  #### k even
  if(k/2 ==floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k)/2+1) dat$respgr[i]<-2
    dat$obs[i]=i}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],2)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    ## design
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    
    Xdis<-cbind(datdis[,preddisp],datdis[,predwithout],0*(datdis[,predawithout]))
    Xa<-cbind(datdis[,preddisp],0*datagr[,predwithout],(datagr[,predawithout]))
    
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    
    ### response matrix
    knew<-k/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
    
    ### only for output needed
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    #for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
  } # end k
  
  
  newList <- list("respmgr"=respmgr, "ngrgr"= ngrgr,"Xgr"= Xgr,"respm"=respm, "ngr"= ngr,"X"= X,"indagr"=indagr, "datmod"=dat, "datdis"=datdis,"datagr"=datagr,
                  "obsind"=obsind, "respmsepneut"=respmsepneut,"ngsepneut"= ngsepneut,"Xsepneut"=Xsepneut)
  
  
  return(newList)}
###################################

equality<- function (hier,tests){
  
  if(tests == "LR") print ("Likelihood ratio tests for equality of effects within disagreement and agreement effects")
  if(tests == "Wald") print ("Wald tests for equality of effects within disagreement and agreement effects")
  print(round(hier$testsequality ,digits=3))
}

dispersion<- function (hier,tests){
  if(tests == "LR") print ("Likelihood ratio tests if variables have dispersion effect")
  if(tests == "Wald") print ("Wald tests if variables have dispersion effect")
  print(round(hier$testsdispersion ,digits=3))
}


#############################################################
MakeDesignEqual <- function(dat,k,nameresp,pred,preda,preddis,equal){  
  
  ## generates design and response fit disagree, agree categories, model <r and -x 
  
  ## dat: data frame
  ## k: number of categories
  ## namesresp: response variable
  ##  pred: names predictors grouped response
  ##  preda: names predictors agreement categories
  ##  preddis: names predictors disagreement categories
  ##  equal: variable that has same effect in disagreement agreement categories
  # output  
  ## $datmod is modified data, contains variables resp0 (response), respgr (grouped response, 3 categories if k odd, 2 if k even)
  
  ## respmgr: respmatrix for grouped response
  ## ngrgr: uninteresting, just ones
  ## Xgr:  explanatory variables for groped response
  
  
  ## respm: respmatrix for  response in non-neutral categories (3 categories if k=6 or 7)
  ## ngr: uninteresting
  ## X: explanatory variables for response in non-neutral categories
  #     with dispersion effects
  ##    order: disagree, agree
  
  preddisp<-equal
  dat$resp0<-as.matrix(dat[,nameresp])  
  dat$resp0<-as.numeric(dat$resp0)  
  dat$respgr<-0
  dat$respneut<-1  ## if non-neutral: 1, if neutral: 2  
  
  ### k odd 
  if(k/2 !=floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k-1)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k+1)/2+1) dat$respgr[i]<-3
    if (dat$resp0[i]==(k+1)/2) dat$respgr[i]<-2
    dat$obs[i]=i
    if (dat$resp0[i]==(k+1)/2)dat$respneut[i]<-2}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],3)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k-1)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k+1)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k+1)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    
    ## design
    
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    preddiswithout<-setdiff(preddis, preddisp)
    
    Xdis<-cbind(datdis[,preddisp],datdis[,predawithout],0*(datdis[,preddiswithout]))
    Xa<-cbind(datagr[,preddisp],0*datagr[,predawithout],(datagr[,preddiswithout]))
    
    dim(Xdis)
    dim(Xa)
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    names(X)
    ### response matrix
    knew<-(k-1)/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
  } # end k
  
  
  #### k even
  if(k/2 ==floor(k/2)) { 
    lim<-dim(dat)[1]
    
    ## grouped
    for (i in 1:lim) {if (dat$resp0[i]<=(k)/2)dat$respgr[i]<-1
    if (dat$resp0[i]>=(k)/2+1) dat$respgr[i]<-2
    dat$obs[i]=i}
    
    Xgr<-as.matrix(dat[,pred])  ### !!!!!!
    respmgr<- matrix(0,dim(Xgr)[1],2)
    for (i in 1:dim(Xgr)[1]) respmgr[i,dat$respgr[i]]<-1
    ngrgr <- rep(1,dim(Xgr)[1])
    
    
    
    ##subsamples
    
    datdis<-dat[dat$resp0<=(k)/2,]  ## disagree
    #datdis$respresc<-(k+1)/2-datdis$resp  ## rescaled response inverted!
    datdis$respresc<-datdis$resp0 ###original
    datdis$indagr<-0   ### indicator agreement categories
    
    
    datagr<-dat[dat$resp0>=(k)/2+1,]  ## agree
    datagr$respresc<-  datagr$resp0-(k)/2  ## 
    datagr$indagr<-1
    
    datagdis<- rbind(datdis,datagr)  ### all observation without neutral
    indagr<-c(datdis$indagr,datagr$indagr)
    
    respdis<-datdis$respresc
    respagr<-datagr$respresc
    resp<-as.matrix(c(datdis$respresc,datagr$respresc))
    obsind<-c(datdis$obs,datagr$obs)
    
    ## design
    
    predwithout<-setdiff(pred, preddisp)
    predawithout<-setdiff(preda, preddisp)
    preddiswithout<-setdiff(preddis, preddisp)
    
    Xdis<-cbind(datdis[,preddisp],datdis[,predawithout],0*(datdis[,preddiswithout]))
    Xa<-cbind(datagr[,preddisp],0*datagr[,predawithout],(datagr[,preddiswithout]))
    
    Xdis<-as.matrix(Xdis)
    Xa<-as.matrix(Xa)
    X <-rbind(Xdis,Xa)
    #X<- as.matrix(X)  
    
    ### response matrix
    knew<-k/2
    k1<- knew-1
    n <-  dim(X)[1]
    p <-  dim(X)[2]
    
    respm<- matrix(0,n,knew)
    for (i in 1:n) respm[i,resp[i]]<-1
    ngr <- rep(1,n)
    
    ### only for output needed
    Xsepneut<-as.matrix(dat[,pred])  ### !!!!!!
    respmsepneut<- matrix(0,dim(Xsepneut)[1],2)
    #for (i in 1:dim(Xsepneut)[1]) respmsepneut[i,dat$respneut[i]]<-1
    ngsepneut <- rep(1,dim(Xsepneut)[1])
    
  } # end k
  
  
  newList <- list("respmgr"=respmgr, "ngrgr"= ngrgr,"Xgr"= Xgr,"respm"=respm, "ngr"= ngr,"X"= X,"indagr"=indagr, "datmod"=dat, "datdis"=datdis,"datagr"=datagr,
                  "obsind"=obsind, "respmsepneut"=respmsepneut,"ngsepneut"= ngsepneut,"Xsepneut"=Xsepneut)
  
  
  return(newList)}
####################




########################################
