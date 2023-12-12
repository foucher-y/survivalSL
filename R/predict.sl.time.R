
predict.sl.time <-function(object, ..., newdata=NULL, newtimes=NULL){
  
  if(is.null(newtimes)==FALSE){
    newtimes<-sort(unique(newtimes))
  }
  M <- length(object$method)
  
  FitVA <- vector("list", M+1)
  
  # names(FitVA) <- c(object$method, "sl")
  names(FitVA)<-c(names(object$model) ,"sl") # CS changement quand keep.pred ==FALSE

  for (i in 1:M) { FitVA[[i]] <- predict(object$models[[i]], newdata=newdata, newtimes=newtimes)$predictions }
  
  FitVA[[M+1]] <- matrix(0, nrow=dim(FitVA[[1]])[1], ncol=dim(FitVA[[1]])[2])
  

  if(object$metric$metric!="loglik.td"){
    w.sl=object$weights$values
    for (i in 1:M) { FitVA[[M+1]]  <- FitVA[[M+1]] + w.sl[i]*FitVA[[i]] }
    if(is.null(newtimes)) {time.pred <- object$times} else {time.pred <- newtimes}
  }
  else{
    if(is.null(newtimes)) {
      time.pred <- object$times
    } 
    else {
      time.pred <- newtimes
    }
    time.pred<-c(0,time.pred)
    
    N=dim(FitVA[[1]])[1]
    #w.sl <- object$weights
    w.sl <- object$weights$coefficients # Attention : A VERIFIER PAR CAMILLE
    
    Mm1=M-1
    
    par.alpha=w.sl[seq(from=1,to=2*Mm1-1, by=2)]
    par.beta=w.sl[seq(from=2,to=2*Mm1, by=2)]
    
    vect.w.times<-matrix(nrow = Mm1, ncol=length(time.pred))
    for( i in 1:Mm1){
      vect.w.times[i,]<-exp(par.alpha[i]+par.beta[i]*time.pred)
    }
    if(max(vect.w.times)==Inf){ #pour eviter des poids infinis
      co<-c()
      for (i in 1:(M-1)){
        co<-c(co,which(vect.w.times[i,]==Inf) )
      }
      # vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,vect.w.times[,co])
      vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,0)
    }
    
    vect.w.times2<-matrix(nrow = M, ncol=length(time.pred))
    for( i in 1:Mm1){
      vect.w.times2[i,]<-vect.w.times[i,]/(1+apply(vect.w.times,2,sum))
    }
    if(M>2){
      vect.w.times2[M,]<-1-apply(vect.w.times2[1:Mm1,],2,sum)
    }
    else{
      vect.w.times2[M,]<-1-vect.w.times2[1,]
    }
    
    for(j in 1:M){
      for (i in 1:dim(FitVA[[1]])[1]){
        if(min(FitVA[[j]][i,])==0){
          FitVA[[j]][i,FitVA[[j]][i,]==0]<-min(FitVA[[j]][i,which(FitVA[[j]][i,]!=0)])
        }
      }
      
    }

    
    
    .haz<-vector("list",M)
    for (i in 1:M){
      .Hazcum <- -log(FitVA[[i]]) #il faudra verifier que FitCV >0
      .Hazcum<-cbind(rep(0,N),.Hazcum)
      # .time <- time.pred[time.pred %in% unique(data.times) ]
      fmono <- apply(.Hazcum, FUN=function(H){mosaic::spliner(H ~ time.pred, monotonic = TRUE)}, MARGIN = 1)
      Dfmono <- lapply(fmono, FUN=function(fonction){mosaicCalc::D(fonction(t) ~ t)})
      .haz[[i]] <-t(sapply(1:N,FUN=function(ni){Dfmono[[ni]](time.pred)} ))
    }
    
    for (i in 1:M){
      if(min(.haz[[i]])<0){
        .haz[[i]]<-apply(.haz[[i]],2,FUN=function(l){ifelse(l<0,0,l)})
      }
      else{
        .haz[[i]]<-.haz[[i]]
      }
    }
    haz.sl<-array(dim = c(dim(.haz[[1]]),M))
    for (i in 1:M){
      haz.sl[,,i]<-.haz[[i]]*vect.w.times2[i,]
    }
    haz.sl<-rowSums(haz.sl, dims=2)
    
    
    # Haz.sl=t(apply(haz.sl, MARGIN=1,integration, x=time.pred))
    
    haz.slmono <- apply(haz.sl, FUN=function(H){mosaic::spliner(H ~ time.pred, monotonic = TRUE)}, MARGIN = 1)
    
    anti.haz.sl.fct=function(lehaz,vec_l,vec_u){
      anti.haz.sl.f <- function(l, u) {cubature::cubintegrate(lehaz, lower = l, upper = u, method = "pcubature")$integral}
      .i<-mapply(anti.haz.sl.f, vec_l, vec_u)
      return(.i)
    }
    
    .l <- time.pred[-length(time.pred)]
    .u <- time.pred[-1]
    
    # system.time(
    aera<-lapply(haz.slmono, anti.haz.sl.fct, .l,.u)
    # )
    
    Haz.sl=apply(cbind(rep(0,N),do.call(rbind,aera)),cumsum,MARGIN = 1)
    
    .pred=exp(-Haz.sl)
    
    
    # survs<-t(.pred)
    # # CS : MAJ pour eviter la boucle
    # timeVector=survfit(Surv(newdata[,times],newdata[,failures])~ 1 )$time
    # idx=findInterval(timeVector, time.event)
    # survs.long=survs[,pmax(1,idx)]
    # 
    
    FitVA$sl=t(.pred)[,-1] #pour enlever le decalage cree
    
    times=time.pred
  }
  
  
  
  
  return(list(predictions=FitVA,
              methods=c(object$methods,"sl"),
              times=time.pred))
}
