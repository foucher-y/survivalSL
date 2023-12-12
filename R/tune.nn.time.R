

# CAMILLE : JE NE COMPRENDS CETTE FONCTION, IL Y A DES ARGUMENTS "METRIC" et "CL" NON-UTILISES -> J'AI SIMPLIFIE  COMME LES AUTRES FONCTIONS "TUNE" 
# CS : metric ne servait a rier 
# CS : option CL non utilisee actuellement car fonction pas parallelisable actuellement 
tune.nn.time <- function(times, failures, group=NULL, 
                           cov.quanti=NULL, cov.quali=NULL, 
                           data, cv=10, 
                           n.nodes,
                           decay,
                           batch.size,
                           epochs){ # CAMILLE : J'AI CHANGE LE NOM DE LA FONCTION POUR HOMOGENEITE

  
  
  
  # CS ici pas besoin de faire un check si group ou cov.quanti ou cov.quali == NULL
  data.nn <- data[,c(times, failures, group, cov.quanti, cov.quali)]
  

# CS : a faire a chaque fois car fonction "manuelle" pour faire la CV 
  sample_id <- sample(nrow(data.nn))
  folds <- cut(seq(1,nrow(data.nn)), breaks=cv, labels=FALSE)
  folds_id <- folds[sample_id]
  data.nn$folds <- folds_id

  
  .f  <- as.formula(paste("Surv(", times, ",", failures, ")", "~."))
  .time <- sort(unique(data.nn[,times]))
  .grid <-  expand.grid(n.nodes=n.nodes, decay=decay, batch.size=batch.size,
                        epochs=epochs)
  # .ibs <- rep(-99, dim(.grid)[1]) #a garder ?



  .CVtune<-vector("list",cv*dim(.grid)[1]) 
  # CAMILLE : Il FAUT DONC ICI VERIFIER ICI QUE DATA$FOLD A LE MEME NB DE PARTITION de 1 A CV
  #  CS on supprime l'option foldid et on creee nous meme a chaque fois donc pas de check a faire
  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      .CVtune[[l]]<-list(train=data.nn[data.nn$folds!=k, ], valid=data.nn[data.nn$folds==k, ], grid=.grid[j,])
      l=l+1
    }
  }

  # CS ici il faut bien garder l'option newtimes
  nn.time.par<-function(xx, times, failures, group, cov.quanti, cov.quali,newtimes){
    
    n.nodes=xx$grid$n.nodes
    decay=xx$grid$decay
    batch_size=xx$grid$batch.size
    epochs=xx$grid$epochs
    
    data=xx$train
    newdata=xx$valid
    
    
    if(!(is.null(group))){
      .data <- data[,c(times, failures, group, cov.quanti, cov.quali)]}   else{
        .data <- data[,c(times, failures, cov.quanti, cov.quali)] }
    
    .f  <- as.formula(paste("Surv(", times, ",", failures, ")", "~."))
    
  
    .deepsurv <- deepsurv(.f, data = .data,  verbose = FALSE, num_nodes=n.nodes,
                          weight_decay=decay, num_workers = 0L,batch_size=as.integer(batch_size),
                          epochs=as.integer(epochs))
    
    # .time <- .time.deepsurv <- sort(unique(.data[,times]))
    .time<-sort(unique(.data[,times]))
    
    
    .newdata <- data.frame(newdata[,c(group, cov.quanti, cov.quali)])
    # .pred <- predict(.deepsurv, newdata=.newdata)
    .pred <- predict(.deepsurv, newdata=newdata)
    .time.deepsurv<-as.numeric(dimnames(.pred)[[2]]) #MAJ CS 20210818

    if(!is.null(newtimes)) { #maj CS 20210818
      .pred.deepsurv <- cbind(rep(1, dim(.pred)[1]), .pred)
      # .pred.deepsurv <-  .pred[,1:(length(.times)+1)] # CS : correction bug
      .time.deepsurv <- c(0, .time.deepsurv)
      
      # CS : MAJ pour eviter la boucle
      idx=findInterval(newtimes, .time.deepsurv)
      .pred=.pred.deepsurv[,pmax(1,idx)]
      
      # .pred <- matrix(-99, nrow = dim(.pred)[1], ncol = length(newtimes))
      # .pred[,1] <- matrix(.pred.deepsurv[ ,.time.deepsurv <= newtimes[1]], ncol= sum(.time.deepsurv<=newtimes[1]) )[,sum(.time.deepsurv<=newtimes[1])]
      # for (i in 1:length(newtimes)) {
      #   .pred[,i] <- .pred.deepsurv[,.time.deepsurv<=newtimes[i]][,sum(.time.deepsurv<=newtimes[i])]
      # }
      
      # .pred<-matrix(nrow=dim(.pred.deepsurv)[1], ncol=length(newtimes))
      # colnames(.pred)<-newtimes
      # .pred[,paste0(.time.deepsurv[which(.time.deepsurv%in%newtimes)])]<-.pred.deepsurv[,.time.deepsurv%in%newtimes]
      # .pred[,dim(.pred)[2]]<-.pred.deepsurv[,.time.deepsurv<=newtimes[length(newtimes)]][,sum(.time.deepsurv<=newtimes[length(newtimes)])]
      # if(is.na(.pred[1,1])==T){
      #   .pred[,1]<-.pred.deepsurv[,1]
      # }
      # .pred<-t(na.approx(t(.pred), method="constant"))
      
      
      .time <- newtimes
    }
    
    
    # .obj <- list(times=.time, model=.deepsurv, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, predictions=.pred)
    #
    # class(.obj) <- "nn.time"
    
    return(.pred)
  }
  
  # if(parallel==T){
  #    .preFIT<-parLapply(cl,.CVtune, nn.time.par,times=times, failures=failures, group=group,
  #                   cov.quanti=cov.quanti, cov.quali=cov.quali,newtimes=.time, seed=seed,chunk.size =1)
  #
  # }
  # else{
  .preFIT<-list()
  for (i in 1:length(.CVtune)){
    .preFIT[[i]]<-nn.time.par(.CVtune[[i]], times=times, failures=failures, group=group,
                                cov.quanti=cov.quanti, cov.quali=cov.quali,newtimes=.time)
  }
    .preFIT<-lapply(.CVtune, nn.time.par, times=times, failures=failures, group=group,
                 cov.quanti=cov.quanti, cov.quali=cov.quali,newtimes=.time)

  # }

  # .FitCV<-vector("list",dim(.grid)[1])
  .FitCV <- replicate(dim(.grid)[1], matrix(NA, nrow = length(data[,times]),
                                            ncol = length(.time)), simplify=F)
  l<-1
  for (k in 1:cv){
    for (j in 1:dim(.grid)[1]){
      # print(j)
      .FitCV[[j]][data.nn$folds==k,] <- .preFIT[[l]]
      l<-l+1
    }
  }

  nn.best.measure<-  function(prediction.matrix, times, failures, data){ 
  # CS : suppression des premieres lignes inutiles
    # rajout d'une boucle qui controle que prediction de la survie !=0 pour eviter les pb avec le log par la suite
    for (i in 1:dim(prediction.matrix)[1]){
      if(min(prediction.matrix[i,])==0){
        prediction.matrix[i,prediction.matrix[i,]==0]<-min(prediction.matrix[i,which(prediction.matrix[i,]!=0)])
      }
    }
    
    time.pred <- sort(unique(data[,times]))
    data.times <- data[,times]
    data.failures <- data[,failures]
    
    .surv <- prediction.matrix[, time.pred %in% unique(data.times[data.failures==1]) ]
    .time <- time.pred[time.pred %in% unique(data.times[data.failures==1]) ]
    .haz <- t(apply(.surv, FUN = function(s) { differentiation(x=.time, fx=-1*log(s)) }, MARGIN=1))
    .indic <- t(sapply(data.times, FUN = function(x) {1*(.time==x)} ))
    .haz.indic <- .haz * .indic
    
    .surv <- prediction.matrix
    .indic <- t(sapply(data.times, FUN = function(x) {1*(time.pred==x)} ))
    .surv.indic <- .surv * .indic
    
    .haz.indiv <- apply(.haz.indic, FUN = "sum", MARGIN=1)[data.failures==1]
    .surv.indiv <- apply(.surv.indic, FUN = "sum", MARGIN=1)
    
    RET=(-1*(sum(log( pmax(.haz.indiv, min(.haz.indiv[.haz.indiv!=0])) )) 
             + sum(log(.surv.indiv))))
    
    return(RET)
    
  }
  
  
  
 # if(parallel==T){
 #   .measure<-parSapply(cl,.FitCV,nn.best.measure,times=times, failures=failures,data=data,
  #                  chunk.size = 1)

 # }
 # else{
    .measure<-sapply(.FitCV, nn.best.measure, times=times, failures=failures, data=data.nn)

 # }

    .res <- data.frame(n.nodes = .grid[,1], decay = .grid[,2], batch.size=.grid[,3],
                       epochs=.grid[,4] , measure = .measure)
    

    .mini<-.res[which(.res$measure==min(.res$measure, na.rm=TRUE) & is.na(.res$measure)==FALSE),]
    .mini<-.mini[1,]

    
    return( list(optimal=list(n.nodes=.mini$n.nodes,
                              decay=.mini$decay,
                              batch.size=.mini$batch.size,
                              epochs=.mini$epochs),
                 results=.res ))
}
