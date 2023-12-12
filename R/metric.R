

metric <- function(times, failures, data, prediction.matrix, prediction.times, metric, pro.time=NULL, ROC.precision=seq(.01, .99, by=.01))
{
  data.times <- data[,times]
  data.failures <- data[,failures]
  timeVector <- survfit(Surv(data[,times],data[,failures])~ 1 )$time
  obj_surv <- Surv(data.times, data.failures)

  time <- obj_surv[, 1]
  ot <- order(time)
  cens <- obj_surv[ot, 2]
  time <- time[ot]

  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE)
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf
  # csurv_btime <- predict(hatcdist, times = timeVector, type = "surv")
  csurv_btime <- predict(hatcdist, times = sort(prediction.times), type = "surv")
  csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
  csurv_btime[csurv_btime == 0] <- Inf

  survs <- t(prediction.matrix)[,ot]

  time.pred <- sort(unique(data[,times]))
  
  # ici ajour:
  timeVector<-sort(unique(prediction.times))
  
  if(is.null(pro.time)) {pro.time <- median(data.times)}
  
  switch(metric,
         ibs={
           bsc <- sapply(1:length(timeVector), FUN = function(j)
           {
             help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
             help2 <- as.integer(time > timeVector[j])
             return(mean((0 - survs[j, ])^2 * help1 * (1/csurv) + (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j])))
           })

           idx <- 2:length(timeVector)
           RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx])/2)
           RET <- RET/diff(range(timeVector))
           RET <- as.matrix(RET)
         },
         bs={
           j <- length(timeVector[which(timeVector<=pro.time)])
           help1 <- as.integer(time <= timeVector[j] &  obj_surv[ot,2] == 1)
           help2 <- as.integer(time > timeVector[j])
           bs <- mean((0 - survs[j, ])^2 * help1 * (1/csurv) + (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j]))
           bs <- as.numeric(bs)
           RET <- bs
         },
         ibll={
           survs[which(survs==0)]<-10**-7
           survs[which(survs==1)]<-1-10**-7
           bll <- sapply(1:length(timeVector), FUN = function(j)
           {
             help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
             help2 <- as.integer(time > timeVector[j])
             RET=-mean(log(1 - survs[j, ]) * help1 * (1/csurv) + log( survs[j, ]) * help2 * (1/csurv_btime[j]))
             RET=ifelse(is.nan(RET)==T,0,RET)
             return(RET)
           })

           idx <- 2:length(timeVector)
           RET <- diff(timeVector) %*% ((bll[idx - 1] + bll[idx])/2)
           RET <- RET/diff(range(timeVector))
           RET <- as.matrix(RET)
         },
         bll={
           survs[which(survs==0)]<-10**-7
           survs[which(survs==1)]<-1-10**-7
           j <- length(timeVector[which(timeVector<=pro.time)])
           help1 <- as.integer(time <= timeVector[j] &  obj_surv[ot,2] == 1)
           help2 <- as.integer(time > timeVector[j])
           bll <- -mean(log(1- survs[j, ]) * help1 * (1/csurv) + log(survs[j, ]) * help2 * (1/csurv_btime[j]))
           bll <- as.numeric(bll)
           RET <- bll
         },
         loglik={  # CAMILLE : Ne fonctionne pas
           if(min(prediction.matrix)==0){
             for (i in 1:dim(prediction.matrix)[1]){
               if(min(prediction.matrix[i,])==0){
                 # print(i)
                 prediction.matrix[i,prediction.matrix[i,]==0]<-min(prediction.matrix[i,which(prediction.matrix[i,]!=0)])
               }
             }
           }

           
           time.pred <- sort(unique(data[,times]))
           # time.pred <- sort(unique(prediction.times))
           data.times <- data[,times]
           data.failures <- data[,failures]
           .surv <- prediction.matrix[, time.pred %in% unique(data.times[data.failures==1]) ]
           .time <- time.pred[time.pred %in% unique(data.times[data.failures==1]) ]
           .haz <- t(apply(.surv, FUN = function(s) { differentiation(x=.time, fx=-1*log(s)) }, MARGIN=1))
           .indic <- t(sapply(data.times, FUN = function(x) {1*(.time==x)} ))
           .haz.indic <- .haz * .indic
           
           .surv <- prediction.matrix
           .indic <- t(sapply(data.times, FUN = function(x) {1*(time.pred==x)} ))
           # .indic <- t(sapply(data.times, FUN = function(x) {1*(sort(unique(prediction.times))==x)} ))
           .surv.indic <- .surv * .indic
           
           .haz.indiv <- apply(.haz.indic, FUN = "sum", MARGIN=1)[data.failures==1]
           .surv.indiv <- apply(.surv.indic, FUN = "sum", MARGIN=1)
           # if(min(.surv.indiv==0)){
           #   for (i in 1:dim(.surv.indiv)[1]){
           #     if(min(.surv.indiv[i,])==0){
           #       # print(i)
           #       .surv.indiv[i,.surv.indiv[i,]==0]<-min(.surv.indiv[i,which(prediction.matrix[i,]!=0)])
           #     }
           #   }
           # }
           RET=((sum(log( pmax(.haz.indiv, min(.haz.indiv[.haz.indiv!=0])) )) + sum(log(.surv.indiv))))
           
         },
         auc={
           
           .data <- data.frame(times=data[,times], failures=data[,failures], variable=1-prediction.matrix[,prediction.times>=pro.time][,1])
           RET <- roc.time(times="times", failures="failures", variable="variable", confounders=~1, data=.data,
                           pro.time=pro.time, precision=ROC.precision)$auc
         },
        ribs={
          timeVector <- timeVector[timeVector<=pro.time]
          bsc <- sapply(1:length(timeVector), FUN = function(j)
          {
            help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
            help2 <- as.integer(time > timeVector[j])
            return(mean((0 - survs[j, ])^2 * help1 * (1/csurv) + (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j])))
          })
          idx <- 2:length(timeVector)
          RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx])/2)
          RET <- RET/diff(range(timeVector))
          RET <- as.matrix(RET)
        },
        ribll={
          timeVector <- timeVector[timeVector<=pro.time]
          bll <- sapply(1:length(timeVector), FUN = function(j)
          {
            help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
            help2 <- as.integer(time > timeVector[j])
            RET=-mean(log(1 - survs[j, ]) * help1 * (1/csurv) + log( survs[j, ]) * help2 * (1/csurv_btime[j]))
            RET=ifelse(is.nan(RET)==T,0,RET)
            return(RET)
          })

          idx <- 2:length(timeVector)
          RET <- diff(timeVector) %*% ((bll[idx - 1] + bll[idx])/2)
          RET <- RET/diff(range(timeVector))
          RET <- as.matrix(RET)
        },
  )
  return(as.numeric(RET))
}

