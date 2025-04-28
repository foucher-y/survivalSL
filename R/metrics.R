metrics <- function(metric, times, failures, data, survivals.matrix, hazards.matrix=NULL,prediction.times, pro.time=NULL, ROC.precision=seq(.01, .99, by=.01))
{


  if (missing(times)) stop("The 'times' argument is required.")
  if (missing(failures)) stop("The 'failures' argument is required.")
  if (missing(survivals.matrix)) stop("The 'prediction.matrix' argument is required.")
  if (missing(prediction.times)) stop("The 'prediction.times' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(metric)) stop("The 'metric' argument is required.")


  if((metric=="ll") & is.null(hazards.matrix))stop("The 'hazards.matrix' argument is required.")


  if (!(is.character(times))) stop("The 'times' argument must be character.")
  if (! (is.character(failures))) stop("The 'failures' argument must be character.")
  if (! (is.character(metric))) stop("The 'metric' argument must be character.")

  if(!(times %in% colnames(data)))stop(paste("The variable",times, "does not exist."))

  if(!(failures %in% colnames(data)))stop(paste("The variable",failures, "does not exist."))

  if(min(metric %in% c("uno_ci","p_ci","bs","ibs","ibll","bll", "ribs","ribll","auc","ll"))==0){
    stop("The argument \"metric\" must be Brier score (bs),
         Pencina concordance index (p_ci),
         Uno concordance index (uno_ci),
         integrated Brier score (ibs), the binomial log-likelihood (bll),
         the integrated binomial log-likelihood (ibll), the restricted ibs (ribs),
         the restricted ibll (ribll), the log-likelihood (loglik),
         the area under the ROC curve (auc), or the log-likelihood (ll)")
  }


  if(is.null(pro.time) & !(metric %in% c("ll","ibll","ibs"))) {pro.time <- median(data[[times]])}

  data.times <- data[[times]]  # our times
  data.failures <- data[[failures]] # our events
  obj_surv <- Surv(data.times, data.failures)

  time <- obj_surv[, 1]
  ot <- order(time)
  cens <- obj_surv[ot, 2]
  time <- time[ot]

  .temp <- survfit(Surv(time, cens==0) ~ 1)
  .csurv <- summary(.temp, times=time, extend=TRUE)$surv  # proba of censorship for all data times

  .csurv[.csurv == 0] <- Inf

  if(!(is.null(pro.time))){
    .csurv.pro.time <- summary(.temp, times=pro.time, extend=TRUE)$surv
    .csurv.pro.time[.csurv.pro.time == 0] <- Inf
  }

  survivals.matrix <- survivals.matrix[ot,] # reorder our survivals matrix

  hazards.matrix <- hazards.matrix[ot,] # same for our hazards matrix

  .data_bis<-data
  data<-data[ot,] # same with our data





  switch(metric,
         p_ci={
           status <- ifelse((data[[failures]]==1) & (data[[times]]>pro.time),0,data[[failures]])
           time <- ifelse(data[[times]]>pro.time,pro.time,data[[times]])
           length<-dim(survivals.matrix[,prediction.times>=pro.time])
           if(is.null(length)){
             predicted<-survivals.matrix[,prediction.times>=pro.time]
           }else{
             predicted<-survivals.matrix[,prediction.times>=pro.time][,1]
           }
           permissible <- 0 # comparable pairs
           concord <- 0
           n <- length(time)
           for (i in 1:n) {
             for (j in (i + 1):n) {
               if(time[i] < time[j] & status[i] == 1 ){
                 permissible<-permissible + 1
                 if(predicted[i] < predicted[j]){
                   concord <- concord + 1
                 }
                 if(predicted[i]==predicted[j]){
                   concord <- concord + 0.5
                 }
               }

               if(status[i] == 1 & status[j] == 0 & time[i]==time[j] ){
                 permissible<-permissible + 1
                 if(predicted[i] < predicted[j]){
                   concord <- concord + 1
                 }
               }



             }

           }
           RET<-ifelse(permissible==0,0,concord/permissible)

         },


         ibs={
           bsc <- sapply(1:length(prediction.times), FUN = function(j)
           {pro.time<-prediction.times[j]
           .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
           .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
           score_risque<-survivals.matrix[,j]
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time ,1,0)
           bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
           bs <- as.numeric(bs)
           return(bs)

           })

           idx <- 2:length(prediction.times)
           RET <- diff(prediction.times) %*% ((bsc[idx - 1] + bsc[idx])/2)
           RET <- RET/diff(range(prediction.times))
           RET <- as.matrix(RET)
         },
         bs={
           length<-dim(survivals.matrix[,prediction.times>=pro.time])
           if(is.null(length)){
             score_risque<-survivals.matrix[,prediction.times>=pro.time]
           }else{
             score_risque<-survivals.matrix[,prediction.times>=pro.time][,1]
           }
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time,1,0)
           bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
           RET <- as.numeric(bs)

         },
         bll={
          length<-dim(survivals.matrix[,prediction.times>=pro.time])
           if(is.null(length)){
             pro.time.surv<-survivals.matrix[,prediction.times>=pro.time]
           }else{
             pro.time.surv<-survivals.matrix[,prediction.times>=pro.time][,1]
           }
           pro.time.surv[which(pro.time.surv==0)]<-10**-7
           pro.time.surv[which(pro.time.surv==1)]<-1-10**-7
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time,1,0)
           bll <- -mean(log(1-pro.time.surv) * help1 * (1/.csurv) + log(pro.time.surv) * help2 * (1/.csurv.pro.time))
           RET <- as.numeric(bll)
         },
         ibll={
           bllc <- sapply(1:length(prediction.times), FUN = function(j)
           {pro.time<-prediction.times[j]
           .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
           .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
           score_risque<-survivals.matrix[,j]
           score_risque[which(score_risque==0)]<-10**-7
           score_risque[which(score_risque==1)]<-1-10**-7
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time ,1,0)
           bll <- -mean(log(1-score_risque) * help1 * (1/.csurv) + log(score_risque) * help2 * (1/.csurv.pro.time))
           return(as.numeric(bll)) })
           idx <- 2:length(prediction.times)
           RET <- diff(prediction.times) %*% ((bllc[idx - 1] + bllc[idx])/2)
           RET <- RET/diff(range(prediction.times))
           RET <- as.matrix(RET)
         },

         ribll={
           prediction.times <- prediction.times[prediction.times <= pro.time]
           survivals.matrix<-survivals.matrix[,prediction.times<=pro.time]
           bllc <- sapply(1:length(prediction.times), FUN = function(j)
           {pro.time<-prediction.times[j]
           .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
           .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
           score_risque<-survivals.matrix[,j]
           score_risque[which(score_risque==0)]<-10**-7
           score_risque[which(score_risque==1)]<-1-10**-7
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time ,1,0)
           bll <- -mean(log(1-score_risque) * help1 * (1/.csurv) + log(score_risque) * help2 * (1/.csurv.pro.time))
           return(as.numeric(bll)) })
           idx <- 2:length(prediction.times)
           RET <- diff(prediction.times) %*% ((bllc[idx - 1] + bllc[idx])/2)
           RET <- RET/diff(range(prediction.times))
           RET <- as.matrix(RET)

         },

         ribs={
           prediction.times <- prediction.times[prediction.times <= pro.time]
           survivals.matrix<-survivals.matrix[,prediction.times<=pro.time]
           bsc <- sapply(1:length(prediction.times), FUN = function(j)
           {pro.time<-prediction.times[j]
           .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
           .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
           score_risque<-survivals.matrix[,j]
           help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
           help2 <- ifelse(data[[times]] > pro.time ,1,0)
           bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
           bs <- as.numeric(bs)
           return(bs)

           })

           idx <- 2:length(prediction.times)
           RET <- diff(prediction.times) %*% ((bsc[idx - 1] + bsc[idx])/2)
           RET <- RET/diff(range(prediction.times))
           RET <- as.matrix(RET)
         },
         auc={
            length<-dim(survivals.matrix[,prediction.times>=pro.time])
           if(is.null(length)){
             pro.time.surv<-survivals.matrix[,prediction.times>=pro.time]
           }else{
             pro.time.surv<-survivals.matrix[,prediction.times>=pro.time][,1]
           }
           .data <- data.frame(times=data[[times]], failures=data[[failures]],
                               variable=(1-pro.time.surv))
           RET <- roc(times="times", failures="failures", variable="variable",
                      confounders=~1, data=.data,
                      pro.time=pro.time, precision=ROC.precision)$auc
         },

         uno_ci={
           status <- ifelse((data[[failures]]==1) & (data[[times]]>pro.time),0,data[[failures]])
           time <- ifelse(data[[times]]>pro.time,pro.time,data[[times]])
           length<-dim(survivals.matrix[,prediction.times>=pro.time])
           if(is.null(length)){
             predicted<-survivals.matrix[,prediction.times>=pro.time]
           }else{
             predicted<-survivals.matrix[,prediction.times>=pro.time][,1]
           }
           permissible <- 0 # comparable pairs
           concord <- 0
           n <- length(time)
           for (i in 1:n) {
             for (j in (i + 1):n) {
               if(time[i] < time[j] & status[i] == 1 ){
                 permissible<-permissible + 1/(.csurv[i])**2
                 if(predicted[i] < predicted[j]){
                   concord <- concord + 1/(.csurv[i])**2
                 }
                 if(predicted[i]==predicted[j]){
                   concord <- concord + 0.5*(1/(.csurv[i])**2)
                 }
               }

               if(status[i] == 1 & status[j] == 0 & time[i]==time[j] ){
                 permissible<-permissible + 1/(.csurv[i])**2
                 if(predicted[i] < predicted[j]){
                   concord <- concord + 1/(.csurv[i])**2
                 }
               }



             }

           }
           RET<-ifelse(permissible==0,0,concord/permissible)

         },

         ll={
           .Data<-data
           .Data$id<-1:nrow(.Data)
           .Data <- .Data[order(.Data[[times]]), ]
           .idx<-findInterval(unique(.Data[[times]]),prediction.times)
           .b<-1
           .result<-sapply(2:nrow(.Data),function(i){
             value<-.Data[[times]][i]
             if(value!=.Data[[times]][i-1]){
               .b<<-.b+1
             }
             return(.b)
           })
           .result<-c(1,.result)
           .Data$rank_time<-.result
           .Data <- .Data[order(.Data$id), ]
           .survivals.matrix_bis<-survivals.matrix[,.idx]
           .hazards.matrix_bis<-hazards.matrix[,.idx]

           .Data$surv<-sapply(1:nrow(.Data), function(i){
             .b<-.survivals.matrix_bis[i,(.Data$rank_time[i])]
             return(.b)
           })

           .Data$haz<-sapply(1:nrow(.Data), function(i){
             .b<-.hazards.matrix_bis[i,(.Data$rank_time[i])]
             return(.b)
           })

           .Data[[failures]]<-ifelse(is.na(.Data$haz),0,.Data[[failures]])
           .Data$haz<-ifelse(is.na(.Data$haz),10**-7,.Data$haz)
           .Data$haz<-ifelse(is.infinite(.Data$haz),10**-7,.Data$haz)




           RET<-sum(.Data[[failures]]*log(.Data$haz)+log(.Data$surv))


         })
  data<-.data_bis
  return(as.numeric(RET))
}




