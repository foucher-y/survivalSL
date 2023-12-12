
predict.rf.time <- function(object, ..., newdata=NULL, newtimes=NULL){  #CAMILLE: Il faudrait verifier pourquoi les courbes ne sont pas tout e fait en excalier, mais avec des contre-marches non-verticales.
  if(is.null(newdata))  {
    .survival <- object$predictions
    .time.interest <- object$times
  }
  else {
    .pred.rf <- predict(object$model, newdata = newdata)
    .survival <- cbind(rep(1, dim(.pred.rf$survival)[1]), .pred.rf$survival)
    .time.interest <- c(0, .pred.rf$time.interest)
    # .pred <- matrix(-99, nrow = dim(.survival)[1], ncol = length(object$times))
    # .pred[,1] <- matrix(.survival[ ,.time.interest <= object$times[1]], ncol= sum(.time.interest<=object$times[1]) )[,sum(.time.interest<=object$times[1])]
    # for (i in 1:length(object$times)) {
    #   .pred[,i] <- .survival[,.time.interest<=object$times[i]][,sum(.time.interest<=object$times[i])]
    # }

    # .pred<-matrix(nrow=dim(.survival)[1], ncol=length(object$times))
    # colnames(.pred)<-object$times
    # .pred[,paste0(.time.interest[which(.time.interest%in%object$times)])]<-.survival[,.time.interest%in%object$times]
    # .pred[,dim(.pred)[2]]<-.survival[,.time.interest<=object$times[length(object$times)]][,sum(.time.interest<=object$times[length(object$times)])]
    # if(is.na(.pred[1,1])==T){
    #   .pred[,1]<-.survival[,1]
    # }
    # .pred<-t(na.approx(t(.pred), method="constant"))
    
    # CS : MAJ pour eviter la boucle
    idx=findInterval(object$times,.time.interest)
    .pred=.survival[,pmax(1,idx)]

    .survival <- .pred
    .time.interest <- object$times
  }


  if(!is.null(newtimes)) {
    .survival <- cbind(rep(1, dim(.survival)[1]), .survival)
    .time.interest <- c(0, .time.interest)
    
    # CS : idem MAJ sans boucle
    idx=findInterval(newtimes,.time.interest)
    .pred=.survival[,pmax(1,idx)]
    
    # .pred <- matrix(-99, nrow = dim(.survival)[1], ncol = length(newtimes))
    # .pred[,1] <- matrix(.survival[ ,.time.interest <= newtimes[1]], ncol= sum(.time.interest<=newtimes[1]) )[,sum(.time.interest<=newtimes[1])]
    # for (i in 1:length(newtimes)) {
    #   .pred[,i] <- .survival[,.time.interest<=newtimes[i]][,sum(.time.interest<=newtimes[i])]
    # }

    # .pred<-matrix(nrow=dim(.survival)[1], ncol=length(newtimes))
    # colnames(.pred)<-newtimes
    # .pred[,paste0(.time.interest[which(.time.interest%in%newtimes)])]<-.survival[,.time.interest%in%newtimes]
    # .pred[,dim(.pred)[2]]<-.survival[,.time.interest<=newtimes[length(newtimes)]][,sum(.time.interest<=newtimes[length(newtimes)])]
    # if(is.na(.pred[1,1])==T){
    #   .pred[,1]<-.survival[,1]
    # }
    # .pred<-t(na.approx(t(.pred), method="constant"))

    .survival <- .pred
    .time.interest <- newtimes
  }

  return(list(times=.time.interest, predictions=.survival))
}

