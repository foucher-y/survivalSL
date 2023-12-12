
#object=.deepsurv
#newdata=.data

predict.nn.time <- function(object, ..., newdata=NULL, newtimes=NULL){
  
  .times <- object$times

  if(is.null(newdata))  { 
    .pred <- object$predictions 
  }
  else {
    # if(!(is.null(object$group))){
    #   .newdata <- newdata[,c(object$group, object$cov.quanti, object$cov.quali)]  }
    # else{
    #   .newdata <- newdata[,c(object$cov.quanti, object$cov.quali)]
    # }
    
    .var <- c(object$cov.quanti, object$cov.quali)
    
    .newdata<-newdata[,.var]

    .pred <- predict(object$model, newdata=.newdata)
    .time.deepsurv<-as.numeric(dimnames(.pred)[[2]])
    
    idx=findInterval(.times, .time.deepsurv)
    .pred=.pred[,pmax(1,idx)]
    # .pred<-.pred[,-c(1,dim(.pred)[2])]

  }

  
  if(!is.null(newtimes)) {
    .pred.deepsurv <- cbind(rep(1, dim(.pred)[1]), .pred)
    # .pred.deepsurv <-  .pred[,1:(length(.times)+1)] # CS : correction bug
    .time.deepsurv <- c(0, .times)
    
    # CS : MAJ pour eviter la boucle
    idx=findInterval(newtimes, .time.deepsurv)
    .pred=.pred.deepsurv[,pmax(1,idx)]

    .times <- newtimes
  }

  return(list(times=.times, predictions=.pred))
}
