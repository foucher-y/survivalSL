

predict.snn <- function(object, ..., newdata=NULL, newtimes=NULL){

  .times <- object$times

  if(is.null(newdata))  {
    .pred <- object$predictions
  }
  else {

    .var <- c(object$cov.quanti, object$cov.quali)

    .newdata<-newdata[,.var]

    .pred <- predict(object$model, newdata=.newdata)
    .time.deepsurv<-as.numeric(dimnames(.pred)[[2]])

    idx=findInterval(.times, .time.deepsurv)
    .pred=.pred[,pmax(1,idx)]
  }


  if(!is.null(newtimes)) {
    .pred.deepsurv <- cbind(rep(1, dim(.pred)[1]), .pred)

        .time.deepsurv <- c(0, .times)

    idx=findInterval(newtimes, .time.deepsurv)
    .pred=.pred.deepsurv[,pmax(1,idx)]

    .times <- newtimes
  }

  return(list(times=.times, predictions=.pred))
}
