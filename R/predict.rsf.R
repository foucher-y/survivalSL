
predict.rsf <- function(object, ..., newdata=NULL, newtimes=NULL){
  if(is.null(newdata))  {
    .survival <- object$predictions
    .time.interest <- object$times
  }
  else {
    .pred.rf <- predict(object$model, newdata = newdata)
    .survival <- cbind(rep(1, dim(.pred.rf$survival)[1]), .pred.rf$survival)
    .time.interest <- c(0, .pred.rf$time.interest)

    idx=findInterval(object$times,.time.interest)
    .pred=.survival[,pmax(1,idx)]

    .survival <- .pred
    .time.interest <- object$times
  }


  if(!is.null(newtimes)) {
    .survival <- cbind(rep(1, dim(.survival)[1]), .survival)
    .time.interest <- c(0, .time.interest)

    idx=findInterval(newtimes,.time.interest)
    .pred=.survival[,pmax(1,idx)]

    .survival <- .pred
    .time.interest <- newtimes
  }

  return(list(times=.time.interest, predictions=.survival))
}

