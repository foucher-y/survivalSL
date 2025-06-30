predict.sltime <-function(object, newdata=NULL, newtimes=NULL, ...){

  M<-length(object$methods)

  Fit<- vector("list", M+1)

  if(is.null(newdata) & is.null(newtimes)){
    for (i in 1:M) { Fit[[i]] <- (predict(object$models[[i]]))$predictions }
    Fit[[M+1]]<-object$predictions[,-1]
    time.pred<-object$times[-1]
    } else {
    time.pred<-newtimes

    if(is.null(newtimes)) {time.pred <- object$times[-1]}

    if(is.null(newdata)){newdata <- object$data}

    for (i in 1:M) { Fit[[i]] <- (predict(object$models[[i]], newdata=newdata, newtimes=time.pred))$predictions}

    w.sl <- object$weights$values

    weighted_matrices <- mapply(function(mat, weight) mat * weight, Fit[-(M+1)], w.sl, SIMPLIFY = FALSE)
    Fit[[M+1]] <- Reduce("+", weighted_matrices)

  }



  names(Fit)<-c(names(object$model), "sl")


  return(list(predictions=Fit,
              times=time.pred))
}



