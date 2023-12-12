

predict.flexsurv <- function(object, ..., newdata=NULL, newtimes=NULL){

  .flex=object$model

  if(is.null(newdata)){
    if(is.null(newtimes)){
      .predlist<-summary(.flex, type = "survival", newdata=object$data, ci = F, se=F )
      .time.temp=.predlist[[1]]$time
    }
    else{
      .predlist<-summary(.flex, type = "survival", newdata=object$data, ci = F, se=F, t=newtimes)
      .time.temp=newtimes
    }
  }
  else{
    if(is.null(newtimes)){
      .predlist<-summary(.flex, type = "survival", newdata=newdata, ci = F, se=F)
      .time.temp=.predlist[[1]]$time
    }
    else{
      .predlist<-summary(.flex, type = "survival", newdata=newdata, ci = F, se=F, t=newtimes)
      .time.temp=newtimes
    }
  }

  .pred=matrix(nrow=length(.predlist), ncol=length(.predlist[[1]]$time))

  for (i in 1:length(.predlist)){
    .pred[i,]=.predlist[[i]]$est
  }


  return(list(times=.time.temp, predictions=.pred))
}
