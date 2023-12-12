
predict.cox <- function(object, ..., newdata=NULL, newtimes=NULL){
  group=object$group
  cov.quanti=object$cov.quanti
  cov.quali=object$cov.quali
  if(is.null(newdata))  {
    .pred.temp <- object$predictions
    .time.temp <- object$times
  }
  else {
    if(is.null(cov.quanti)==T & is.null(cov.quali)==T){
      .coxphsurv<-survfit(object$model, newdata = newdata,se.fit = F)
      
      .sumcoxphsurv<-summary(.coxphsurv, times=sort(unique(object$times)))
      .pred.temp <- t(.sumcoxphsurv$surv)
      .time.temp <- .sumcoxphsurv$time
    }
    else{
      .bs=NULL
      .bin=NULL
      if(is.null(cov.quanti)==F){
        .bs <- eval(parse(text=paste("cbind(",
                                     paste("bs(newdata$",cov.quanti,",df=3)", collapse = ", ")
                                     ,")") ) )
      }
      if(is.null(cov.quali)==F){
        .bin <- eval(parse(text=paste("cbind(",  paste("newdata$", cov.quali, collapse = ", "), ")") ) )
      }
      # .bs <- eval(parse(text=paste("cbind(",paste("newdata$", object$cov.quanti, collapse = ", "),")") ) )
      # .bin <- eval(parse(text=paste("cbind(",  paste("newdata$", object$cov.quali, collapse = ", "), ")") ) )
      
      .cov <- cbind(.bs,.bin)
      if(!(is.null(object$group))){
        .x <- cbind(newdata[,object$group], .cov, .cov * newdata[,object$group])
      }
      else{
        .x <- .cov
      }
      if(class(object$model)[1]=="coxph"){
        .lp.coxph <- predict(object$model, newdata = newdata, type="lp")
        .pred.temp <- exp(matrix(exp(.lp.coxph)) %*% t(as.matrix(-1*object$hazard)))
        .time.temp <- object$times
      }
      else{
        .lp.lasso <- predict(object$model, newx = .x)
        .pred.temp <- exp(matrix(exp(.lp.lasso)) %*% t(as.matrix(-1*object$hazard)))
        .time.temp <- object$times
      }

      
    }
  
  }
  # CS Sortir la partie nouveau temps qui est indep de newdata
  if(!is.null(newtimes)){
    .pred.temp <- cbind(rep(1, dim(.pred.temp )[1]), .pred.temp)
    .time.temp <- c(0, .time.temp)
    
    # CS : MAJ du calcul pour eviter la boucle
    idx=findInterval(newtimes,.time.temp)
    .pred=.pred.temp[,pmax(1,idx)]
    # .pred <- matrix(-99, nrow = dim(.pred.temp)[1], ncol = length(newtimes))
    # .pred[,1] <- matrix(.pred.temp[ ,.time.temp <= newtimes[1]], ncol= sum(.time.temp<=newtimes[1]) )[,sum(.time.temp<=newtimes[1])]
    # for (i in 1:length(newtimes)) {
    #   .pred[,i] <- .pred.temp[,.time.temp<=newtimes[i]][,sum(.time.temp<=newtimes[i])]
    # }
    # .pred.temp <- .pred
    .time.temp <- newtimes
    .pred.temp=.pred
  }
  
  
  return(list(times=.time.temp, predictions=.pred.temp))
}
