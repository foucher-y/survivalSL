

plot.sl.time <- function(x, ..., col=1, lty=1, lwd=1, type="b", pch = 16, ylab=NULL, xlab=NULL,
                         ylim=NULL, xlim=NULL, cex=1, cex.lab=1, cex.axis=1, cex.main=1,
                         method="sl", n.groups=5, pro.time=NULL, newdata=NULL,times=NULL, failures=NULL)
{
  pred.times <- x$times
  
  if(is.null(newdata))
  {
  times <-  x$data$times
  failures <- x$data$failures

  if(method=="sl" & is.list(x$predictions)){
    pred.matrix<-x$predictions[[length(x$predictions)]]
  }
  if(method=="sl" & !is.list(x$predictions)){
    pred.matrix<-x$predictions
  }
  if(method!="sl"){
    pred.matrix <- x$models[[method]]$predictions
  }
  
  
  if(is.null(pro.time)) {pro.time <- median(times)}
  
  .pred <- pred.matrix[,pred.times<pro.time][,sum(pred.times<pro.time)]
  
  .grps <- as.numeric(cut(.pred, 
                          breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))), 
                          labels = 1:n.groups))
  
  #.grps <- as.numeric(quant_groups(.pred, groups = n.groups))
  
  .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
  
  .survfit <- summary(survfit(Surv(times, failures) ~ as.factor(.grps)))
  
  .obs <- sapply(1:n.groups, FUN = function(x) { last(.survfit$surv[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
  .lower <- sapply(1:n.groups, FUN = function(x) { last(.survfit$lower[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
  .upper <- sapply(1:n.groups, FUN = function(x) { last(.survfit$upper[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
  
  if(is.null(ylab)) {ylab <- "Observed survival"}
  if(is.null(xlab)) {xlab <- "Predicted survival"}
  
  if(is.null(ylim)) {ylim <- c(0,1)}
  if(is.null(xlim)) {xlim  <- c(0,1)}
  
  plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
       type = type, col = col, lty = lty, lwd = lwd,
       pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
  
  abline(c(0,1), lty=2)
  
  segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
  }
  
  else
  {
    .t <- newdata[,times]
    .f <- newdata[,failures]
    
    if(is.null(pro.time)) {pro.time <- median(.t)}
    
    pred.valid <- predict(x, newdata=newdata, newtimes=pred.times)
    
    prediction.matrix <- pred.valid$predictions[[method]]
    
    # pred.times<-c(0,pred.times)
    .pred <- prediction.matrix[,pred.times<pro.time][,sum(pred.times<pro.time)]
    # .pred <- prediction.matrix[,pred.valid$times<pro.time][,sum(pred.valid$times<pro.time)]
    
    .grps <- as.numeric(cut(.pred, 
               breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))), 
               labels = 1:n.groups))
    
    #.grps <- as.numeric(dvmisc::quant_groups(.pred, groups = n.groups))
    
    .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
    
    .survfit <- summary(survfit(Surv(.t, .f) ~ as.factor(.grps)))
    
    .obs <- sapply(1:n.groups, FUN = function(x) { last(.survfit$surv[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .lower <- sapply(1:n.groups, FUN = function(x) { last(.survfit$lower[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .upper <- sapply(1:n.groups, FUN = function(x) { last(.survfit$upper[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    
    if(is.null(ylab)) {ylab <- "Observed survival"}
    if(is.null(xlab)) {xlab <- "Predicted survival"}
    
    if(is.null(ylim)) {ylim <- c(0,1)}
    if(is.null(xlim)) {xlim  <- c(0,1)}
    
    plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         type = type, col = col, lty = lty, lwd = lwd, pch = pch,
         ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
    
    abline(c(0,1), lty=2)
    
    segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
  }
  
  
}






