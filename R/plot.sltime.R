plot.sltime <- function(x, method="sl", n.groups=5, pro.time=NULL, newdata=NULL, ...)
{

  variables_formula <- all.vars(x$formula)
  .times <- variables_formula[1]
  .failures <- variables_formula[2]

  if(is.null(newdata))
  {

    pred.times <- x$times

    if(is.null(pro.time)) {pro.time <- median(pred.times)}
    times <- ifelse(x$data[[.times]]>=pro.time,pro.time,x$data[[.times]])
    failures <- ifelse(x$data[[.failures]]==1 & x$data[[.times]]>pro.time,0,x$data[[.failures]])


    if(method!="sl"){
      idx<-which(x$methods==method)
      object<-x$models[[idx]]
      .pred<-predict(object=object,newtimes=pro.time)$predictions
    }
    else{
      idx<-length(x$methods)+1
      .pred<-predict(object=x,newtimes=pro.time)$predictions[[idx]]
    }



    .grps <- as.numeric(cut(.pred,
                            breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))

    .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )

    .survfit <- summary(survfit(Surv(times, failures) ~ as.factor(.grps)))

    .obs <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$surv[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .lower <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$lower[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .upper <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$upper[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )

    if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
    if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
    if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
    if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
    if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
    if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
    if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
    if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
    if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}

    if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
    if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}

    if(hasArg(ylab)==FALSE) {ylab <- "Mean of survival predictions"} 
    else {ylab <- list(...)$ylab}
    if(hasArg(xlab)==FALSE) {xlab <- "Kaplan-Meier estimations"} else 
    {xlab <- list(...)$xlab}
    
    plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         type = type, col = col, lty = lty, lwd = lwd,
         pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)

    abline(c(0,1), lty=2)

    segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
    text(x = 0.2, y = 0.02, labels = paste("Calculated at time =", 
                                           round(pro.time, 2)), pos = 4, cex = 1)
  }

  else
  {

    .t <- newdata[[.times]]
    .f <- newdata[[.failures]]

    if(is.null(pro.time)) {pro.time <- median(.t)}

    if(method!="sl"){
      idx<-which(x$methods==method)
      object<-x$models[[idx]]
      .pred<-predict(object=object,newdata=newdata,newtimes=pro.time)$predictions
    }
    else{
      idx<-length(x$methods)+1
      .pred<-predict(object=x,newdata=newdata,newtimes=pro.time)$predictions[[idx]]
    }


    times <- ifelse(.t>=pro.time,pro.time,.t)
    failures <- ifelse(.f==1 & .t>pro.time,0,.f)

    .grps <- as.numeric(cut(.pred,
                            breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))

    .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )

    .survfit <- summary(survfit(Surv(times, failures) ~ as.factor(.grps)))

    .obs <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$surv[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .lower <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$lower[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )
    .upper <- sapply(1:n.groups, FUN = function(x) {
      last(.survfit$upper[ as.numeric(.survfit$strata)==x & .survfit$time<=pro.time ]) } )

    if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
    if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
    if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
    if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
    if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
    if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
    if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
    if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
    if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}

    if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
    if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}

    if(hasArg(ylab)==FALSE) {ylab <- "Mean of survival predictions"} 
    else {ylab <- list(...)$ylab}
    if(hasArg(xlab)==FALSE) {xlab <- "Kaplan-Meier estimations"} else 
    {xlab <- list(...)$xlab}
    plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         type = type, col = col, lty = lty, lwd = lwd, pch = pch,
         ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)

    abline(c(0,1), lty=2)

    segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
    text(x = 0.2, y = 0.02, labels = paste("Calculated at time =", 
                                           round(pro.time, 2)), pos = 4, cex = 1)
  }
}


