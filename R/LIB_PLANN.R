
LIB_PLANN <- function(times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL,
                    data, inter, size, decay, maxit, MaxNWts){

  .outcome <- paste("Surv(", times, ",", failures, ")")

  if(!(is.null(group))) {
    if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", group, "+", paste( cov.quanti,  collapse = " + "),
                              " + ", paste(cov.quali, collapse = " + "),
                              collapse = " ") )
    }
    if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", group, "+",
                              paste( cov.quanti, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", group, "+",
                              paste(cov.quali, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", group) )
    }
  }   else {
    if(is.null(cov.quanti)==F & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~", paste( cov.quanti,  collapse = " + "), " + ", paste(cov.quali, collapse = " + "),
                              collapse = " ") )
    }
    if(is.null(cov.quanti)==F & is.null(cov.quali)==T) {
      .f <- as.formula( paste(.outcome, "~", paste( cov.quanti, collapse = " + "),collapse = " ") )
    }
    if(is.null(cov.quanti)==T & is.null(cov.quali)==F) {
      .f <- as.formula( paste(.outcome, "~",  paste(cov.quali, collapse = " + "),collapse = " ") )
    }
  }

  .plann <- survivalPLANN(formula=.f, data=data, inter, size = size, decay = decay,  maxit = maxit, MaxNWts = MaxNWts)

  .time <- sort(unique(data[,times]))

  .pred.plann <- predict(.plann)

  .survival <- cbind(rep(1, dim(.pred.plann$predictions)[1]), .pred.plann$predictions)
  .time.interest <- c(0, .pred.plann$times)

  idx=findInterval(.time,.time.interest)
  .pred=.survival[,pmax(1,idx)]

  .obj <- list(model=.plann,
               library="LIB_PLANN",
               group=group, cov.quanti=cov.quanti, cov.quali=cov.quali,
               data=data.frame(times=data[,times], failures=data[,failures],
                               data[, !(dimnames(data)[[2]] %in% c(times, failures))]),
               times=.time,  predictions=.pred)

  class(.obj) <- "libsl"

  return(.obj)
}
