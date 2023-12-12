
nn.time <- function(times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL,
                      data, n.nodes, decay,
                      batch.size, epochs){

  if(!(is.null(group))){
    .data <- data[,c(times, failures, group, cov.quanti, cov.quali)] } else{
    .data <- data[,c(times, failures, cov.quanti, cov.quali)]   }

  .f  <- as.formula(paste("Surv(", times, ",", failures, ")", "~."))
  
  .deepsurv <- deepsurv(.f, data = .data,  verbose = FALSE, num_nodes=n.nodes, 
                        weight_decay=decay, num_workers =0L,
                        batch_size=batch.size,
                        epochs=epochs)
  
  # .time <- .time.deepsurv <- sort(unique(.data[,times]))
  .time <- sort(unique(.data[,times]))

  .pred <- predict(.deepsurv, newdata = .data)  
  .time.deepsurv<-as.numeric(dimnames(.pred)[[2]])
  
  idx=findInterval(.time, .time.deepsurv)
  .pred=.pred[,pmax(1,idx)]
  
  # .pred<-.pred[,-c(1,dim(.pred)[2])]
  
  .obj <- list(model=.deepsurv, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali,
          data=data.frame(times=data[,times], failures=data[,failures],
          data[, !(dimnames(data)[[2]] %in% c(times, failures))]), times=.time, predictions=.pred)

  class(.obj) <- "nn.time"

  return(.obj)
}
