LIB_RSF <- function(formula,
                    data, nodesize, mtry, ntree, seed=NULL){

  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(nodesize)) stop("The 'nodesize' argument is required.")
  if (missing(mtry)) stop("The 'mtry' argument is required.")
  if (missing(ntree)) stop("The 'ntree' argument is required.")

  variables_formula <- all.vars(formula)

  times <- variables_formula[1]
  failures <- variables_formula[2]


  if("." %in% variables_formula){
    vars<-setdiff(names(data),c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
    variables_formula <- all.vars(formula)
  }


  variables_existent <- all(variables_formula %in% names(data))
  if (!variables_existent) stop("One or more variables from the formula do not exist in the data.")

  rm(variables_existent)

  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  ns_terms <- grep("ns\\(", all_terms, value = TRUE)
  bs_terms <- grep("bs\\(", all_terms, value = TRUE)
  if(length(strata_terms) >= 1) stop("The 'randomForestSRC' package does not support the use of 'strata()' in the formula.")


  if((length(ns_terms) >= 1)|(length(bs_terms) >= 1)|(length(all_terms)>(length(variables_formula)-2))){
    vars<-setdiff(variables_formula,c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
  }
  rm(all_terms,strata_terms,ns_terms,bs_terms)



  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")


  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }



  options(rf.cores=1, mc.cores=1)
  .rf <- rfsrc(formula, data = data, nodesize = nodesize, mtry = mtry, ntree = ntree, splitrule="logrank",seed=-seed)

  .time <- sort(unique(data[[times]]))

  .pred.rf <- predict(.rf,)
  .survival <- cbind(rep(1, dim(.pred.rf$survival.oob)[1]), .pred.rf$survival.oob) # We add a 1 because the function is step-like, so it's possible that an element
  #of .time is between 0 and x (time of the first event), and in that case, we assign it a survival of 1.
  .time.interest <- c(0, .pred.rf$time.interest)

  .idx=findInterval(.time,.time.interest)
  # indInterval: it returns which interval the values of .time belong to within .time.interest.
  .pred=cbind(rep(1, dim(.survival[,.idx])[1]),.survival[,.idx] )


  .obj <- list(model=.rf,
               library="LIB_RSF",
               formula=formula,
               data=data,
               times=c(0,.time),predictions=.pred)

  class(.obj) <- "libsl"

  return(.obj)
}




