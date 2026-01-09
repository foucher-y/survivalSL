LIB_PHgompertz <- function(formula,data){

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")

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
  if(length(strata_terms) >= 1) stop("The 'flexsurv' package does not support the use of 'strata()' in the formula.")

  rm(all_terms,strata_terms)


  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")

  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }



  .flex<-flexsurvreg(formula, data = data,
                     dist = "gompertz",
                     inits=c(-1,1/mean(data[[times]])),
                     hessian=F, method="Nelder-Mead")


  .survivalist<-predict(
    .flex,
    newdata=data,
    type = "survival",
    times = unique(sort(c(0,data[[times]])))
  )




  .survivals <- t(sapply(.survivalist$.pred, function(x) x[[2]]))


  .obj <- list(model=.flex,
               library="LIB_PHgompertz",
               formula=formula,
               data=data,
               times=unique(sort(c(0,data[[times]]))), predictions=.survivals)

  class(.obj) <- "libsl"

  return(.obj)
}




