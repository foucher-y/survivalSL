LIB_COXall<- function(formula, data){


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
  if(length(strata_terms) > 1) stop("More than one 'strata' term found in the formula. Only one stratified variable is allowed.")

  rm(all_terms,strata_terms)


  if(any(sapply(data,is.character)))stop("Error : some columns are of type character. Only numeric or factor variables are allowed.")


  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }


  .coxph <- coxph(formula, data=data)


  .coxphsurv<-survfit(.coxph, newdata = data,se.fit = F)


  .sumcoxphsurv<-summary(.coxphsurv, times=sort(unique(data[[times]])))
  .pred <- t(.sumcoxphsurv$surv)
  .survivals<-cbind(rep(1, dim(.pred)[1]), .pred)


  .obj <- list(model=.coxph,
               library="LIB_COXall",
               formula=formula,
               data=data,
               times=c(0,sort(unique(data[[times]]))),predictions=.survivals )

  class(.obj) <- "libsl"

  return(.obj)
}





