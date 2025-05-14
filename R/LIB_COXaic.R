LIB_COXaic<- function(formula, data, penalty=NULL){

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
  if(length(strata_terms) >= 1) stop("More than one 'strata' term found in the formula. No stratified variable is allowed.")
  
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
  
  
  variables_formula_bis<- variables_formula[-c(1,2)]

  
  if(!(is.null(penalty))){
    
    if(length(penalty)!=length(variables_formula_bis))stop("Penalty length does not equal the number of variables.")
    if(!all(unique(penalty) %in% c(0,1)))stop("Penalty must be numeric and have only 0 or 1.")
    var<- variables_formula_bis[ifelse(penalty==0,TRUE,FALSE)]}


  .outcome <- paste("Surv(", times, ",", failures, ")")

  if(!(is.null(penalty))){ # Case one : we want absolutely some variables in our model !
    .f0<- as.formula(paste(.outcome, "~", paste(var, collapse = " + ")))
  }
  else{ # Else : we don't care
    .f0<-as.formula( paste(.outcome, "~1" ))
  }



  .fit0<-coxph(.f0, data=data)

  .fit<-coxph(formula, data=data)

  .final.formula<-stepAIC(.fit0, scope=formula(.fit), direction="forward", k=2, trace=FALSE)$formula

  .result<-LIB_COXall(.final.formula,data)


  .obj <- list(model=.result$model,
               library="LIB_COXaic",
               formula=.result$formula,
               data=.result$data,
               times=.result$times,predictions=.result$predictions )





  class(.obj) <- "libsl"

  return(.obj)
}

