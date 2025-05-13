tuneCOXridge <- function(formula, data, penalty = NULL, cv = 10, parallel =
                           FALSE, lambda, seed = NULL){

  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(lambda)) stop("The 'lambda' argument is required.")



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

  if(any(sapply(data[,variables_formula],is.character)))stop("Error : some columns are of type character. Only numeric or factor variables are allowed.")


  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms) >= 1) stop("The 'glmnet' package does not support the use of 'strata()' in the formula.")

  rm(all_terms,strata_terms)

  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)

  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }

  .y <- Surv(data[[times]], data[[failures]])
  .x <- model.matrix(formula,data)[,-1]
  set.seed(seed)
  foldid <- sample(rep(seq(cv), length.out = nrow(.x)))


  if(!(is.null(penalty))) {
    .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = parallel, alpha=0,keep=F, foldid = foldid,
                           penalty.factor=penalty,
                           lambda=lambda
    )

  }

  else{
    .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = parallel, alpha=0,keep=F, foldid = foldid,
                           lambda=lambda
    )
  }



  return(list(optimal=list(lambda=.cv.ridge$lambda.min),
              results = data.frame( lambda=.cv.ridge$lambda, deviance=.cv.ridge$cvm)))
}

