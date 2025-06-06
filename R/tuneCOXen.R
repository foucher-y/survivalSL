tuneCOXen<- function(formula, data, penalty = NULL, cv = 10, parallel =
                       FALSE, alpha=seq(.1,.9,.1), lambda=NULL, seed = NULL){




  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }


  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(alpha)) stop("The 'alpha' argument is required.")
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

  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")


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
  .results<-c()
  
  if(!(is.null(penalty))){
    
    if(length(penalty)!=length(variables_formula[-c(1,2)]))stop("Penalty length does not equal the number of variables.")
    if(!all(unique(penalty) %in% c(0,1)))stop("Penalty must be numeric and have only 0 or 1.")}
  

  set.seed(seed)
  foldid <- sample(rep(seq(cv), length.out = nrow(.x)))

  if(!(is.null(penalty))) {
    for( a in 1:length(alpha)){
      .cv.en<-glmnet::cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance", foldid=foldid,
                                foldsid="folds", parallel = parallel, alpha=alpha[a],
                                penalty.factor = penalty,
                                lambda=lambda)
      .results<-rbind(.results,
                      cbind(rep(alpha[a],length(.cv.en$lambda)),.cv.en$lambda,.cv.en$cvm))

    }
  }

  else{
    for(a in 1:length(alpha)){
      .cv.en<-glmnet::cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance", foldid=foldid,
                                foldsid="folds", parallel = parallel, alpha=alpha[a],
                                lambda=lambda)
      .results<-rbind(.results,
                      cbind(rep(alpha[a],length(.cv.en$lambda)),.cv.en$lambda,.cv.en$cvm))
    }
  }



  colnames(.results)=c("alpha","lambda","cvm")
  .results=data.frame(.results)


  return(list(optimal=list(alpha=.results[which(.results$cvm==min(.results$cvm)),1] ,
                           lambda=.results[which(.results$cvm==min(.results$cvm)),2] ),
              results = .results))
}



