tunePHspline<- function(formula, data, cv = 10, metric = "auc",k=1:4, pro.time
                        = NULL, seed = NULL, ROC.precision = seq(0.01, 0.99,
                                                                 by = 0.01)){

  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(k)) stop("The 'k' argument is required.")
  if (missing(metric)) stop("The 'metric' argument is required.")


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
  if(length(strata_terms) >= 1) stop("The 'flexsurv' package does not support the use of 'strata()' in the formula.")

  rm(all_terms,strata_terms)

  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }


  if(is.null(pro.time)){
    pro.time=median(data[[times]])
  }


  .data_bis<-data
  .time <- unique(sort(c(0,pro.time,data[[times]])))




  set.seed(seed)
  sample_id <- sample(nrow(data))
  folds <- cut(seq(1,nrow(data)), breaks=cv, labels=FALSE)
  folds_id <- folds[sample_id]
  data$folds <- folds_id
  data$id<-1:nrow(data)

  CVtune <- lapply(1:cv, function(i) {
    list(
      train = data[data$folds != i, ],
      valid = data[data$folds == i, ]
    )})




  if(any(sapply(data,is.factor))){
    factor_vars<-names(data)[sapply(data,is.factor)]
    inside<-function(factor,train,valid){
      result<-all(unique(valid[,factor]) %in% unique(train[,factor]))
      return(result)
    }
    check_CVtune<-function(factors,CV){
      result<-unlist(lapply(factors,inside,train=CV$train,valid=CV$valid))
      return(all(result))
    }

    i<-0

    while(check_CVtune(factor_vars,CVtune)==FALSE){
      i<-i+1
      seed<-sample(1:1000,1)
      warning("The seed has been changed.")
      set.seed(seed)
      sample_id <- sample(nrow(data))
      folds <- cut(seq(1,nrow(data)), breaks=cv, labels=FALSE)
      folds_id <- folds[sample_id]
      data$folds <- folds_id
      data$id<-1:nrow(data)

      CVtune <- lapply(1:cv, function(i) {
        list(
          train = data[data$folds != i, ],
          valid = data[data$folds == i, ]
        )})

      if(i>=3 & check_CVtune(factor_vars,CVtune)==FALSE)stop("Certain levels of some factor variables in the validation sample are not present in the training sample. Please change the seed.")

    }


  }


  spline_function<-function(x,knot){
    .flex<-flexsurvspline(formula, data = x$train, scale="hazard", k=knot,
                          hessian=F, method="Nelder-Mead")


    .haz=NULL

    if(metric=="ll"){
      .hazlist <- predict(
        .flex,
        newdata=x$valid,
        type = "haz",
        times = .time
      )
      .haz<-t(sapply(.hazlist$.pred, function(x) x[[2]]))
    }


    .survivalist<-predict(
      .flex,
      newdata=x$valid,
      type = "survival",
      times = .time
    )

    .survivals <- t(sapply(.survivalist$.pred, function(x) x[[2]]))
    return(predictions=list(survivals=.survivals, haz=.haz, id=x$valid$id))

  }

  knot_function<-function(knot){
    result<-lapply(CVtune,spline_function,knot=knot)
    surv_list<-lapply(result,function(x)(return(x$survivals)))
    id_list<- lapply(result,function(x)(return(x$id)))
    survivals<-do.call(rbind,surv_list)
    id<-unlist(id_list)
    haz<-NULL
    if(metric=="ll"){
      haz_list<-lapply(result,function(x)(return(x$haz)))
      haz<-do.call(rbind,haz_list)
    }

    return(predictions=list(survivals=survivals, haz=haz, id=id))

  }


  result<-lapply(k,knot_function)


  metric_function<-function(x){
    data<-data[x$id, ]
    survivals.matrix<-x$survivals
    hazards.matrix<-NULL
    if(metric=="ll"){
      hazards.matrix<-x$haz
    }
    resultat<-metrics(metric,times,failures,data,survivals.matrix,hazards.matrix,.time,pro.time,ROC.precision)
    return(resultat)

  }


  metric_results<-unlist(lapply(result,metric_function))


  if(metric %in% c("bs","ibs","ribs","bll","ibll","ribll")){
    .idx<-which.min(metric_results)
  }

  else{
    .idx<-which.max(metric_results)
  }

  data<-.data_bis


  return( list(optimal=list(k=k[.idx]), results=data.frame(k=k,metric_results=metric_results)))
}


