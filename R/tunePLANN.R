tunePLANN <- function(formula, data, cv=10, inter=1, size=c(2, 4, 6, 8, 10), decay=c(0.001, 0.01, 0.02, 0.05), maxit=100, MaxNWts=10000,maxtime=NULL,seed=NULL,metric="auc",pro.time=NULL,ROC.precision=seq(.01, .99, by=.01)){


  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(inter)) stop("The 'inter' argument is required.")
  if (missing(decay)) stop("The 'decay' argument is required.")
  if (missing(maxit)) stop("The 'maxit' argument is required.")
  if (missing(MaxNWts)) stop("The 'MaxNWts' argument is required.")


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
  ns_terms <- grep("ns\\(", all_terms, value = TRUE)
  bs_terms <- grep("bs\\(", all_terms, value = TRUE)

  if(length(strata_terms) >= 1) stop("The 'survivalPLANN' package does not support the use of 'strata()' in the formula.")


  if((length(ns_terms) >= 1)|(length(bs_terms) >= 1)|(length(all_terms)>(length(variables_formula)-2))){
    vars<-setdiff(variables_formula,c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
  }

  rm(all_terms,strata_terms,ns_terms,bs_terms)



  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }



  if(metric=="ll"){
    haz_function<-function(surv,times){
      x<-1
      result<-sapply(2:length(surv),function(i){
        value<-surv[i]
        if(value!=surv[i-1]){
          x<<- x+1
        }
        return(x)
      })


      df<-data.frame(temps=times,value=-log(surv),result=c(1,result))
      df_unique <- df[!duplicated(df$result), ]
      if(nrow(df_unique)>1){
        diff_1<-diff(df_unique$temps)
        diff_2 <- diff(df_unique$value)
        resultat <- diff_2/diff_1
        resultat<-c(Inf,resultat,NA)
        idx=findInterval(times,c(0,df_unique$temps))
        bj<-c(Inf,resultat[idx])
      }
      else{
        bj<-c(rep(Inf,(length(surv)-1)),NA)
      }


      return(bj)

    }


  }


  if(is.null(pro.time)){
    pro.time=median(data[[times]])
  }


  .data_bis<-data
  .time <-unique(sort(c(0,pro.time,data[[times]])))

  if (is.null(maxtime) || maxtime < max(.time)) {
    maxtime <- max(.time) + 1
  }



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





  hp<-expand.grid(inter, decay, size, maxit,MaxNWts)
  colnames(hp)<-c("inter", "decay", "size", "maxit","MaxNWts")
  hp_list<-apply(hp, 1, as.list)

  plann_function<-function(x,y){

    model_matrix <- model.matrix(formula, data = rbind(x$train,x$valid))[1:nrow(x$train),-1]

    data_bis<-data.frame(times=x$train[[times]],failures=x$train[[failures]])
    colnames(data_bis)<-c(times,failures)

    vars<-setdiff(names(cbind(data_bis,model_matrix)),c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    .formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))

    .plann <- sPLANN(formula=.formula, data=cbind(data_bis,model_matrix), inter=y$inter, size = y$size, decay = y$decay,  maxit = y$maxit, MaxNWts = y$MaxNWts, pro.time=maxtime)

    newdata <- model.matrix(formula, rbind(x$valid,x$train))[1:nrow(x$valid),]
    .survivals <- predict(.plann, newdata = data.frame(newdata), newtimes = .time)$predictions
    return(predictions=list(survivals=.survivals, id=x$valid$id))

  }



  y_function<-function(y){
    result<-lapply(CVtune,plann_function,y=y)
    surv_list<-lapply(result,function(x)(return(x$survivals)))
    id_list<- lapply(result,function(x)(return(x$id)))
    survivals<-do.call(rbind,surv_list)
    id<-unlist(id_list)
    return(predictions=list(survivals=survivals, id=id))

  }



  result<-lapply(hp_list,y_function)

  metric_function<-function(x){
    data<-data[x$id, ]
    survivals.matrix<-x$survivals
    hazards.matrix<-NULL
    if(metric=="ll"){
      hazards.matrix<-t(apply(x$survivals[,-1],1,haz_function,times=.time))
    }
    resultat<-metrics(metric=metric,times=times,failures=failures,data=data,survivals.matrix=survivals.matrix,hazards.matrix=hazards.matrix,prediction.times=.time,pro.time=pro.time,ROC.precision=ROC.precision)
    return(resultat)

  }


  metric_results<-unlist(lapply(result,metric_function))


  if(metric %in% c("bs","ibs","ribs","bll","ibll","ribll")){
    .idx<-which.min(metric_results)
  }else{
    .idx<-which.max(metric_results)
  }

  data<-.data_bis


  return( list(optimal=list(
    inter=hp[.idx,1],
    decay=hp[.idx,2],
    size=hp[.idx,3],
    maxit=hp[.idx,4],
    MaxNWts=hp[.idx,5]

  ), results=cbind(hp,metric_results) ))

}





