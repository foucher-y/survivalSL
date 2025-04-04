tuneRSF <- function(formula, data, nodesize, mtry, ntree, seed=NULL){


  if(any(sapply(data,is.character)))stop("Error : some columns are of type character. Only numeric or factor variables are allowed.")


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


  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)

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

  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }


  old <- options()
  on.exit(options(old))

  options(rf.cores=1, mc.cores=1)

  find.tune.rf.fast<-function(param.test, f, data ){

    res.rsf <- rfsrc(f, data = data, nodesize = param.test[1], mtry = param.test[2] ,
                     ntree = param.test[3] , seed=-seed, splitrule="logrank")

    res<-c(param.test[1], param.test[2], param.test[3],  tail(res.rsf$err.rate, 1))

    return(res)
  }

  .grid <-  expand.grid(nodesize=nodesize, mtry=mtry, ntree=ntree)
  .grid=cbind(.grid[,1],.grid[,2], .grid[,3])


  .tune.rf<-apply(.grid,MARGIN=1, FUN=find.tune.rf.fast, f=formula, data=data)

  .res=t(.tune.rf)
  colnames(.res)=c("nodesize","mtry","ntree","error")
  .res=data.frame(.res)

  .mini<-.res[which(.res$error==min(.res$error, na.rm=TRUE) & is.na(.res$error)==FALSE),]
  .mini<-.mini[1,]

  .optimal=list(nodesize=as.numeric(.mini$nodesize),
                mtry=as.numeric(.mini$mtry),
                ntree=as.numeric(.mini$ntree))

  return(list(optimal=.optimal, results = t(.tune.rf)))
}


