
survivalSL <- function(methods, metric="ci",  data, times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL,
                     cv=10, param.tune=NULL, pro.time=NULL,  optim.local.min=FALSE,
                     ROC.precision=seq(.01,.99,.01), param.weights.fix=NULL,
                     param.weights.init=NULL, keep.predictions=TRUE,
                     progress=TRUE) {

  #####################
  ### Quality tests ###
  #####################

  if(length(methods)<=1)
  { stop("Number of methods need to be greater or equal to 2 to estimate a SuperLearner")   }

  if(length(metric)>1){
    warning(paste0("SuperLearner is currently developped for one metric. Results for metric ",metric[1]))
    metric=metric[1]
  }

  if(min(metric%in%c("ci","bs","loglik","ibs","ibll","bll", "ribs","ribll","auc"))==0){
    stop("The argument \"metric\" must be Brier score (bs),
         Concordance index (ci),
         Integrated Brier score (ibs), the binomilar log-likelihood (bll),
         the Integrated binomial log-likelihood (ibll), the restricted ibs (ribs),
         the restricted ibll (ribll), the log-likelihood (loglik), or
         the area under the ROC curve (auc)")
  }

  if(!is.data.frame(data) & !is.matrix(data)){
    stop("The argument \"data\" need to be a data.frame or a matrix") }


  if( is.null(group)==F){
    if(length(group)>1){
      stop("Only one variable can be use as group")
    }
    if(min(group %in%colnames(data))==0 & is.character(group)==T){
      stop("Group name is not present in data")
    }
  }

  if( is.null(cov.quanti)==F){
    if(min(group %in%colnames(data))==0 & is.character(cov.quanti)==T){
      stop("At least one name of quantitative covariate is not present in data")
    }
  }

  if( is.null(cov.quali)==F){
    if(min(cov.quali %in%colnames(data))==0 & is.character(cov.quali)==T){
      stop("At least one name of qualitative covariate is not present in data")
    }
  }

  if(is.null(group)==T&is.null(cov.quanti)==T&is.null(cov.quali)==T){
    stop("SuperLearner need at least one group or one quantitative or one qualitative covariate")
  }

  if(!is.null(group)==T){
    data<-data[,c(times,failures,group,cov.quanti,cov.quali)]
  }
  if(is.null(group)==T){
    data<-data[,c(times,failures,cov.quanti,cov.quali)]
  }

  if (any(is.na(data))){
    data<-na.omit(data)
    warning("Data need to be without NA. NA is removed")
  }

  if(!(is.null(group))){
    if(!is.character(group) & !is.numeric(group) ){
      stop("The argument \"group\" need to be scalar or a character string") }

    mod <- unique(data[,group])
    if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
      stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the argument \"group\" ")
    }

  }

  if(length(data[,times])!=length(data[,failures])){
    stop("The length of the times must be equaled to the length of the events in the training data") }

  mod2 <- unique(data[,failures])
  if(length(mod2) != 2 | ((mod2[1] != 0 & mod2[2] != 1) & (mod2[1] != 1 & mod2[2] != 0))){
    stop("Two modalities encoded 0 (for censored patients) and 1 (for dead patients) are required in the argument \"failures\" ")
  }

  if (!is.numeric(data[,times])){
    stop("Time variable is not numeric")}

  if (min(data[,times])<=0){
    stop("Time variable need to be positive")
  }

  if (cv < 3 | !is.numeric(cv)) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  .meth_rm=c()
  if(sum(methods %in% "AFTgamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="AFTgamma")[-1])
    warning("SuperLearner can use only one AFTgamma method. We remove the others.")
  }
  if(sum(methods %in% "AFTllogis")>=2){
    .meth_rm=c(.meth_rm,which(methods=="AFTllogis")[-1])
    warning("SuperLearner can use only one AFTllogis method. We remove the others.")
  }
  if(sum(methods %in% "AFTggamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="AFTggamma")[-1])
    warning("SuperLearner can use only one AFTggamma method. We remove the others.")
  }
  if(sum(methods %in% "AFTweibull")>=2){
    .meth_rm=c(.meth_rm,which(methods=="AFTweibull")[-1])
    warning("SuperLearner can use only one AFTweibull method. We remove the others.")
  }
  if(sum(methods %in% "PHexponential")>=2){
    .meth_rm=c(.meth_rm,which(methods=="PHexponential")[-1])
    warning("SuperLearner can use only one PHexponential method. We remove the others.")
  }
  if(sum(methods %in% "PHgompertz")>=2){
    .meth_rm=c(.meth_rm,which(methods=="PHgompertz")[-1])
    warning("SuperLearner can use only one PHgompertz method. We remove the others.")
  }
  if(sum(methods %in% "PHspline")>=2){
    .meth_rm=c(.meth_rm,which(methods=="PHspline")[-1])
    warning("SuperLearner can use only one PHspline method. We remove the others.")
  }

  if(sum(methods %in% "COXlasso")==1){
    if(!(is.null(param.tune[[which(methods=="COXlasso")]]))){
      if(!is.list(param.tune[[which(methods=="COXlasso")]])){
        stop("Argument param.tune for COXlasso need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="COXlasso")]])%in%"lambda"))==0){
        stop("Tune parameters for COXlasso need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="COXlasso")]]$lambda)|
           is.null(param.tune[[which(methods=="COXlasso")]]$lambda))){
        stop("Lambda tune parameters for COXlasso need to be a scalar or a vector or NULL")
      }
    }
  }

  if(sum(methods %in% "COXlasso")>=2){
    if(length(param.tune[which(methods=="COXlasso")])!=length(unique(param.tune[which(methods=="COXlasso")]))){
      stop("Tune parameters for COXlasso methods need to be unique")
    }
    for (i in 1:sum(methods %in% "COXlasso")){
      if(!(is.null(param.tune[[which(methods=="COXlasso")[i]]]))){
        if(!is.list(param.tune[[which(methods=="COXlasso")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th COXlasso need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="COXlasso")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th COXlasso need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="COXlasso")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="COXlasso")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th COXlasso need to be a scalar or a vector or NULL"))
        }
      }
    }
  }



  if(sum(methods %in% "PHspline")==1){
    if(!(is.null(param.tune[[which(methods=="PHspline")]]))){
      if(!is.list(param.tune[[which(methods=="PHspline")]])){
        stop("Argument param.tune for PHspline need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="PHspline")]])%in%"k"))==0){
        stop("Tune parameters for PHspline need to have k")
      }
      if(!(is.numeric(param.tune[[which(methods=="PHspline")]]$k)|
           is.null(param.tune[[which(methods=="PHspline")]]$k))){
        stop("Lambda tune parameters for PHspline need to be a scalar or a vector or NULL")
      }
    }
  }

  if(sum(methods %in% "PHspline")>=2){
    if(length(param.tune[which(methods=="PHspline")])!=length(unique(param.tune[which(methods=="PHspline")]))){
      stop("Tune parameters for PHspline methods need to be unique")
    }
    for (i in 1:sum(methods %in% "PHspline")){
      if(!(is.null(param.tune[[which(methods=="PHspline")[i]]]))){
        if(!is.list(param.tune[[which(methods=="PHspline")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th PHspline need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="PHspline")[i]]])%in%"k"))==0){
          stop(paste("Tune parameters for the ",i,"th PHspline need to have k"))
        }
        if(!(is.numeric(param.tune[[which(methods=="PHspline")[i]]]$k)|
             is.null(param.tune[[which(methods=="PHspline")[i]]]$k))){
          stop(paste("Lambda tune parameters for the ",i,"th PHspline need to be a scalar or a vector or NULL"))
        }
      }
    }
  }


  if(sum(methods %in% "COXridge")==1){
    if(!(is.null(param.tune[[which(methods=="COXridge")]]))){
      if(!is.list(param.tune[[which(methods=="COXridge")]])){
        stop("Argument param.tune for COXridge need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="COXridge")]])%in%"lambda"))==0){
        stop("Tune parameters for COXridge need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="COXridge")]]$lambda)|
           is.null(param.tune[[which(methods=="COXridge")]]$lambda))){
        stop("Lambda tune parameters for COXridge need to be a scalar or a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "COXridge")>=2){
    if(length(param.tune[which(methods=="COXridge")])!=length(unique(param.tune[which(methods=="COXridge")]))){
      stop("Tune parameters for COXridge methods need to be unique")
    }
    for (i in 1:sum(methods %in% "COXridge")){
      if(!(is.null(param.tune[[which(methods=="COXridge")[i]]]))){
        if(!is.list(param.tune[[which(methods=="COXridge")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th COXridge need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="COXridge")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th COXridge need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="COXridge")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="COXridge")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th COXridge need to be a scalar or a vector or NULL"))
        }
      }
    }
  }

  if(sum(methods %in% "COXen")==1){
    if(!(is.null(param.tune[[which(methods=="COXen")]]))){
      if(!is.list(param.tune[[which(methods=="COXen")]])){
        stop("Argument param.tune for COXen need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="COXen")]])%in%"lambda"))==0){
        stop("Tune parameters for COXen need to have lambda")
      }
      if(sum((names(param.tune[[which(methods=="COXen")]])%in%"alpha"))==0){
        stop("Tune parameters for COXen need to have alpha")
      }
      if(!(is.numeric(param.tune[[which(methods=="COXen")]]$lambda)|
           is.null(param.tune[[which(methods=="COXen")]]$lambda))){
        stop("Lambda tune parameters for COXen need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="COXen")]]$alpha)|
           is.null(param.tune[[which(methods=="COXen")]]$alpha))){
        stop("alpha tune parameters for COXen need to be a scalar or a vector or NULL")
      }
      if(min(param.tune[[which(methods=="COXen")]]$alpha)<0 | max(param.tune[[which(methods=="COXen")]]$alpha)>1){
        stop("tune parameters for COXen alpha need to be in ]0;1[")
      }
    }
  }
  if(sum(methods %in% "COXen")>=2){
    if(length(param.tune[which(methods=="COXen")])!=length(unique(param.tune[which(methods=="COXen")]))){
      stop("Tune parameters for COXen methods need to be unique")
    }
    for (i in 1:sum(methods %in% "COXen")){
      if(!(is.null(param.tune[[which(methods=="COXen")[i]]]))){
        if(!is.list(param.tune[[which(methods=="COXen")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th COXen need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="COXen")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th COXen need to have lambda"))
        }
        if(sum((names(param.tune[[which(methods=="COXen")[i]]])%in%"alpha"))==0){
          stop(paste("Tune parameters for the ",i,"th COXen need to have alpha"))
        }
        if(!(is.numeric(param.tune[[which(methods=="COXen")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="COXen")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th COXen need to be a scalar or a vector or NULL"))
        }
        if(!(is.numeric(param.tune[[which(methods=="COXen")[i]]]$alpha)|
             is.null(param.tune[[which(methods=="COXen")[i]]]$alpha))){
          stop(paste("Alpha tune parameters for the ",i,"th COXen need to be a scalar or a vector or NULL"))
        }
        if(min(param.tune[[which(methods=="COXen")[i]]]$alpha)<0 | max(param.tune[[which(methods=="COXen")[i]]]$alpha)>1){
          stop("tune parameters for COXen alpha need to be in ]0;1[")
        }
      }
    }
  }

  if(sum(methods %in% "COXaic")==1){
    if(!(is.null(param.tune[[which(methods=="COXaic")]]))){
      if(!is.list(param.tune[[which(methods=="COXaic")]])){
        stop("Argument param.tune for COXaic need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="COXaic")]])%in%"final.model"))==0){
        stop("Tune parameters for COXaic need to have final.model")
      }
      if(sum((names(param.tune[[which(methods=="COXaic")]])%in%"model.min"))==0){
        stop("Tune parameters for COXaic need to have model.min")
      }
      if(sum((names(param.tune[[which(methods=="COXaic")]])%in%"model.max"))==0){
        stop("Tune parameters for COXaic need to have model.max")
      }
    }
  }
  if(sum(methods %in% "COXaic")>=2){
    if(length(param.tune[which(methods=="COXaic")])!=length(unique(param.tune[which(methods=="COXaic")]))){
      stop("Tune parameters for COXaic methods need to be unique")
    }
    for (i in 1:sum(methods %in% "COXaic")){
      if(!(is.null(param.tune[[which(methods=="COXaic")[i]]]))){
        if(!is.list(param.tune[[which(methods=="COXaic")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th COXaic need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="COXaic")[i]]])%in%"finl.model.cov"))==0){
          stop(paste("Tune parameters for the ",i,"th COXaic need to have finl.model.cov"))
        }
        if(sum((names(param.tune[[which(methods=="COXaic")[i]]])%in%"model.min"))==0){
          stop(paste("Tune parameters for the ",i,"th COXaic need to have model.min"))
        }
        if(sum((names(param.tune[[which(methods=="COXaic")[i]]])%in%"model.max"))==0){
          stop("Tune parameters for COXaic need to have model.max")
        }
      }
    }
  }


  if(sum(methods %in% "RSF")==1){
    if(!(is.null(param.tune[[which(methods=="RSF")]]))){
      if(!is.list(param.tune[[which(methods=="RSF")]])){
        stop("Argument param.tune for RSF need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="RSF")]])%in%"nodesize"))==0){
        stop("Tune parameters for RSF need to have nodesize")
      }
      if(sum((names(param.tune[[which(methods=="RSF")]])%in%"mtry"))==0){
        stop("Tune parameters for RSF need to have mtry")
      }
      if(sum((names(param.tune[[which(methods=="RSF")]])%in%"ntree"))==0){
        stop("Tune parameters for RSF need to have ntree")
      }
      if(!(is.numeric(param.tune[[which(methods=="RSF")]]$nodesize)|
           is.null(param.tune[[which(methods=="RSF")]]$nodesize))){
        stop("nodesize tune parameters for RSF need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="RSF")]]$mtry)|
           is.null(param.tune[[which(methods=="RSF")]]$mtry))){
        stop("mtry tune parameters for RSF need to be a scalar or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="RSF")]]$ntree)|
           is.null(param.tune[[which(methods=="RSF")]]$ntree))){
        stop("ntree tune parameters for RSF need to be a scalar or NULL")
      }
    }
  }
  if(sum(methods %in% "RSF")>=2){
    if(length(param.tune[which(methods=="RSF")])!=length(unique(param.tune[which(methods=="RSF")]))){
      stop("Tune parameters for RSF methods need to be unique")
    }
    for (i in 1:sum(methods %in% "RSF")){
      if(!is.list(param.tune[[which(methods=="RSF")[i]]])){
        stop(paste("Argument param.tune for the ",i,"th RSF need to be a list"))
      }
      if(!(is.null(param.tune[[which(methods=="RSF")[i]]]))){
        if(sum((names(param.tune[[which(methods=="RSF")[i]]])%in%"nodesize"))==0){
          stop("Tune parameters for RSF need to have nodesize")
        }
        if(sum((names(param.tune[[which(methods=="RSF")[i]]])%in%"mtry"))==0){
          stop("Tune parameters for RSF need to have mtry")
        }
        if(sum((names(param.tune[[which(methods=="RSF")[i]]])%in%"ntree"))==0){
          stop("Tune parameters for RSF need to have ntree")
        }
        if(!(is.numeric(param.tune[[which(methods=="RSF")[i]]]$nodesize)|
             is.null(param.tune[[which(methods=="RSF")[i]]]$nodesize))){
          stop("nodesize tune parameters for RSF need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="RSF")[i]]]$mtry)|
             is.null(param.tune[[which(methods=="RSF")[i]]]$mtry))){
          stop("mtry tune parameters for RSF need to be a scalar or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="RSF")[i]]]$ntree)|
             is.null(param.tune[[which(methods=="RSF")[i]]]$ntree))){
          stop("ntree tune parameters for RSF need to be a scalar or NULL")
        }
      }
    }
  }

  if(sum(methods %in% "SNN")==1){
    if(!(is.null(param.tune[[which(methods=="SNN")]]))){
      if(!is.list(param.tune[[which(methods=="SNN")]])){
        stop("Argument param.tune for SNN need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="SNN")]])%in%"n.nodes"))==0){
        stop("Tune parameters for SNN need to have n.nodes")
      }
      if(sum((names(param.tune[[which(methods=="SNN")]])%in%"decay"))==0){
        stop("Tune parameters for SNN need to have decay")
      }
      if(sum((names(param.tune[[which(methods=="SNN")]])%in%"batch.size"))==0){
        stop("Tune parameters for SNN need to have batch.size")
      }
      if(sum((names(param.tune[[which(methods=="SNN")]])%in%"epochs"))==0){
        stop("Tune parameters for SNN need to have epochs")
      }
      if(!(is.numeric(param.tune[[which(methods=="SNN")]]$n.nodes)|
           is.null(param.tune[[which(methods=="SNN")]]$n.nodes))){
        stop("n.nodes tune parameters for SNN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="SNN")]]$decay)|
           is.null(param.tune[[which(methods=="SNN")]]$decay))){
        stop("decay tune parameters for SNN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="SNN")]]$batch.size)|
           is.null(param.tune[[which(methods=="SNN")]]$batch.size))){
        stop("batch.size tune parameters for SNN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="SNN")]]$epochs)|
           is.null(param.tune[[which(methods=="SNN")]]$epochs))){
        stop("epochs tune parameters for SNN need to be a scalar, a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "SNN")>=2){
    if(length(param.tune[which(methods=="SNN")])!=length(unique(param.tune[which(methods=="SNN")]))){
      stop("Tune parameters for SNN methods need to be unique")
    }
    for (i in 1:sum(methods %in% "SNN")){
      if(!(is.null(param.tune[[which(methods=="SNN")[i]]]))){
        if(!is.list(param.tune[[which(methods=="SNN")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th SNN need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="SNN")[i]]])%in%"n.nodes"))==0){
          stop("Tune parameters for SNN need to have n.nodes")
        }
        if(sum((names(param.tune[[which(methods=="SNN")[i]]])%in%"decay"))==0){
          stop("Tune parameters for SNN need to have decay")
        }
        if(sum((names(param.tune[[which(methods=="SNN")[i]]])%in%"batch.size"))==0){
          stop("Tune parameters for SNN need to have batch.size")
        }
        if(sum((names(param.tune[[which(methods=="SNN")[i]]])%in%"epochs"))==0){
          stop("Tune parameters for SNN need to have epochs")
        }
        if(!(is.numeric(param.tune[[which(methods=="SNN")[i]]]$n.nodes)|
             is.null(param.tune[[which(methods=="SNN")[i]]]$n.nodes))){
          stop("nodesize tune parameters for SNN need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="SNN")[i]]]$decay)|
             is.null(param.tune[[which(methods=="SNN")[i]]]$decay))){
          stop("decay tune parameters for SNN need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="SNN")[i]]]$batch.size)|
             is.null(param.tune[[which(methods=="SNN")[i]]]$batch.size))){
          stop("batch.size tune parameters for SNN need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="SNN")[i]]]$epochs)|
             is.null(param.tune[[which(methods=="SNN")[i]]]$epochs))){
          stop("epochs tune parameters for SNN need to be a scalar, a vector or NULL")
        }
      }
    }
  }

  if(length(.meth_rm)>=1){
    methods=methods[-.meth_rm]
    param.tune=param.tune[-.meth_rm]
  }

  if((max(ROC.precision)==1) | (min(ROC.precision)==0)){
    stop("values for ROC.precision need to be in ]0;1[")
  }


  if(is.null(param.weights.fix)==FALSE & is.null(param.weights.init)==FALSE){
    warning("Weights can not be fix and initial at the same time. SuperLearner ignored initial values")
    param.weights.init<-NULL
  }

  if(is.null(param.weights.fix)==FALSE | is.null(param.weights.init)==FALSE){
    if(is.null(param.weights.fix)==FALSE){
      if(is.numeric(param.weights.fix)==FALSE){
        stop("param.weights.fix need to be numeric")
      }
    }

    if(is.null(param.weights.init)==FALSE){
      if(is.numeric(param.weights.init)==FALSE){
        stop("param.weights.init need to be numeric")
      }

        if(length(param.weights.init)!=(length(methods)-1)){
          stop("wrong lenth for param.weights.init")
        }

    }
  }
  if(is.null(param.weights.fix)==TRUE & is.null(param.weights.init)==TRUE){
      param.weights.init<-rep(0,length(methods)-1)
  }

  ######################################################
  ### Loss functions used for the weigths estimation ###
  ######################################################

  ibs<-function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))


    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    survs <- t(.pred)[,ot]


    bsc<-sapply(1:length(timeVector), FUN = function(j)
    {
      help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
      help2 <- as.integer(time > timeVector[j])
      return(mean((0 - survs[j, ])^2 * help1 * (1/csurv) +
                    (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j])))
    })

    idx <- 2:length(timeVector)
    RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx])/2)
    RET <- RET/diff(range(timeVector))
    RET=as.matrix(RET)

    return(RET)
  }

  brs<-function(par, FitCV, timeVector,
                obj_surv, ot, csurv, csurv_btime, time, pro.time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))

    .par<-exp(c(par))/(1+sum(exp(par)))

    .par<-c(.par,1-sum(.par))

    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    survs <- t(.pred)[,ot]


    j=length(timeVector[which(timeVector<=pro.time)])


    help1 <- as.integer(time <= timeVector[j] &  obj_surv[ot,2] == 1)
    help2 <- as.integer(time > timeVector[j])
    bs=mean((0 - survs[j, ])^2 * help1 * (1/csurv) +
              (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j]))

    bs=as.numeric(bs)
    return(bs)
  }

  minus.roc <- function(par, FitCV, timeVector,
                        obj_surv, ot, time, pro.time, ROC.precision) {

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))

    .par<-exp(c(par))/(1+sum(exp(par)))

    .par<-c(.par,1-sum(.par))

    for (i in 1:length(FitCV)){ .pred[,,i]<-FitCV[[i]]*.par[i] }
    .pred<-rowSums(.pred, dims=2)

    survs <- t(.pred)[,ot]

    j=length(timeVector[which(timeVector<=pro.time)])

    .data <- data.frame(times = time, failures = obj_surv[ot,2], predictions=1-survs[j, ])

    return(-1*roc(times="times", failures="failures", variable="predictions", confounders=~1, data=.data,
             pro.time=pro.time, precision=ROC.precision)$auc)

  }


  ibll <- function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))


    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    survs <- t(.pred)[,ot]
    survs[which(survs==0)]<-10**-7
    survs[which(survs==1)]<-1-10**-7

    bll<-sapply(1:length(timeVector), FUN = function(j)
    {
      help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
      help2 <- as.integer(time > timeVector[j])
      return(-mean(log(1 - survs[j, ]) * help1 * (1/csurv) +
                     log( survs[j, ]) * help2 * (1/csurv_btime[j])))
    })

    idx <- 2:length(timeVector)
    RET <- diff(timeVector) %*% ((bll[idx - 1] + bll[idx])/2)
    RET <- RET/diff(range(timeVector))
    RET=as.matrix(RET)

    return(RET)
  }

  bll<-function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time, pro.time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))

    .par<-exp(c(par))/(1+sum(exp(par)))

    .par<-c(.par,1-sum(.par))

    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    survs <- t(.pred)[,ot]
    survs[which(survs==0)]<-10**-7
    survs[which(survs==1)]<-1-10**-7

    j=length(timeVector[which(timeVector<=pro.time)])


    help1 <- as.integer(time <= timeVector[j] &  obj_surv[ot,2] == 1)
    help2 <- as.integer(time > timeVector[j])
    bll=-mean(log(1- survs[j, ]) * help1 * (1/csurv) +
                log(survs[j, ]) * help2 * (1/csurv_btime[j]))

    bll=as.numeric(bll)
    return(bll)
  }

  ribs<-function(par, FitCV, timeVector, obj_surv, ot, csurv, csurv_btime, time, pro.time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))


    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    .pred=.pred[,timeVector<=pro.time]

    survs <- t(.pred)[,ot]

    timeVector=timeVector[timeVector<=pro.time]

    bsc<-sapply(1:length(timeVector), FUN = function(j)
    {
      help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
      help2 <- as.integer(time > timeVector[j])
      return(mean((0 - survs[j, ])^2 * help1 * (1/csurv) +
                    (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j])))
    })

    idx <- 2:length(timeVector)
    RET <- diff(timeVector) %*% ((bsc[idx - 1] + bsc[idx])/2)
    RET <- RET/diff(range(timeVector))
    RET=as.matrix(RET)

    return(RET)
  }


  ribll<-function(par, FitCV, timeVector,  obj_surv, ot, csurv, csurv_btime, time, pro.time){

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))


    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)

    .pred=.pred[,timeVector<=pro.time]

    survs <- t(.pred)[,ot]

    timeVector=timeVector[timeVector<=pro.time]

    bll<-sapply(1:length(timeVector), FUN = function(j)
    {
      help1 <- as.integer(time <= timeVector[j] & obj_surv[ot,2] == 1)
      help2 <- as.integer(time > timeVector[j])
      return(-mean(log(1 - survs[j, ]) * help1 * (1/csurv) +
                     log( survs[j, ]) * help2 * (1/csurv_btime[j])))
    })

    idx <- 2:length(timeVector)
    RET <- diff(timeVector) %*% ((bll[idx - 1] + bll[idx])/2)
    RET <- RET/diff(range(timeVector))
    RET=as.matrix(RET)

    return(RET)
  }

  minus.ci <- function(par, FitCV, data.times, data.failures, time.pred,
                        pro.time) {

    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))

    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred <- rowSums(.pred, dims=2)

    predicted <- .pred[,time.pred>=pro.time][,1]

    timeVector<-sort(unique(time.pred))

    obj_surv <- Surv(data.times, data.failures)
    time <- obj_surv[, 1]
    status <- obj_surv[, 2]

    permissible <- 0
    concord <- 0
    par_concord <- 0

    n <- length(time)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if ((time[i] < time[j] &
             status[i] == 0) | (time[j] < time[i] & status[j] == 0)) {
          next
        }

        if (time[i] == time[j] & status[i] == 0 & status[j] == 0) {
          next
        }

        permissible <- permissible + 1

        if (time[i] != time[j]) {
          if ((time[i] < time[j] &
               predicted[i] < predicted[j]) |
              (time[j] < time[i] & predicted[j] < predicted[i])) {
            concord <- concord + 1
          } else if (predicted[i] == predicted[j]) {
            par_concord <- par_concord + 0.5
          }
        }

        if (time[i] == time[j] & status[i] == 1 & status[j] == 1) {
          if (predicted[i] == predicted[j]) {
            concord <- concord + 1
          } else {
            par_concord <- par_concord + 0.5
          }
        }

        if (time[i] == time[j] &
            ((status[i] == 1 &
              status[j] == 0) | (status[i] == 0 & status[j] == 1))) {
          if ((status[i] == 1 &
               predicted[i] < predicted[j]) |
              (status[j] == 1 & predicted[j] < predicted[i])) {
            concord <- concord + 1
          } else {
            par_concord <- par_concord + 0.5
          }
        }
      }
    }

    return( -1 * (concord + par_concord) / permissible)
  }


  ###################################################
  ### Initialisation et recuperation param.tune ###
  ###################################################

  if(sum(!(methods %in% c("COXlasso", "COXridge", "RSF", "SNN", "COXen",
                          "AFTweibull","AFTweibull","AFTggamma","AFTgamma",
                          "PHgompertz","PHexponential",
                          "AFTllogis","COXaic","COXall", "PHspline")))>=1){
    stop("New method is not yet implemented") }

  M<-length((methods))
  N <- length(data[,times])

  max.progess <- M + cv * M + 4
  pb <- txtProgressBar(min = 0, max = max.progess, style = 3, width = 50, char = "=")
  ip <- 0
  setTxtProgressBar(pb, ip)

  names.meth=c(rep(NA,M))
  for(i in unique(methods)){
    if(length(which(methods==i))==1){
      names.meth[which(methods==i)]=i
    }
    else{
      names.meth[which(methods==i)]=paste0(i,1:length(which(methods==i)))
    }
  }

  time.pred <- sort(unique(data[,times]))

  if(is.null(param.tune)==FALSE){
    if(length(param.tune)!=M){
      stop("Param.tune need to have one element per method. Please modifiy param.tune or set it = NULL")
    }
    for (me in 1:M){
      if(is.null(param.tune[[me]])==T & !(methods[me] %in% c("AFTgamma",
                            "AFTggamma","AFTweibull","AFTllogis","PHexponential","PHgompertz"))){
        if(methods[me] %in%"PHspline"){
          param.tune[[me]]=list(k=1:4)
        }
        if(methods[me] %in%"COXen"){
          param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
        }
        if(methods[me] %in%"COXaic"){
          param.tune[[me]]=list(final.model=NA, model.min=NULL, model.max=NULL)
        }
        if(methods[me] %in%"SNN"){
          param.tune[[me]]=list(n.nodes=c(2, 3, 4, 6, 10, 20),
                                decay=c(0, 0.01, 0.1),
                                batch.size=256L,
                                epochs=1L)
        }
        if(methods[me] %in% "COXlasso"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"COXridge"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"RSF"){
          param.tune[[me]]=list(mtry=(length(group)+length(cov.quanti)+length(cov.quali))/2+2,
                                nodesize=c(2, 4, 6, 10, 20, 30, 50, 100),
                                ntree=500)
        }
      }
    }
  }


  if(is.null(param.tune)){
    param.tune=vector("list",M)
    for (me in 1:M){
      if(methods[me] %in%"COXen"){
        param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
      }
      if(methods[me] %in%"COXaic"){
        param.tune[[me]]=list(final.model=NA, model.min=NULL,model.max=NULL)
      }
      if(methods[me] %in%"SNN"){
        param.tune[[me]]=list(n.nodes=c(2, 3, 4, 6, 10, 20),
                              decay=c(0, 0.01, 0.1),
                              batch.size=256L,
                              epochs=1L)
      }
      if(methods[me] %in% "COXlasso"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in% "COXridge"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in% "PHspline"){
        param.tune[[me]]=list(k=1:4)
      }
      if(methods[me] %in% "RSF"){
        param.tune[[me]]=list(mtry=seq(1,(length(group)+length(cov.quanti)+length(cov.quali))/2+2),
                              nodesize=c(2, 4, 6, 10, 20, 30, 50, 100), ntree=500)
      }
    }
  }

  if(is.null(pro.time)==T & sum(metric %in% c("ci","bs","bll","ribs","ribll","roc"))){
    pro.time=median(data[,times])
  }

  .model<-vector("list",M)
  names(.model)<-names.meth

  .tune.optimal<-vector("list",M)
  names(.tune.optimal)<-names.meth

  .tune.results<-vector("list",M)
  names(.tune.results)<-names.meth


  for (me in 1:M){

    if(methods[me] == "AFTweibull" ){
      .tune.optimal[[me]]=NA

      .AFTweibull <- AFTweibull(times=times, failures=failures, group=group,
                                  cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.AFTweibull
      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.AFTweibull) }

    if(methods[me] == "AFTggamma"){
      .tune.optimal[[me]]=NA

      .AFTggamma <- AFTggamma(times=times, failures=failures, group=group,
                                cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.AFTggamma

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.AFTggamma) }

    if(methods[me] == "AFTgamma" ){
      .tune.optimal[[me]]=NA

      .AFTgamma <- AFTgamma(times=times, failures=failures, group=group,
                              cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.AFTgamma

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.AFTgamma)  }


    if(methods[me] == "AFTllogis" ){
      .tune.optimal[[me]]=NA

      .AFTllogis <- AFTllogis(times=times, failures=failures, group=group,
                              cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.AFTllogis
      rm(.AFTllogis)    }

    if(methods[me] == "PHgompertz" ){
      .tune.optimal[[me]]=NA

      .PHgompertz <- PHgompertz(times=times, failures=failures, group=group,
                                  cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.PHgompertz

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.PHgompertz) }

    if(methods[me] == "PHexponential"){
      .tune.optimal[[me]]=NA

      .PHexponential <- PHexponential(times=times, failures=failures, group=group,
                                        cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.PHexponential

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.PHexponential) }

    if(methods[me] == "COXall"){
      .tune.optimal[[me]]=NA

      .coxall <- COXall(times=times, failures=failures, group=group,
                                        cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)

      .model[[me]]<-.coxall

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.coxall) }

    if(methods[me] == "PHspline"){

      if(is.null(param.tune[[me]]$k)==T | length(param.tune[[me]]$k)>1){

        .tune<- tunePHspline(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                               cov.quali=cov.quali, data=data, cv=cv,
                               k=param.tune[[me]]$k)
        .tune.optimal[[me]]=.tune$optimal
        .tune.results[[me]]=.tune$results
        rm(.tune)  }

      else{ .tune.optimal[[me]]=list(lambda=param.tune[[me]]$k) }


      .PHspline <- PHspline(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                              cov.quali=cov.quali, data=data, k=.tune.optimal[[me]]$k)

      .model[[me]]<-.PHspline

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.PHspline)    }

    if(methods[me] == "COXlasso"){

      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){

        .tune<- tuneCOXlasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                               cov.quali=cov.quali, data=data, cv=cv,
                               parallel=FALSE, lambda=param.tune[[me]]$lambda)

        .tune.optimal[[me]]=.tune$optimal
        .tune.results[[me]]=.tune$results
        rm(.tune)  }

      else{ .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda) }


      .COXlasso <- COXlasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                              cov.quali=cov.quali, data=data, lambda=.tune.optimal[[me]]$lambda)

      .model[[me]]<-.COXlasso

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.COXlasso)  }

    if(methods[me] == "COXridge"){

      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){
        .tune<- tuneCOXridge(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                               cov.quali=cov.quali,  data=data, cv=cv,
                               parallel = FALSE, lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda)
      }

      .COXridge <- COXridge(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                              cov.quali=cov.quali, data=data,
                              lambda=.tune.optimal[[me]]$lambda)
      .model[[me]]<-.COXridge

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.COXridge)    }

    if(methods[me] == "COXen"){

      if(length(param.tune[[me]]$alpha)==1 & length(param.tune[[me]]$lambda)==1){
        .tune.optimal[[me]]=list(alpha=param.tune[[me]]$alpha, lambda=param.tune[[me]]$lambda)
      }
      else{
        .tune <- tuneCOXen(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                             cov.quali=cov.quali, data=data, cv=cv,
                             parallel=FALSE,
                             alpha=param.tune[[me]]$alpha,
                             lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results

        rm(.tune) }

      .COXen <- COXen(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                        cov.quali=cov.quali, data=data, alpha=.tune.optimal[[me]]$alpha,
                        lambda=.tune.optimal[[me]]$lambda)

      .model[[me]]<-.COXen

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.COXen) }

    if(methods[me] == "COXaic"){

      if(is.na(param.tune[[me]]$final.model)==FALSE){
        .tune.optimal[[me]]=list(final.model=param.tune[[me]]$final.model)
      }
      else{

        .tune <- tuneCOXaic(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                             cov.quali=cov.quali, data=data,
                             model.min=param.tune[[me]]$model.min,
                             model.max=param.tune[[me]]$model.max)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results

        rm(.tune) }

      .COXaic<- COXaic(times=times, failures=failures, group=group, data=data,
                       cov.quanti=cov.quanti, cov.quali=cov.quali,
                       final.model = .tune.optimal[[me]]$final.model)

      .model[[me]]<-.COXaic

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.COXaic)  }

    if (methods[me] == "RSF"){

      if(length(param.tune[[me]]$nodesize)!=1 | length(param.tune[[me]]$mtry)!=1 |
         length(param.tune[[me]]$ntree)!=1){
        .tune<-tuneRSF(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                            cov.quali=cov.quali, data=data,
                            nodesize=param.tune[[me]]$nodesize,
                            mtry=param.tune[[me]]$mtry,
                            ntree=param.tune[[me]]$ntree)

        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]<-list(nodesize=param.tune[[me]]$nodesize,
                                  mtry=param.tune[[me]]$mtry,
                                  ntree=param.tune[[me]]$ntree)
      }
      .RSF <-RSF(times=times, failures=failures,
                         group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=data,
                         nodesize=.tune.optimal[[me]]$nodesize,
                         mtry=.tune.optimal[[me]]$mtry, ntree=.tune.optimal[[me]]$ntree)

      .model[[me]]<-.RSF

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.RSF) }


    if (methods[me] == "SNN"){
      torch<-reticulate::import("torch")
      torch$set_num_threads(1L)

      if(length(param.tune[[me]]$n.nodes)!=1 | length(param.tune[[me]]$decay)!=1 |
         length(param.tune[[me]]$batch.size)!=1 |length(param.tune[[me]]$epochs)!=1 ){

        .tune <- tuneSNN(times=times, failures=failures, group=group,
                               cov.quanti=cov.quanti, cov.quali=cov.quali,
                               data=data, cv=cv,
                               n.nodes=param.tune[[me]]$n.nodes,
                               decay=param.tune[[me]]$decay,
                               batch.size=param.tune[[me]]$batch.size,
                               epochs=param.tune[[me]]$epochs)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]<-list(n.nodes=param.tune[[me]]$n.nodes,
                                  decay=param.tune[[me]]$decay,
                                  batch.size=param.tune[[me]]$batch.size,
                                  epochs=param.tune[[me]]$epochs)
      }

      .SNN <-SNN(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                         cov.quali=cov.quali, data=data,
                         n.nodes=as.numeric(.tune.optimal[[me]]$n.nodes),
                         decay=as.numeric(.tune.optimal[[me]]$decay),
                         batch.size=as.integer(.tune.optimal[[me]]$batch.size),
                         epochs=as.integer(.tune.optimal[[me]]$epochs))


      .model[[me]]<-.SNN

      ip <- ip+1
      setTxtProgressBar(pb, ip)
      rm(.SNN)    }
  }

  ########################
  ### Cross-Validation  ##
  ########################

  sample_id=sample(nrow(data))
  folds <- cut(seq(1,nrow(data)), breaks=cv, labels=FALSE)
  folds_id=folds[sample_id]
  data$folds=folds_id

  CV<-vector("list",cv*M)
  j<-1
  for(m in 1:M){
    for (k in 1:cv){
      CV[[j]]<-list(train= data[data$folds!=k, ],valid=data[data$folds==k, ], num_method=m)
      j<-j+1
    }
  }

  CV_all_method<-function(CV, method, Tune,
                          times, failures, group, cov.quanti, cov.quali,time.pred){
    num_method<-CV$num_method
    meth<-method[num_method]
    if(meth == "AFTweibull"){
      fit<-AFTweibull(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                       cov.quali=cov.quali, data=CV$train)
      pred=predict(fit,  newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "AFTggamma"){
      fit<-AFTggamma(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                      cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "AFTgamma"){
      fit<-AFTgamma(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "AFTllogis"){
      fit<-AFTllogis(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "PHgompertz"){
      fit<-PHgompertz(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                       cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "PHexponential"){
      fit<-PHexponential(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                          cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "PHspline"){
      fit<-PHspline(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train,
                     k=Tune[[num_method]]$k)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "COXlasso"){
      fit<-COXlasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train,
                     lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "COXen"){
      fit<-COXen(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                 cov.quali=cov.quali, data=CV$train,
                  alpha=Tune[[num_method]]$alpha, lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth =="COXridge"){
      fit<-COXridge(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                    cov.quali=cov.quali, data=CV$train, lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="COXaic"){
      fit<-COXaic(times=times, failures=failures, group=group, data=data,
                   final.model = Tune[[num_method]]$final.model, cov.quanti=cov.quanti,
                   cov.quali=cov.quali)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="COXall"){
      fit<-COXall(times=times, failures=failures, group=group,
                   cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="RSF"){
      fit<-RSF(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                   cov.quali=cov.quali,  data=CV$train,
                   nodesize=Tune[[num_method]]$nodesize, mtry=Tune[[num_method]]$mtry,
                   ntree=Tune[[num_method]]$ntree)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="SNN"){
      fit<-SNN(times=times, failures=failures, group=group,  cov.quanti=cov.quanti, cov.quali=cov.quali, data=CV$train,
                     n.nodes=as.numeric(Tune[[num_method]]$n.nodes),
                     decay=as.numeric(Tune[[num_method]]$decay),
                     batch.size=as.integer(Tune[[num_method]]$batch.size),
                     epochs=as.integer(Tune[[num_method]]$epochs))
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    return(pred)
  }

  preFitCV<-lapply(CV, CV_all_method,method=methods,
                   Tune=.tune.optimal, times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                   cov.quali=cov.quali, time.pred=time.pred)

  FitCV<-vector("list", M)
  for(m in 1:M){
    FitCV[[m]]<-matrix(nrow=dim(.model[[1]]$predictions)[1], ncol= dim(.model[[1]]$predictions)[2])
    for (j in 1:cv){
      FitCV[[m]][data$folds==j,]<-preFitCV[[(m-1)*cv+j]]

      ip <- ip+1
      setTxtProgressBar(pb, ip)
    }
  }
  names(FitCV)<-names.meth

  rm(preFitCV, CV)

  ################
  # OPTIMISATION #
  ################

  data.times <- data[,times]
  data.failures <- data[,failures]
  timeVector <- survfit(Surv(data[,times],data[,failures])~ 1 )$time

  obj_surv <- Surv(data.times, data.failures)

  time <- obj_surv[, 1]
  ot <- order(time)
  cens <- obj_surv[ot, 2]
  time <- time[ot]

  .temp <- survfit(Surv(time, cens==0) ~ 1)
  csurv <- summary(.temp, times=time, extend=TRUE)$surv

  csurv[csurv == 0] <- Inf

  csurv_btime <- summary(.temp, times=timeVector, extend=TRUE)$surv

  csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
  csurv_btime[csurv_btime == 0] <- Inf

  .optim.method="Nelder-Mead"
  if(M==2){  .optim.method="BFGS"  }

  if(is.null(param.weights.fix)==TRUE){

      switch(metric,
             ibs={
               estim<-optim(par=param.weights.init, fn=ibs, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=ibs, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime, time=time,hessian=F,
                              method=.optim.method)
               }
             },
             bs={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<-optim(par=param.weights.init, fn=brs, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=brs, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                              method=.optim.method)
               }
             },
             auc={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<-optim(par=param.weights.init, fn=minus.roc, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, ROC.precision = ROC.precision,
                            time=time,pro.time=pro.time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=minus.roc, FitCV = FitCV,  ROC.precision = ROC.precision,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot,
                              time=time,pro.time=pro.time,hessian=F,
                              method=.optim.method)
               }
             },
             ibll={
               estim<-optim(par=param.weights.init, fn=ibll, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=ibll, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime, time=time,hessian=F,
                              method=.optim.method)
               }
             },
             bll={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<-optim(par=param.weights.init, fn=bll, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=bll, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                              method=.optim.method)
               }
             },
             ribs={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<-optim(par=param.weights.init, fn=ribs, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=ribs, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                              method=.optim.method)
               }
             },
             ribll={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<-optim(par=param.weights.init, fn=ribll, FitCV = FitCV,
                            timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                            method=.optim.method)

               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=ribll, FitCV = FitCV,
                              timeVector=timeVector,  obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime,time=time,pro.time=pro.time,hessian=F,
                              method=.optim.method)
               }
             },
             ci={
               if(is.null(pro.time)){
                 pro.time=median(data[,times])
               }
               estim<- optim(par=param.weights.init, fn=minus.ci, FitCV = FitCV,
                             data.times =  data.times, data.failures = data.failures,
                             time.pred=time.pred, pro.time=pro.time, hessian=F, method=.optim.method)
               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par,fn=minus.ci, FitCV = FitCV, data.times =  data.times,
                              data.failures = data.failures, time.pred=time.pred,
                              pro.time=pro.time, hessian=F, method=.optim.method)
               }
             }
      )
  }

  ip <- ip+3
  setTxtProgressBar(pb, ip)

  ############################
  # Compute Survival from SL #
  ############################

  if(is.null(param.weights.fix)==FALSE){
    estim=list()
    estim$par=param.weights.fix
  }

  FitALL<-vector("list",M)
  names(FitALL)<-names.meth

  for(me in 1:M){
    FitALL[[me]]<-.model[[me]]$predictions
  }

w.sl <- c(exp(c(estim$par,0)) / ( 1+sum(exp(estim$par))) )
.SL<-array(dim = c(dim(FitALL[[1]]),length(FitALL)))
for (i in 1:length(FitCV)){ .SL[,,i]<-.model[[i]]$predictions*w.sl[i]  }
surv.SL <-rowSums(.SL, dims=2)


  ######################
  # PREPARATION RETURN #
  ######################

  if(keep.predictions==TRUE) {
    FitALL<-vector("list",M)
    names(FitALL)<-names.meth

    for(me in 1:M){
      FitALL[[me]]<-.model[[me]]$predictions
    }

    FitALL$sl<-surv.SL
    temp.predictions <- FitALL}

  else {
    temp.predictions <- surv.SL
    temp.predictions<-as.data.frame(temp.predictions)
  }

ip <- ip+1
setTxtProgressBar(pb, ip)

  res<-list(times=time.pred,
            predictions=temp.predictions,
            data=data.frame(times=data[,times], failures=data[,failures],
                            data[, !(dimnames(data)[[2]] %in% c(times, failures))]),
            outcomes=list(times=times, failures=failures),
            predictors=list(group=group, cov.quanti=cov.quanti, cov.quali=cov.quali),
            ROC.precision=ROC.precision,
            cv=cv,
            pro.time=pro.time,
            methods=methods,
            models=.model,
            weights=list(coefficients=estim$par, values=w.sl),
            metric=list(metric=metric, value=estim$value),
            param.tune=list(optimal=.tune.optimal, results=.tune.results))

  class(res) <- "sltime"

  close(pb)

  return(res)
}

