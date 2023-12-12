
sl.time <- function( methods=c("cox.lasso", "aft.ggamma"),
                     metric="bs",  data, times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL,
                     cv=10, param.tune=NULL, pro.time=NULL,  optim.local.min=FALSE, 
                     ROC.precision=seq(.01,.99,.01), param.weights.fix=NULL,
                     param.weights.init=NULL, keep.predictions=TRUE,
                     verbose=TRUE) {
  
  ###################################################
  ### Initialisation et autres tests de coherence ###
  ###################################################
  
  if(length(methods)<=1)
  { stop("Number of methods need to be greater or equal to 2 to estimate a SuperLearner")   }
  
  if(length(metric)>1){
    warning(paste0("SuperLearner is currently developped for one metric. Results for metric ",metric[1]))
    metric=metric[1]
  }
  
  if(min(metric%in%c("bs","loglik","ibs","ibll","bll",
                     "ribs","ribll","auc",
                     "loglik.td"))==0){
    stop("The argument \"metric\" must be Brier score (bs),
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
  
  # keep only needed column of data
  if(!is.null(group)==T){
    data<-data[,c(times,failures,group,cov.quanti,cov.quali)]
  }
  if(is.null(group)==T){
    data<-data[,c(times,failures,cov.quanti,cov.quali)]
  }
  
  # check valeur manquante dans data
  if (any(is.na(data))){
    data<-na.omit(data)   # CAMILLE : J'ai remplace par na.omit car je ne trouve pas cette fonction na.rm()
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
  
  # check sur les temps d'evenement et indicatrice d'event
  #a supprimer car check valeur manquantes avant
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
  
  
  # check CV >2
  if (cv < 3 | !is.numeric(cv)) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  
  
  # CAMILLE TODO e revoir entierement
  # faire les autres checks des param.tune
  # check methods pour detecter des incoherences
  # ici si plusieurs methodes AFT ou PH , impossible, il faut en sortir
  .meth_rm=c()
  if(sum(methods %in% "aft.gamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="aft.gamma")[-1])
    warning("SuperLearner can use only one aft.gamma method. We remove the others.")
  }
  if(sum(methods %in% "aft.llogis")>=2){
    .meth_rm=c(.meth_rm,which(methods=="aft.llogis")[-1])
    warning("SuperLearner can use only one aft.llogis method. We remove the others.")
  }
  if(sum(methods %in% "aft.ggamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="aft.ggamma")[-1])
    warning("SuperLearner can use only one aft.ggamma method. We remove the others.")
  }
  if(sum(methods %in% "aft.weibull")>=2){
    .meth_rm=c(.meth_rm,which(methods=="aft.weibull")[-1])
    warning("SuperLearner can use only one aft.weibull method. We remove the others.")
  }
  if(sum(methods %in% "ph.exponential")>=2){
    .meth_rm=c(.meth_rm,which(methods=="ph.exponential")[-1])
    warning("SuperLearner can use only one ph.exponential method. We remove the others.")
  }
  if(sum(methods %in% "ph.gompertz")>=2){
    .meth_rm=c(.meth_rm,which(methods=="ph.gompertz")[-1])
    warning("SuperLearner can use only one ph.gompertz method. We remove the others.")
  }
  
  if(sum(methods %in% "cox.lasso")==1){ #si une seule methode lasso, check si param.tune OK
    if(!(is.null(param.tune[[which(methods=="cox.lasso")]]))){
      if(!is.list(param.tune[[which(methods=="cox.lasso")]])){
        stop("Argument param.tune for cox.lasso need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="cox.lasso")]])%in%"lambda"))==0){
        stop("Tune parameters for cox.lasso need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="cox.lasso")]]$lambda)|
           is.null(param.tune[[which(methods=="cox.lasso")]]$lambda))){
        stop("Lambda tune parameters for cox.lasso need to be a scalar or a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "cox.lasso")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="cox.lasso")])!=length(unique(param.tune[which(methods=="cox.lasso")]))){
      stop("Tune parameters for cox.lasso methods need to be unique")
    }
    for (i in 1:sum(methods %in% "cox.lasso")){
      if(!(is.null(param.tune[[which(methods=="cox.lasso")[i]]]))){
        if(!is.list(param.tune[[which(methods=="cox.lasso")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th cox.lasso need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="cox.lasso")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.lasso need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="cox.lasso")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="cox.lasso")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th cox.lasso need to be a scalar or a vector or NULL"))
        }
      }
    }
  }
  
  if(sum(methods %in% "cox.ridge")==1){ #si une seule methode lasso, check si param.tune OK
    if(!(is.null(param.tune[[which(methods=="cox.ridge")]]))){
      if(!is.list(param.tune[[which(methods=="cox.ridge")]])){
        stop("Argument param.tune for cox.ridge need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="cox.ridge")]])%in%"lambda"))==0){
        stop("Tune parameters for cox.ridge need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="cox.ridge")]]$lambda)|
           is.null(param.tune[[which(methods=="cox.ridge")]]$lambda))){
        stop("Lambda tune parameters for cox.ridge need to be a scalar or a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "cox.ridge")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="cox.ridge")])!=length(unique(param.tune[which(methods=="cox.ridge")]))){
      stop("Tune parameters for cox.ridge methods need to be unique")
    }
    for (i in 1:sum(methods %in% "cox.ridge")){
      if(!(is.null(param.tune[[which(methods=="cox.ridge")[i]]]))){
        if(!is.list(param.tune[[which(methods=="cox.ridge")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th cox.ridge need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="cox.ridge")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.ridge need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="cox.ridge")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="cox.ridge")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th cox.ridge need to be a scalar or a vector or NULL"))
        }
      }
    }
  }
  
  if(sum(methods %in% "cox.en")==1){ #si une seule methode lasso, check si param.tune OK
    if(!(is.null(param.tune[[which(methods=="cox.en")]]))){
      if(!is.list(param.tune[[which(methods=="cox.en")]])){
        stop("Argument param.tune for cox.en need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="cox.en")]])%in%"lambda"))==0){
        stop("Tune parameters for cox.en need to have lambda")
      }
      if(sum((names(param.tune[[which(methods=="cox.en")]])%in%"alpha"))==0){
        stop("Tune parameters for cox.en need to have alpha")
      }
      if(!(is.numeric(param.tune[[which(methods=="cox.en")]]$lambda)|
           is.null(param.tune[[which(methods=="cox.en")]]$lambda))){
        stop("Lambda tune parameters for cox.en need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="cox.en")]]$alpha)|
           is.null(param.tune[[which(methods=="cox.en")]]$alpha))){
        stop("alpha tune parameters for cox.en need to be a scalar or a vector or NULL")
      }
      if(min(param.tune[[which(methods=="cox.en")]]$alpha)<0 | max(param.tune[[which(methods=="cox.en")]]$alpha)>1){
        stop("tune parameters for cox.en alpha need to be in ]0;1[")
      }
    }
  }
  if(sum(methods %in% "cox.en")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="cox.en")])!=length(unique(param.tune[which(methods=="cox.en")]))){
      stop("Tune parameters for cox.en methods need to be unique")
    }
    for (i in 1:sum(methods %in% "cox.en")){
      if(!(is.null(param.tune[[which(methods=="cox.en")[i]]]))){
        if(!is.list(param.tune[[which(methods=="cox.en")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th cox.en need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="cox.en")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.en need to have lambda"))
        }
        if(sum((names(param.tune[[which(methods=="cox.en")[i]]])%in%"alpha"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.en need to have alpha"))
        }
        if(!(is.numeric(param.tune[[which(methods=="cox.en")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="cox.en")[i]]]$lambda))){
          stop(paste("Lambda tune parameters for the ",i,"th cox.en need to be a scalar or a vector or NULL"))
        }
        if(!(is.numeric(param.tune[[which(methods=="cox.en")[i]]]$alpha)|
             is.null(param.tune[[which(methods=="cox.en")[i]]]$alpha))){
          stop(paste("Alpha tune parameters for the ",i,"th cox.en need to be a scalar or a vector or NULL"))
        }
        if(min(param.tune[[which(methods=="cox.en")[i]]]$alpha)<0 | max(param.tune[[which(methods=="cox.en")[i]]]$alpha)>1){
          stop("tune parameters for cox.en alpha need to be in ]0;1[")
        }
      }
    }
  }
  
  if(sum(methods %in% "cox.aic")==1){ 
    if(!(is.null(param.tune[[which(methods=="cox.aic")]]))){
      if(!is.list(param.tune[[which(methods=="cox.aic")]])){
        stop("Argument param.tune for cox.aic need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="cox.aic")]])%in%"final.model"))==0){
        stop("Tune parameters for cox.aic need to have final.model")
      }
      if(sum((names(param.tune[[which(methods=="cox.aic")]])%in%"model.min"))==0){
        stop("Tune parameters for cox.aic need to have model.min")
      }
      if(sum((names(param.tune[[which(methods=="cox.aic")]])%in%"model.max"))==0){
        stop("Tune parameters for cox.aic need to have model.max")
      }
      # CS YF voir pour lec checks sur les autres parametres
      # if(!(is.numeric(param.tune[[which(methods=="cox.univ")]]$minscreen)|
      #      is.null(param.tune[[which(methods=="cox.univ")]]$minscreen))){
      #   stop("minscreen tune parameters for cox.univ need to be a scalar or a vector or NULL")
      # }
      # if(!(is.numeric(param.tune[[which(methods=="cox.univ")]]$min.p)|
      #      is.null(param.tune[[which(methods=="cox.univ")]]$min.p))){
      #   stop("min.p tune parameters for cox.univ need to be a scalar or a vector or NULL")
      # }
      # if(min(param.tune[[which(methods=="cox.univ")]]$min.p)<0 | max(param.tune[[which(methods=="cox.univ")]]$min.p)>1){
      #   stop("tune parameters for cox.univ min.p need to be in ]0;1[")
      # }
      # if(all(param.tune[[which(methods=="cox.univ")]]$cov.select%in%c(cov.quanti,cov.quali))==FALSE & 
      #    is.na(param.tune[[which(methods=="cox.univ")]]$cov.select)==FALSE){
      #   stop("At least one cov.select for tune parameters for cox.univ is not present in cov.quanti or cov.quali")
      # }
    }
  }
  if(sum(methods %in% "cox.aic")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="cox.aic")])!=length(unique(param.tune[which(methods=="cox.aic")]))){
      stop("Tune parameters for cox.aic methods need to be unique")
    }
    for (i in 1:sum(methods %in% "cox.aic")){
      if(!(is.null(param.tune[[which(methods=="cox.aic")[i]]]))){
        if(!is.list(param.tune[[which(methods=="cox.aic")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th cox.aic need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="cox.aic")[i]]])%in%"finl.model.cov"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.aic need to have finl.model.cov"))
        }
        if(sum((names(param.tune[[which(methods=="cox.aic")[i]]])%in%"model.min"))==0){
          stop(paste("Tune parameters for the ",i,"th cox.aic need to have model.min"))
        }
        if(sum((names(param.tune[[which(methods=="cox.aic")[i]]])%in%"model.max"))==0){
          stop("Tune parameters for cox.aic need to have model.max")
        }
      #   if(!(is.numeric(param.tune[[which(methods=="cox.univ")[i]]]$minscreen)|
      #        is.null(param.tune[[which(methods=="cox.univ")[i]]]$minscreen))){
      #     stop(paste("minscreen tune parameters for the ",i,"th cox.univ need to be a scalar or a vector or NULL"))
      #   }
      #   if(!(is.numeric(param.tune[[which(methods=="cox.univ")[i]]]$min.p)|
      #        is.null(param.tune[[which(methods=="cox.univ")[i]]]$min.p))){
      #     stop(paste("min.p tune parameters for the ",i,"th cox.univ need to be a scalar or a vector or NULL"))
      #   }
      #   if(min(param.tune[[which(methods=="cox.univ")[i]]]$min.p)<0 | max(param.tune[[which(methods=="cox.univ")[i]]]$min.p)>1){
      #     stop("tune parameters for cox.univ min.p need to be in ]0;1[")
      #   }
      #   if(all(param.tune[[which(methods=="cox.univ")[i]]]$cov.select%in%c(cov.quanti,cov.quali))==FALSE & 
      #      is.na(param.tune[[which(methods=="cox.univ")[i]]]$cov.select)==FALSE){
      #     stop("At least one cov.select for tune parameters for cox.univ is not present in cov.quanti or cov.quali")
      #   }
      }
    }
  }
  
  
  if(sum(methods %in% "rf.time")==1){ #si une seule methode lasso, check si param.tune OK
    if(!(is.null(param.tune[[which(methods=="rf.time")]]))){
      if(!is.list(param.tune[[which(methods=="rf.time")]])){
        stop("Argument param.tune for rf.time need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="rf.time")]])%in%"nodesize"))==0){
        stop("Tune parameters for rf.time need to have nodesize")
      }
      if(sum((names(param.tune[[which(methods=="rf.time")]])%in%"mtry"))==0){
        stop("Tune parameters for rf.time need to have mtry")
      }
      if(sum((names(param.tune[[which(methods=="rf.time")]])%in%"ntree"))==0){
        stop("Tune parameters for rf.time need to have ntree")
      }
      if(!(is.numeric(param.tune[[which(methods=="rf.time")]]$nodesize)|
           is.null(param.tune[[which(methods=="rf.time")]]$nodesize))){
        stop("nodesize tune parameters for rf.time need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="rf.time")]]$mtry)|
           is.null(param.tune[[which(methods=="rf.time")]]$mtry))){
        stop("mtry tune parameters for rf.time need to be a scalar or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="rf.time")]]$ntree)|
           is.null(param.tune[[which(methods=="rf.time")]]$ntree))){
        stop("ntree tune parameters for rf.time need to be a scalar or NULL")
      }
    }
  }
  if(sum(methods %in% "rf.time")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="rf.time")])!=length(unique(param.tune[which(methods=="rf.time")]))){
      stop("Tune parameters for rf.time methods need to be unique")
    }
    for (i in 1:sum(methods %in% "rf.time")){
      if(!is.list(param.tune[[which(methods=="rf.time")[i]]])){
        stop(paste("Argument param.tune for the ",i,"th rf.time need to be a list"))
      }
      if(!(is.null(param.tune[[which(methods=="rf.time")[i]]]))){
        if(sum((names(param.tune[[which(methods=="rf.time")[i]]])%in%"nodesize"))==0){
          stop("Tune parameters for rf.time need to have nodesize")
        }
        if(sum((names(param.tune[[which(methods=="rf.time")[i]]])%in%"mtry"))==0){
          stop("Tune parameters for rf.time need to have mtry")
        }
        if(sum((names(param.tune[[which(methods=="rf.time")[i]]])%in%"ntree"))==0){
          stop("Tune parameters for rf.time need to have ntree")
        }
        if(!(is.numeric(param.tune[[which(methods=="rf.time")[i]]]$nodesize)|
             is.null(param.tune[[which(methods=="rf.time")[i]]]$nodesize))){
          stop("nodesize tune parameters for rf.time need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="rf.time")[i]]]$mtry)|
             is.null(param.tune[[which(methods=="rf.time")[i]]]$mtry))){
          stop("mtry tune parameters for rf.time need to be a scalar or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="rf.time")[i]]]$ntree)|
             is.null(param.tune[[which(methods=="rf.time")[i]]]$ntree))){
          stop("ntree tune parameters for rf.time need to be a scalar or NULL")
        }
      }
    }
  }
  
  if(sum(methods %in% "nn.time")==1){ #si une seule methode lasso, check si param.tune OK
    if(!(is.null(param.tune[[which(methods=="nn.time")]]))){
      if(!is.list(param.tune[[which(methods=="nn.time")]])){
        stop("Argument param.tune for nn.time need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="nn.time")]])%in%"n.nodes"))==0){
        stop("Tune parameters for nn.time need to have n.nodes")
      }
      if(sum((names(param.tune[[which(methods=="nn.time")]])%in%"decay"))==0){
        stop("Tune parameters for nn.time need to have decay")
      }
      if(sum((names(param.tune[[which(methods=="nn.time")]])%in%"batch.size"))==0){
        stop("Tune parameters for nn.time need to have batch.size")
      }
      if(sum((names(param.tune[[which(methods=="nn.time")]])%in%"epochs"))==0){
        stop("Tune parameters for nn.time need to have epochs")
      }
      if(!(is.numeric(param.tune[[which(methods=="nn.time")]]$n.nodes)|
           is.null(param.tune[[which(methods=="nn.time")]]$n.nodes))){
        stop("n.nodes tune parameters for nn.time need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="nn.time")]]$decay)|
           is.null(param.tune[[which(methods=="nn.time")]]$decay))){
        stop("decay tune parameters for nn.time need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="nn.time")]]$batch.size)|
           is.null(param.tune[[which(methods=="nn.time")]]$batch.size))){
        stop("batch.size tune parameters for nn.time need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="nn.time")]]$epochs)|
           is.null(param.tune[[which(methods=="nn.time")]]$epochs))){
        stop("epochs tune parameters for nn.time need to be a scalar, a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "nn.time")>=2){ #si plusieurs methodes, check si unique puis check comme une methode pour chacune
    if(length(param.tune[which(methods=="nn.time")])!=length(unique(param.tune[which(methods=="nn.time")]))){
      stop("Tune parameters for nn.time methods need to be unique")
    }
    for (i in 1:sum(methods %in% "nn.time")){
      if(!(is.null(param.tune[[which(methods=="nn.time")[i]]]))){
        if(!is.list(param.tune[[which(methods=="nn.time")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th nn.time need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="nn.time")[i]]])%in%"n.nodes"))==0){
          stop("Tune parameters for nn.time need to have n.nodes")
        }
        if(sum((names(param.tune[[which(methods=="nn.time")[i]]])%in%"decay"))==0){
          stop("Tune parameters for nn.time need to have decay")
        }
        if(sum((names(param.tune[[which(methods=="nn.time")[i]]])%in%"batch.size"))==0){
          stop("Tune parameters for nn.time need to have batch.size")
        }
        if(sum((names(param.tune[[which(methods=="nn.time")[i]]])%in%"epochs"))==0){
          stop("Tune parameters for nn.time need to have epochs")
        }
        if(!(is.numeric(param.tune[[which(methods=="nn.time")[i]]]$n.nodes)|
             is.null(param.tune[[which(methods=="nn.time")[i]]]$n.nodes))){
          stop("nodesize tune parameters for nn.time need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="nn.time")[i]]]$decay)|
             is.null(param.tune[[which(methods=="nn.time")[i]]]$decay))){
          stop("decay tune parameters for nn.time need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="nn.time")[i]]]$batch.size)|
             is.null(param.tune[[which(methods=="nn.time")[i]]]$batch.size))){
          stop("batch.size tune parameters for nn.time need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="nn.time")[i]]]$epochs)|
             is.null(param.tune[[which(methods=="nn.time")[i]]]$epochs))){
          stop("epochs tune parameters for nn.time need to be a scalar, a vector or NULL")
        }
      }
    }
  }
  
  if(length(.meth_rm)>=1){ #CS a garder ?
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
      if(metric=="loglik.td"){
        if(length(param.weights.fix)!=(length(methods))*2-2){
          stop("wrong lenth for param.weights.fix")
        }
      }
      if(metric!="loglik.td"){
        if(length(param.weights.fix)!=(length(methods)-1)){
          stop("wrong lenth for param.weights.fix")
        }
      }
    }
    if(is.null(param.weights.init)==FALSE){
      if(is.numeric(param.weights.init)==FALSE){
        stop("param.weights.init need to be numeric")
      }
      if(metric=="loglik.td"){
        if(length(param.weights.init)!=((length(methods)-1)*2)){
          stop("wrong lenth for param.weights.init")
        }
      }
      if(metric!="loglik.td"){
        if(length(param.weights.init)!=(length(methods)-1)){
          stop("wrong lenth for param.weights.init")
        }
      }
    }
  }
  if(is.null(param.weights.fix)==TRUE & is.null(param.weights.init)==TRUE){
    if(metric!="loglik.td"){
      param.weights.init<-rep(0,length(methods)-1)
    }
    if(metric=="loglik.td"){
      param.weights.init<-rep(0,((length(methods)-1)*2))
    }
  }
  
  
  if(verbose==T){
    print("check data ok")
  }
  
  

  ################################################
  ### Fonctions utilisees dans optim des poids ###
  ################################################
  minus.partial.loglik <- function(par, FitCV, data.times, data.failures, differentiation) {
    
    time.pred <- sort(unique(data.times))
    
    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))
    
    
    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)
    
    for (i in 1:dim(.pred)[1]){
      if(min(.pred[i,])==0){
        .pred[i,.pred[i,]==0]<-min(.pred[i,which(.pred[i,]!=0)])
      }
    }
    
    .surv <- .pred[, time.pred %in% unique(data.times[data.failures==1]) ]
    .time <- time.pred[time.pred %in% unique(data.times[data.failures==1]) ]
    .haz <- t(apply(.surv, FUN = function(s) { differentiation(x=.time, fx=-1*log(s)) }, MARGIN=1))
    .indic <- t(sapply(data.times, FUN = function(x) {1*(.time==x)} ))
    .haz.indic <- .haz * .indic
    
    .surv <- .pred
    .indic <- t(sapply(data.times, FUN = function(x) {1*(time.pred==x)} ))
    .surv.indic <- .surv * .indic
    
    .haz.indiv <- apply(.haz.indic, FUN = "sum", MARGIN=1)[data.failures==1]
    .surv.indiv <- apply(.surv.indic, FUN = "sum", MARGIN=1)
    
    return(-1*(sum(log( pmax(.haz.indiv, min(.haz.indiv[.haz.indiv!=0])) )) + sum(log(.surv.indiv))))
  }
  
  
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
  
  minus.roc <- function(par, FitCV, data.times, data.failures, time.pred,
                        pro.time, precision) {
    .pred<-array(dim = c(dim(FitCV[[1]]),length(FitCV)))
    .par<-exp(c(par))/(1+sum(exp(par)))
    
    
    .par<-c(.par,1-sum(.par))
    for (i in 1:length(FitCV)){
      .pred[,,i]<-FitCV[[i]]*.par[i]
    }
    .pred<-rowSums(.pred, dims=2)
    
    # .pred=.pred[,timeVector<=pro.time]
    
    
    .pred=1-.pred[,time.pred>=pro.time][,1]
    
    .data <- data.frame(times = data.times, failures = data.failures, predictions=.pred)
    
    # .roc <- roc.time(times="times", failures="failures", variable="predictions",
    #                  confounders=~1, data=.data, pro.time=pro.time,      
    #                  precision=precision)
    # 
    # return(-.roc$auc)
    .roc<-timeROC::timeROC(T = data.times, delta = data.failures,
                            marker = .pred,
                            cause = 1, iid = TRUE,
                            times = pro.time)
    
    return(-.roc$AUC[2]*1)
  }
  
  # CS -> YF : 3 fonctions a minimiser : deux loglik totales 1 avec poids SL dep du temps, 
  # l autre avec poids SL non dep du temps
  # troisieme fonction qui fait du IBS avec poids SL dep du temps : attention c est tres long pour ibs.timedep
  minus.loglik.timedep <- function(par, haz,M, data.times, data.failures,time.event) {
    
    
    par.alpha=par[seq(from=1,to=2*(M-1)-1, by=2)]
    par.beta=par[seq(from=2,to=2*(M-1), by=2)]
    
    vect.w.times<-matrix(nrow = (M-1), ncol=length(time.event))
    for( i in 1:(M-1)){
      vect.w.times[i,]<-exp(par.alpha[i]+par.beta[i]*time.event)
    }
    if(max(vect.w.times)==Inf){ #pour eviter des poids infinis
      co<-c()
      for (i in 1:(M-1)){
        co<-c(co,which(vect.w.times[i,]==Inf) )
      }
      vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,0)
    }
    
    vect.w.times2<-matrix(nrow = M, ncol=length(time.event))
    for( i in 1:(M-1)){
      vect.w.times2[i,]<-vect.w.times[i,]/(1+apply(vect.w.times,2,sum))
    }
    
    if(M>2){
      vect.w.times2[M,]<-1-apply(vect.w.times2[1:(M-1),],2,sum)
    }
    else{
      vect.w.times2[M,]<-1-vect.w.times2[1,]
    }
    
    
    haz.sl<-array(dim = c(dim(haz[[1]]),M))
    for (i in 1:M){
      haz.sl[,,i]<-haz[[i]]*vect.w.times2[i,]
    }
    haz.sl<-rowSums(haz.sl, dims=2)
    
    
    
    .time=time.event
    
    indicti <- t(sapply(data.times, FUN = function(x) {1*(.time==x)} ))
    indictjsupti<-t(sapply(data.times, FUN = function(x) {1*(.time>=x)} ))
    
    num=haz.sl*indicti
    
    denum=haz.sl*indictjsupti
    
    
    num.s= apply(num, FUN = "sum", MARGIN=1)[data.failures==1]
    denum.s= apply(denum, FUN = "sum", MARGIN=1)[data.failures==1]
    
    num.s<-ifelse(num.s<0,0,num.s)
    denum.s<-ifelse(denum.s<0,0,denum.s)
    
    
    li=num.s/denum.s
    
    logli=sum(log(li))
    
    min.logli=-logli
    return(min.logli)
  }
  
  minus.loglik<- function(par, haz,M, data.times, data.failures,time.event) {
    
    
    # par.alpha=par[seq(from=1,to=2*(M-1)-1, by=2)]
    # par.beta=par[seq(from=2,to=2*(M-1), by=2)]
    # 
    # vect.w.times<-matrix(nrow = (M-1), ncol=length(time.event))
    # for( i in 1:(M-1)){
    #   vect.w.times[i,]<-exp(par.alpha[i]+par.beta[i]*time.event)
    # }
    # if(max(vect.w.times)==Inf){ #pour eviter des poids infinis
    #   co<-c()
    #   for (i in 1:(M-1)){
    #     co<-c(co,which(vect.w.times[i,]==Inf) )
    #   }
    #   vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,0)
    # }
    # 
    # vect.w.times2<-matrix(nrow = M, ncol=length(time.event))
    # for( i in 1:(M-1)){
    #   vect.w.times2[i,]<-vect.w.times[i,]/(1+apply(vect.w.times,2,sum))
    # }
    # 
    # if(M>2){
    #   vect.w.times2[M,]<-1-apply(vect.w.times2[1:(M-1),],2,sum)
    # }
    # else{
    #   vect.w.times2[M,]<-1-vect.w.times2[1,]
    # }
    # 
    .par<-exp(c(par))/(1+sum(exp(par)))
    .par<-c(.par,1-sum(.par))
    
    
    haz.sl<-array(dim = c(dim(haz[[1]]),M))
    for (i in 1:M){
      # haz.sl[,,i]<-haz[[i]]*vect.w.times2[i,]
      haz.sl[,,i]<-haz[[i]]*.par[i]
    }
    haz.sl<-rowSums(haz.sl, dims=2)
    
    
    
    .time=time.event
    
    
    indicti <- t(sapply(data.times, FUN = function(x) {1*(.time==x)} ))
    indictjsupti<-t(sapply(data.times, FUN = function(x) {1*(.time>=x)} ))
    
    num=haz.sl*indicti
    
    denum=haz.sl*indictjsupti
    
    
    num.s= apply(num, FUN = "sum", MARGIN=1)[data.failures==1]
    denum.s= apply(denum, FUN = "sum", MARGIN=1)[data.failures==1]
    
    num.s<-ifelse(num.s<0,0,num.s)
    denum.s<-ifelse(denum.s<0,0,denum.s)
    
    
    li=num.s/denum.s
    
    logli=sum(log(li))
    
    min.logli=-logli
    return(min.logli)
  }
  
  ibs.timedep<-function(par, haz, M,time.event, obj_surv, ot, csurv, csurv_btime, time){
    
    
    par.alpha=par[seq(from=1,to=2*(M-1)-1, by=2)]
    par.beta=par[seq(from=2,to=2*(M-1), by=2)]
    
    vect.w.times<-matrix(nrow = (M-1), ncol=length(time.event))
    for( i in 1:(M-1)){
      vect.w.times[i,]<-exp(par.alpha[i]+par.beta[i]*time.event)
    }
    if(max(vect.w.times)==Inf){ #pour eviter des poids infinis
      co<-c()
      for (i in 1:(M-1)){
        co<-c(co,which(vect.w.times[i,]==Inf) )
      }
      vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,0)
    }
    
    vect.w.times2<-matrix(nrow = M, ncol=length(time.event))
    for( i in 1:(M-1)){
      vect.w.times2[i,]<-vect.w.times[i,]/(1+apply(vect.w.times,2,sum))
    }
    
    if(M>2){
      vect.w.times2[M,]<-1-apply(vect.w.times2[1:(M-1),],2,sum)
    }
    else{
      vect.w.times2[M,]<-1-vect.w.times2[1,]
    }
    
    
    haz.sl<-array(dim = c(dim(haz[[1]]),M))
    for (i in 1:M){
      haz.sl[,,i]<-haz[[i]]*vect.w.times2[i,]
    }
    haz.sl<-rowSums(haz.sl, dims=2)
    
    
    haz.slmono <- apply(haz.sl, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
    
    anti.haz.sl.fct=function(lehaz,vec_l,vec_u){
      anti.haz.sl.f <- function(l, u) {cubature::cubintegrate(lehaz, lower = l, upper = u, method = "pcubature")$integral}
      .i<-mapply(anti.haz.sl.f, vec_l, vec_u)
      return(.i)
    }
    
    .l <- time.event[-length(time.event)]
    .u <- time.event[-1]
    
    # system.time(
    aera<-lapply(haz.slmono, anti.haz.sl.fct, .l,.u)
    # )
    
    Haz.sl=apply(cbind(rep(0,N),do.call(rbind,aera)),cumsum,MARGIN = 1)
    
    .pred=exp(-Haz.sl)
    
    # survs <- t(.pred)[,ot]
    # survs <- t(.pred)[ot,]
    
    survs<-.pred[,ot]
    
    bsc<-sapply(1:length(time.event), FUN = function(j)
    {
      help1 <- as.integer(time <= time.event[j] & obj_surv[ot,2] == 1)
      help2 <- as.integer(time > time.event[j])
      return(mean((0 - survs[j, ])^2 * help1 * (1/csurv) +
                    (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j])))
    })
    
    idx <- 2:length(time.event)
    RET <- diff(time.event) %*% ((bsc[idx - 1] + bsc[idx])/2)
    RET <- RET/diff(range(time.event))
    RET=as.matrix(RET)
    
    return(RET)
  }
  
  ###################################################
  ### Initialisation et recuperation param.tune ###
  ###################################################
  
  # CS on peut avoir plusieurs methodes
  # methods<-unique(methods)
  
  #
  
  
  if(sum(!(methods %in% c("cox.lasso", "cox.ridge", "rf.time", "nn.time", "cox.en",
                          "aft.weibull","aft.weibull","aft.ggamma","aft.gamma","ph.gompertz","ph.exponential",
                          "aft.llogis","cox.aic","cox.all")))>=1){
    stop("New method is not yet implemented")
    # a implementer
    # recuperation code r tune.names_methode
    # recuperation code r fit.names_methode
    # recuperamtion  code r predict.names_methode
    if(verbose==T){print("new method")}
  }
  
  M<-length((methods))
  N <- length(data[,times])
  
  names.meth=c(rep(NA,M))
  for(i in unique(methods)){
    if(length(which(methods==i))==1){
      names.meth[which(methods==i)]=i
    }
    else{
      names.meth[which(methods==i)]=paste0(i,1:length(which(methods==i)))
    }
  }
  
  
  time.pred <- sort(unique(data[,times])) #supprime bien les donnees manquantes
  
  
  # mettre les param.tune par default si manquants 
  # YF : inversion des tests sur param.tune 
  if(is.null(param.tune)==FALSE){ #si certains parametres sont deja rentrees
    if(length(param.tune)!=M){
      stop("Param.tune need to have one element per method. Please modifiy param.tune or set it = NULL")
    }
    for (me in 1:M){
      if(is.null(param.tune[[me]])==T & !(methods[me] %in% c("aft.gamma","aft.ggamma","aft.weibull","ph.exponential","ph.gompoertz","aft.llogis"))){
        if(methods[me] %in%"cox.en"){
          param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
        }
        if(methods[me] %in%"cox.aic"){
          param.tune[[me]]=list(final.model=NA, model.min=NULL, model.max=NULL)
        }
        if(methods[me] %in%"nn.time"){
          param.tune[[me]]=list(n.nodes=c(2, 3, 4, 6, 10, 20),
                                decay=c(0, 0.01, 0.1),
                                batch.size=256L,
                                epochs=1L)
        }
        if(methods[me] %in% "cox.lasso"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"cox.ridge"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"rf.time"){
          param.tune[[me]]=list(mtry=(length(group)+length(cov.quanti)+length(cov.quali))/2+2,
                                nodesize=c(2, 4, 6, 10, 20, 30, 50, 100),
                                ntree=500)
        }
      }
    }
  }
  if(is.null(param.tune)){
    if(verbose==T){print("tune parameters are set to default values for all methods")}
    param.tune=vector("list",M)
    for (me in 1:M){
      if(methods[me] %in%"cox.en"){
        param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
      }
      if(methods[me] %in%"cox.aic"){
        param.tune[[me]]=list(final.model=NA, model.min=NULL,model.max=NULL)
      }
      if(methods[me] %in%"nn.time"){
        param.tune[[me]]=list(n.nodes=c(2, 3, 4, 6, 10, 20),
                              decay=c(0, 0.01, 0.1),
                              batch.size=256L,
                              epochs=1L)
      }
      if(methods[me] %in% "cox.lasso"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in%"cox.ridge"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in%"rf.time"){
        param.tune[[me]]=list(mtry=seq(1,(length(group)+length(cov.quanti)+length(cov.quali))/2+2),
                              nodesize=c(2, 4, 6, 10, 20, 30, 50, 100),
                              ntree=500)
      }
      if (methods[me] %in% c("aft.gamma","aft.ggamma","aft.weibull","ph.exponential","ph.gompoertz","aft.llogis",
                             "cox.all")){
        # param.tune[[me]]=NULL
      }
    }
  }

  
  if(is.null(pro.time)==T & sum(metric %in%c("bs","bll","ribs","ribll","roc"))){
    if(verbose==T){print("pro.time is needed for at least one metric. pro.time was assign to median time value")}
    pro.time=median(data[,times])
  }
  
  
  if(metric%in%c("ibs.td","loglik.td")){
    time.dep.weights<-TRUE
  }
  else{
    time.dep.weights<-FALSE
  }
  
  # FitALL<-vector("list",M)
  # names(FitALL)<-names.meth # cela ne peut plus marche
  
  .model<-vector("list",M)
  names(.model)<-names.meth # cela ne peut plus marche
  
  .tune.optimal<-vector("list",M)
  names(.tune.optimal)<-names.meth
  
  .tune.results<-vector("list",M)
  names(.tune.results)<-names.meth
  
  
  # # CS j'ai supprime car plus besoin
  # # pour que chaque method utilise les memes donnees pour la CV utilisee pour tune.parameter
  # data$id=1:N
  # #set.seed(seed)
  # 
  # sample_id=sample(nrow(data))
  # folds <- cut(seq(1,nrow(data)), breaks=cv, labels=FALSE)
  # folds_id=folds[sample_id]
  # data$folds=folds_id
  
  
  
  ########################
  # STEP 1
  # ########################
  

  
  for (me in 1:M){
    
    if(methods[me] == "aft.weibull" ){
      if(verbose==T){print("start aft.weibull")}
      # .tune.aft.weibull=list()
      # .tune.aft.weibull$optimal=NA
      .tune.optimal[[me]]=NA
      
      .aft.weibull <- aft.weibull(times=times, failures=failures, group=group,
                                  cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.aft.weibull
      # FitALL[[me]]<- .aft.weibull$predictions
      rm(.aft.weibull)
    }
    if(methods[me] == "aft.ggamma"){
      if(verbose==T){print("start aft.ggamma")}
      # .tune.aft.ggamma=list()
      # .tune.aft.ggamma$optimal=NA
      .tune.optimal[[me]]=NA
      
      .aft.ggamma <- aft.ggamma(times=times, failures=failures, group=group,
                                cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.aft.ggamma
      # FitALL[[me]] <- .aft.ggamma$predictions
      rm(.aft.ggamma)
    }
    if(methods[me] == "aft.gamma" ){
      if(verbose==T){print("start aft.gamma")}
      # .tune.aft.gamma=list()
      # .tune.aft.gamma$optimal=NA
      .tune.optimal[[me]]=NA
      
      .aft.gamma <- aft.gamma(times=times, failures=failures, group=group,
                              cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.aft.gamma
      # FitALL[[me]] <- .aft.gamma$predictions
      rm(.aft.gamma)
    }
    if(methods[me] == "aft.llogis" ){
      if(verbose==T){print("start aft.llogis")}
      .tune.optimal[[me]]=NA
      
      .aft.llogis <- aft.llogis(times=times, failures=failures, group=group,
                              cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.aft.llogis
      # FitALL[[me]] <- .aft.gamma$predictions
      rm(.aft.llogis)
    }
    
    if(methods[me] == "ph.gompertz" ){
      if(verbose==T){print("start ph.gompertz")}
      # .tune.ph.gompertz=list()
      # .tune.ph.gompertz$optimal=NA
      .tune.optimal[[me]]=NA
      
      .ph.gompertz <- ph.gompertz(times=times, failures=failures, group=group,
                                  cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.ph.gompertz
      # FitALL[[me]] <- .ph.gompertz$predictions
      rm(.ph.gompertz)
    }
    if(methods[me] == "ph.exponential"){
      if(verbose==T){print("start ph.exponential")}
      # .tune.ph.exponential=list()
      # .tune.ph.exponential$optimal=NA
      .tune.optimal[[me]]=NA
      
      .ph.exponential <- ph.exponential(times=times, failures=failures, group=group,
                                        cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.ph.exponential
      # FitALL[[me]] <- .ph.exponential$predictions
      rm(.ph.exponential)
    }
    
    if(methods[me] == "cox.all"){
      if(verbose==T){print("start cox.all")}
      # .tune.ph.exponential=list()
      # .tune.ph.exponential$optimal=NA
      .tune.optimal[[me]]=NA
      
      .coxall <- cox.all(times=times, failures=failures, group=group,
                                        cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      
      .model[[me]]<-.coxall
      # FitALL[[me]] <- .ph.exponential$predictions
      rm(.coxall)
    }
    

    if(methods[me] == "cox.lasso"){
      if(verbose==T){print("start cox.lasso")}
      
      # soit user ne specifie pas de lambda et le package cherche l'optimal
      # soit user veut chercher parmi un set de lambda
      # if(is.null(param.tune$cox.lasso$lambda)==T | length(param.tune$cox.lasso$lambda)>1){
      # CS param.tune indexe sur position methode
      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){ 
        if(verbose==T){print("estimation of tune parameters for cox.lasso")}
        .tune<- tune.cox.lasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                               cov.quali=cov.quali, data=data, cv=cv,
                               parallel=FALSE, lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]=.tune$optimal
        .tune.results[[me]]=.tune$results
        rm(.tune)
        
      }
      else{ #si on a qu'un seul lambda donc on a pas eu besoin de le chercher
        .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda)
      }
      
      
      .cox.lasso <- cox.lasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=data,
                              lambda=.tune.optimal[[me]]$lambda)
      
      .model[[me]]<-.cox.lasso
      # .model[["cox.lasso"]]$predictions<-NULL
      # FitALL[[me]] <- .cox.lasso$predictions
      rm(.cox.lasso)
    }
    if(methods[me] == "cox.ridge"){
      if(verbose==T){print("start cox.ridge.r")}
      
      # soit user ne specifie pas de lambda et le package cherche l'optimal
      # soit user veut chercher parmi un set de lambda
      # CS idem changement parma.tune
      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){
        if(verbose==T){print("estimation of tune parameters for cox.ridge")}
        .tune<- tune.cox.ridge(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali,
                               data=data, cv=cv,
                               parallel = FALSE, lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda)
      }
      
      .cox.ridge <- cox.ridge(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali,
                              data=data,
                              lambda=.tune.optimal[[me]]$lambda)
      .model[[me]]<-.cox.ridge
      # FitALL[[me]]<- .cox.ridge$predictions
      rm(.cox.ridge)
      
    }
    if(methods[me] == "cox.en"){
      if(verbose==T){print("start cox.en.r")}
      # CS idem param.tune ont ete indexe sur me
      #si un seul alpha et un seul lambda aucune recherche param.tune a faire
      if(length(param.tune[[me]]$alpha)==1 & length(param.tune[[me]]$lambda)==1){
        .tune.optimal[[me]]=list(alpha=param.tune[[me]]$alpha, lambda=param.tune[[me]]$lambda)
      }
      else{
        if(verbose==T){print("estimation of tune parameters for cox.en")}
        .tune <- tune.cox.en(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                             cov.quali=cov.quali, data=data, cv=cv,
                             # foldsid="folds",  # YOHANN # CS
                             parallel=FALSE,
                             alpha=param.tune[[me]]$alpha,
                             lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        
        rm(.tune)
      }
      
      
      
      
      .cox.en <- cox.en(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=data,
                        alpha=.tune.optimal[[me]]$alpha,
                        lambda=.tune.optimal[[me]]$lambda)
      
      .model[[me]]<-.cox.en
      # .model[["cox.en"]]$predictions<-NULL
      # FitALL[[me]] <- .cox.en$predictions
      rm(.cox.en)
    }
    
    if(methods[me] == "cox.aic"){
      if(verbose==T){print("start cox.aic")}
      # CS idem param.tune ont ete indexe sur me
      #si un seul alpha et un seul lambda aucune recherche param.tune a faire
      if(is.na(param.tune[[me]]$final.model)==FALSE){
        .tune.optimal[[me]]=list(final.model=param.tune[[me]]$final.model)
      }
      else{
        if(verbose==T){print("estimation of tune parameters for cox.aic")}
        .tune <- tune.cox.aic(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                             cov.quali=cov.quali, data=data, 
                             # foldsid="folds",  # YOHANN # CS
                             model.min=param.tune[[me]]$model.min,
                             model.max=param.tune[[me]]$model.max)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        
        rm(.tune)
      }
      
      
      
      
      .cox.aic<- cox.aic(times=times, failures=failures, group=group, data=data,
                       cov.quanti=cov.quanti, cov.quali=cov.quali, final.model = .tune.optimal[[me]]$final.model)
      
      .model[[me]]<-.cox.aic
      # .model[["cox.en"]]$predictions<-NULL
      # FitALL[[me]] <- .cox.en$predictions
      rm(.cox.aic)
    }
    
    ### rf.time
    if (methods[me] == "rf.time"){
      if(verbose==T){print("start rf.time")}
      # CS idem param.tune indexe sur me
      if(length(param.tune[[me]]$nodesize)!=1 | length(param.tune[[me]]$mtry)!=1 | length(param.tune[[me]]$ntree)!=1){
        if(verbose==T){print("estimation of tune parameters for rf.time")}
        .tune<-tune.rf.time(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                            cov.quali=cov.quali, data=data, cv=cv,
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
      .rf.time <-rf.time(times=times, failures=failures, 
                         group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=data,
                         nodesize=.tune.optimal[[me]]$nodesize, 
                         mtry=.tune.optimal[[me]]$mtry, ntree=.tune.optimal[[me]]$ntree)
      
      .model[[me]]<-.rf.time
      # .model[["rf.time"]]$predictions<-NULL
      # FitALL[[me]] <- .rf.time$predictions
      rm(.rf.time)
    }
    
    #### nn.time
    if (methods[me] == "nn.time"){
      if(verbose==T){print("start nn.time.r")}
      torch<-reticulate::import("torch")
      torch$set_num_threads(1L)
      # if(parallel_tune==F){
      #   cl=NULL
      # }
      # CS idem param.tune indexe sur me
      if(length(param.tune[[me]]$n.nodes)!=1 | length(param.tune[[me]]$decay)!=1 |
         length(param.tune[[me]]$batch.size)!=1 |length(param.tune[[me]]$epochs)!=1 ){ #si il faut chercher param.tune
        if(verbose==T){print("estimation of tune parameters for nn.time")}
        
        #pour la recherche des parametres tunes xx est une liste avec 3 elements: data.train / data.valid
        # et grid=nodes+decay fixe parmi tous les k folds et les nodes/decays en a tester
        
        .tune<- tune.nn.time(times=times, failures=failures, group=group,
                               cov.quanti=cov.quanti, cov.quali=cov.quali,
                               data=data, cv=cv, 
                               n.nodes=param.tune[[me]]$n.nodes,
                               decay=param.tune[[me]]$decay,
                               batch.size=param.tune[[me]]$batch.size,
                               epochs=param.tune[[me]]$epochs)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        #  parallel=parallel_tune, cl=cl) # YOHANN
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]<-list(n.nodes=param.tune[[me]]$n.nodes,
                                  decay=param.tune[[me]]$decay,
                                  batch.size=param.tune[[me]]$batch.size,
                                  epochs=param.tune[[me]]$epochs)
      }
      
      
      .nn.time <-nn.time(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=data,
                             n.nodes=as.numeric(.tune.optimal[[me]]$n.nodes),
                             decay=as.numeric(.tune.optimal[[me]]$decay),
                             batch.size=as.integer(.tune.optimal[[me]]$batch.size),
                             epochs=as.integer(.tune.optimal[[me]]$epochs))
      
      
      .model[[me]]<-.nn.time
      # .model[["nn.time"]]$predictions<-NULL
      # FitALL[[me]] <- .nn.time$predictions
      # attributes(FitALL$nn.time)=NULL
      rm(.nn.time)
    }
    
  }


  ########################
  ### Cross-Validation  ##
  ########################
  
  if(verbose==T){print("cross-validation initiation")}
  #set.seed(seed)
  sample_id=sample(nrow(data))
  folds <- cut(seq(1,nrow(data)), breaks=cv, labels=FALSE)
  folds_id=folds[sample_id]
  data$folds=folds_id
  
  # sample_id2=sample(nrow(data))
  # data$folds2=folds[sample_id2]
  
  CV<-vector("list",cv*M)
  j<-1
  for(m in 1:M){
    for (k in 1:cv){
      CV[[j]]<-list(train= data[data$folds!=k, ],valid=data[data$folds==k, ], num_method=m)
      j<-j+1
    }
  }
  
  # CS modification de toutes les fonctions pour n'utiliser que predict vu que les modeles ont deja ete fittes
  CV_all_method<-function(CV, method, Tune,
                          times, failures, group, cov.quanti, cov.quali,time.pred){
    num_method<-CV$num_method
    meth<-method[num_method]
    if(meth == "aft.weibull"){
      fit<-aft.weibull(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                       cov.quali=cov.quali, data=CV$train)
      pred=predict(fit,  newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "aft.ggamma"){
      fit<-aft.ggamma(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                      cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "aft.gamma"){
      fit<-aft.gamma(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "aft.llogis"){
      fit<-aft.llogis(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "ph.gompertz"){
      fit<-ph.gompertz(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                       cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "ph.exponential"){
      fit<-ph.exponential(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                          cov.quali=cov.quali, data=CV$train)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "cox.lasso"){
      fit<-cox.lasso(times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                     cov.quali=cov.quali, data=CV$train,
                     lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "cox.en"){
      fit<-cox.en(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali, data=CV$train,
                  alpha=Tune[[num_method]]$alpha,
                  lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth =="cox.ridge"){
      fit<-cox.ridge(times=times, failures=failures, group=group,  cov.quanti=cov.quanti, cov.quali=cov.quali, data=CV$train,
                     lambda=Tune[[num_method]]$lambda)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
      
    }
    if(meth =="cox.aic"){
      fit<-cox.aic(times=times, failures=failures, group=group, data=data,
                   final.model = Tune[[num_method]]$final.model, cov.quanti=cov.quanti, cov.quali=cov.quali)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
      
    }
    if(meth =="cox.all"){
      fit<-cox.all(times=times, failures=failures, group=group,
                   cov.quanti=cov.quanti, cov.quali=cov.quali, data=data)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
      
    }
    if(meth =="rf.time"){
      fit<-rf.time(times=times, failures=failures, group=group, cov.quanti=cov.quanti, cov.quali=cov.quali,  data=CV$train,
                   nodesize=Tune[[num_method]]$nodesize, mtry=Tune[[num_method]]$mtry, 
                   ntree=Tune[[num_method]]$ntree)
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
      
    }
    if(meth =="nn.time"){
      fit<-nn.time(times=times, failures=failures, group=group,  cov.quanti=cov.quanti, cov.quali=cov.quali, data=CV$train,
                     n.nodes=as.numeric(Tune[[num_method]]$n.nodes),
                     decay=as.numeric(Tune[[num_method]]$decay),
                     batch.size=as.integer(Tune[[num_method]]$batch.size),
                     epochs=as.integer(Tune[[num_method]]$epochs))
      pred<-predict(fit,newtimes=time.pred, newdata=CV$valid)$predictions
    }
    return(pred)
  }
  
  
  
  # CV POUR CHAQUE METHODE POUR CHAQUE k en sequentiel / parallel
  if(verbose==T){print("start CV")}
  preFitCV<-lapply(CV, CV_all_method,method=methods,
                   Tune=.tune.optimal, times=times, failures=failures, group=group, cov.quanti=cov.quanti,
                   cov.quali=cov.quali, time.pred=time.pred)
  
  FitCV<-vector("list", M)
  for(m in 1:M){
    FitCV[[m]]<-matrix(nrow=dim(.model[[1]]$predictions)[1], ncol= dim(.model[[1]]$predictions)[2])
    for (k in 1:cv){
      FitCV[[m]][data$folds==k,]<-preFitCV[[(m-1)*cv+k]]
    }
  }
  names(FitCV)<-names.meth
  
  # CS YF pb de memoire
  rm(preFitCV, CV) 
  
  #############################
  # OPTIMISATION
  #########################
  
  data.times <- data[,times]
  data.failures <- data[,failures]
  timeVector=survfit(Surv(data[,times],data[,failures])~ 1 )$time
  
  obj_surv=Surv(data.times, data.failures)
  
  time <- obj_surv[, 1]
  ot <- order(time)
  cens <- obj_surv[ot, 2]
  time <- time[ot]
  
  hatcdist <- prodlim::prodlim(Surv(time, cens) ~ 1, reverse = TRUE)
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf
  csurv_btime <- predict(hatcdist, times = timeVector, type = "surv")
  csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
  csurv_btime[csurv_btime == 0] <- Inf
  
  
  
  # calcul des poids pour chaque metrique
  # CS : ajout d'une option pour optiom et optimParallel qui change la methode d'estimation suivant le nb de poids a estimer
  # ie si 1 poids / 2 learning candidates "BGFS" (pas Brent car il faut preciser lower / upper en +)
  # sinon "Nelder-Mead" 
  .optim.method="Nelder-Mead"
  if(M==2 & time.dep.weights==F){
    .optim.method="BFGS"
  }
  #  CS : MAJ de tous les optim / optimParallel suivant :^
  if(is.null(param.weights.fix)==TRUE){
      # ESTIM<-vector("list",length(metric))
      # w.sl<-vector("list",length(metric))
      # for(me in 1:length(metric)){
      if(verbose==T){print(paste0("start optim for metric ",metric))}
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
             partial.loglik={
               estim<- optim(par=param.weights.init, fn=minus.partial.loglik, FitCV = FitCV, data.times =  data.times,
                             data.failures = data.failures,
                             differentiation = differentiation,hessian=F,
                             method=.optim.method)
               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par, fn=minus.loglik, FitCV = FitCV, data.times =  data.times,
                              data.failures = data.failures,
                              differentiation = differentiation,hessian=F,
                              method=.optim.method)
               }
             }, 
             auc={
               estim<- optim(par=param.weights.init, fn=minus.roc, FitCV = FitCV, data.times =  data.times,
                             data.failures = data.failures,
                             time.pred=time.pred,
                             pro.time=pro.time, precision=ROC.precision,hessian=F,
                             method=.optim.method)
               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par,fn=minus.roc, FitCV = FitCV, data.times =  data.times,
                              data.failures = data.failures,
                              time.pred=time.pred,
                              pro.time=pro.time, precision=ROC.precision,hessian=F,
                              method=.optim.method)
               }
             },
             # CS -> YF ici les trois noubelles fonctions a minimiser
             # loglik et loglik.td 
             loglik={
               
               for(j in 1:M){
                 for (i in 1:dim(FitCV[[1]])[1]){
                   if(min(FitCV[[j]][i,])==0){
                     FitCV[[j]][i,FitCV[[j]][i,]==0]<-min(FitCV[[j]][i,which(FitCV[[j]][i,]!=0)])
                   }
                 }
                 
               }
               
               indic.event.t=time.pred %in% data.times[data.failures==1]
               time.event<-c(0,time.pred[indic.event.t])
               .haz<-vector("list",M)
               
               for (i in 1:M){
                 .Hazcum <- -log(FitCV[[i]][, indic.event.t]) #il faudra verifier que FitCV >0
                 .Hazcum<-cbind(rep(0,N),.Hazcum)
                 # .time <- time.pred[time.pred %in% unique(data.times) ]
                 fmono <- apply(.Hazcum, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
                 Dfmono <- lapply(fmono, FUN=function(fonction){mosaicCalc::D(fonction(t) ~ t)})
                 .haz[[i]] <-t(sapply(1:N,FUN=function(ni){Dfmono[[ni]](time.event)} ))
               }
               for (i in 1:M){
                 if(min(.haz[[i]])<0){
                   .haz[[i]]<-apply(.haz[[i]],2,FUN=function(l){ifelse(l<0,0,l)})
                 }
                 else{
                   .haz[[i]]<-.haz[[i]]
                 }
               }
               
               estim<-optim(par=param.weights.init, fn=minus.loglik, 
                            haz = .haz,
                            M=M, data.times=data.times, 
                            data.failures=data.failures,time.event=time.event,
                            hessian=F,
                            method=.optim.method)
               
               if(optim.local.min==T){
                 start_par=estim$par
                 estim2<-optim(par=start_par, fn=minus.loglik, 
                               haz = .haz,
                               M=M, data.times=data.times, 
                               data.failures=data.failures,time.event=time.event,
                               hessian=F,
                               method=.optim.method)
               }
               
               
             },
             loglik.td={
               
               for(j in 1:M){
                 for (i in 1:dim(FitCV[[1]])[1]){
                   if(min(FitCV[[j]][i,])==0){
                     FitCV[[j]][i,FitCV[[j]][i,]==0]<-min(FitCV[[j]][i,which(FitCV[[j]][i,]!=0)])
                   }
                 }
                 
               }
               
               indic.event.t=time.pred %in% data.times[data.failures==1]
               time.event<-c(0,time.pred[indic.event.t])
               .haz<-vector("list",M)
               
               for (i in 1:M){
                 .Hazcum <- -log(FitCV[[i]][, indic.event.t]) #il faudra verifier que FitCV >0
                 .Hazcum<-cbind(rep(0,N),.Hazcum)
                 # .time <- time.pred[time.pred %in% unique(data.times) ]
                 fmono <- apply(.Hazcum, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
                 Dfmono <- lapply(fmono, FUN=function(fonction){mosaicCalc::D(fonction(t) ~ t)})
                 .haz[[i]] <-t(sapply(1:N,FUN=function(ni){Dfmono[[ni]](time.event)} ))
               }
               for (i in 1:M){
                 if(min(.haz[[i]])<0){
                   .haz[[i]]<-apply(.haz[[i]],2,FUN=function(l){ifelse(l<0,0,l)})
                 }
                 else{
                   .haz[[i]]<-.haz[[i]]
                 }
               }
               
               estim<-optim(par=param.weights.init, fn=minus.loglik.timedep, 
                            haz = .haz,
                            M=M, data.times=data.times, 
                            data.failures=data.failures,time.event=time.event,
                            hessian=F,
                            method=.optim.method)
               
               if(optim.local.min==T){
                 start_par=estim$par
                 estim2<-optim(par=start_par, fn=minus.loglik.timedep, 
                               haz = .haz,
                               M=M, data.times=data.times, 
                               data.failures=data.failures,time.event=time.event,
                               hessian=F,
                               method=.optim.method)
               }
               
               
             },
             ibs.td={
               # if(is.null(pro.time)){
               #   pro.time=median(data[,times])
               # }
               time.pred <- sort(unique(data.times))
               
               for(j in 1:M){
                 for (i in 1:dim(FitCV[[1]])[1]){
                   if(min(FitCV[[j]][i,])==0){
                     FitCV[[j]][i,FitCV[[j]][i,]==0]<-min(FitCV[[j]][i,which(FitCV[[j]][i,]!=0)])
                   }
                 }
                 
               }
               
               indic.event.t=time.pred %in% data.times[data.failures==1]
               time.event<-c(0,time.pred[indic.event.t])
               .haz<-vector("list",M)
               
               for (i in 1:M){
                 .Hazcum <- -log(FitCV[[i]][, indic.event.t]) #il faudra verifier que FitCV >0
                 .Hazcum<-cbind(rep(0,N),.Hazcum)
                 # .time <- time.pred[time.pred %in% unique(data.times) ]
                 fmono <- apply(.Hazcum, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
                 Dfmono <- lapply(fmono, FUN=function(fonction){mosaicCalc::D(fonction(t) ~ t)})
                 .haz[[i]] <-t(sapply(1:N,FUN=function(ni){Dfmono[[ni]](time.event)} ))
               }
               for (i in 1:M){
                 if(min(.haz[[i]])<0){
                   .haz[[i]]<-apply(.haz[[i]],2,FUN=function(l){ifelse(l<0,0,l)})
                 }
                 else{
                   .haz[[i]]<-.haz[[i]]
                 }
               }
               
               estim<-optim(par=param.weights.init, fn=ibs.timedep, haz=.haz, M=M,time.event=time.event,
                            obj_surv=obj_surv, ot=ot, csurv=csurv,
                            csurv_btime=csurv_btime,time=time,
                            hessian=F,
                            method=.optim.method)
               
               if(optim.local.min==T){
                 start_par=estim$par
                 estim<-optim(par=start_par,  fn=ibs.timedep, haz=.haz, M=M,time.event=time.event,
                              obj_surv=obj_surv, ot=ot, csurv=csurv,
                              csurv_btime=csurv_btime,time=time,
                              hessian=F,
                              method=.optim.method)
               }
               
             }
      )
  }
  

  ######################################
  # Calcul Survie SL
  ###########
  
  if(verbose==T){print("pool results")}
  
  if(is.null(param.weights.fix)==FALSE){
    estim=list()
    estim$par=param.weights.fix
  }
  
  # CS YF pb de memoire
  FitALL<-vector("list",M)
  names(FitALL)<-names.meth # cela ne peut plus marche
  
  for(me in 1:M){
    FitALL[[me]]<-.model[[me]]$predictions
    
  }
  
  if(time.dep.weights==T){
    w.sl <- estim$par
    time.pred <- sort(unique(data.times))
    Mm1=M-1
    
    par.alpha=estim$par[seq(from=1,to=2*Mm1-1, by=2)]
    par.beta=estim$par[seq(from=2,to=2*Mm1, by=2)]
    
    vect.w.times<-matrix(nrow = Mm1, ncol=length(data.times))
    for( i in 1:Mm1){
      vect.w.times[i,]<-exp(par.alpha[i]+par.beta[i]*data.times)
    }
    if(max(vect.w.times)==Inf){ #pour eviter des poids infinis
      co<-c()
      for (i in 1:(M-1)){
        co<-c(co,which(vect.w.times[i,]==Inf) )
      }
      vect.w.times[,co]<-ifelse(vect.w.times[,co]==Inf,10**10,0)
    }
    
    vect.w.times2<-matrix(nrow = M, ncol=length(data.times))
    for( i in 1:Mm1){
      vect.w.times2[i,]<-vect.w.times[i,]/(1+apply(vect.w.times,2,sum))
    }
    if(M>2){
      vect.w.times2[M,]<-1-apply(vect.w.times2[1:Mm1,],2,sum)
    }
    else{
      vect.w.times2[M,]<-1-vect.w.times2[1,]
    }
    
    # View(vect.w.times2)
    indic.event.t=time.pred %in% data.times[data.failures==1]
    time.event<-c(0,time.pred[indic.event.t])
    
    for(j in 1:M){
      for (i in 1:dim(.model[[1]]$predictions)[1]){
        if(min(.model[[j]]$predictions[i,])==0){
          .model[[j]]$predictions[i,.model[[j]]$predictions[i,]==0]<-min(.model[[j]]$predictions[i,which(.model[[j]]$predictions[i,]!=0)])
        }
      }
      
    }
    
    .haz<-vector("list",M)
    for (i in 1:M){
      .Hazcum <- -log(.model[[i]]$predictions[, indic.event.t]) #il faudra verifier que FitCV >0
      .Hazcum<-cbind(rep(0,N),.Hazcum)
      # .time <- time.pred[time.pred %in% unique(data.times) ]
      fmono <- apply(.Hazcum, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
      Dfmono <- lapply(fmono, FUN=function(fonction){mosaicCalc::D(fonction(t) ~ t)})
      .haz[[i]] <-t(sapply(1:N,FUN=function(ni){Dfmono[[ni]](time.event)} ))
    }
    for (i in 1:M){
      if(min(.haz[[i]])<0){
        .haz[[i]]<-apply(.haz[[i]],2,FUN=function(l){ifelse(l<0,0,l)})
      }
      else{
        .haz[[i]]<-.haz[[i]]
      }
    }
    
    haz.sl<-array(dim = c(dim(.haz[[1]]),M))
    for (i in 1:M){
      haz.sl[,,i]<-.haz[[i]]*vect.w.times2[i,]
    }
    haz.sl<-rowSums(haz.sl, dims=2)
    
    
    # Haz.sl=t(apply(haz.sl, MARGIN=1,integration, x=time.pred))
    
    haz.slmono <- apply(haz.sl, FUN=function(H){mosaic::spliner(H ~ time.event, monotonic = TRUE)}, MARGIN = 1)
    
    anti.haz.sl.fct=function(lehaz,vec_l,vec_u){
      anti.haz.sl.f <- function(l, u) {cubature::cubintegrate(lehaz, lower = l, upper = u, method = "pcubature")$integral}
      .i<-mapply(anti.haz.sl.f, vec_l, vec_u)
      return(.i)
    }
    
    .l <- time.event[-length(time.event)]
    .u <- time.event[-1]
    
    # system.time(
    aera<-lapply(haz.slmono, anti.haz.sl.fct, .l,.u)
    # )
    
    Haz.sl=apply(cbind(rep(0,N),do.call(rbind,aera)),cumsum,MARGIN = 1)
    
    .pred=exp(-Haz.sl)
    survs<-t(.pred)
    # CS : MAJ pour eviter la boucle
    idx=findInterval(timeVector, time.event)
    survs.long=survs[,pmax(1,idx)]
    
    
    surv.SL=survs.long
    
  }
  if(time.dep.weights==F){
    w.sl <- c(exp(c(estim$par,0)) / ( 1+sum(exp(estim$par))) )
    .SL<-array(dim = c(dim(FitALL[[1]]),length(FitALL)))
    for (i in 1:length(FitCV)){
      .SL[,,i]<-.model[[i]]$predictions*w.sl[i]
    }
    surv.SL <-rowSums(.SL, dims=2)
    # names(w.sl) <- names.meth
  }

  ##################################################
  # PREPARATION RETURN
  ########################################
  
  # Il s'agit ici de la fin de la fonction d'apprentissage (par ex: superlearner.surv ),
  # il faut y placer un return d'une liste qui possedent a minima les arguments suivant :
  # method, chacun des modeles estimes (unr liste), w.sl, FitALL, et time.pred, en renomant
  # des objets comprehansible et en coherence avec les sortie des fonctions cox.lasso et autres.
  # Proposer les phrases de l'aide qui seront dans le package. Comme pour les fonctions cox.lasso
  # et autres, il faudrait prevoir deux arguments supplementaires newdata et newtimes.
  # Il faut aussi prevoir une fonction "predict" qui appel l'objet stocke a la sortie de la fonction
  # d'apprentissage, par exemple predict.superlearner.surv, avec ces deux options newdata et newtimes.
  
  
  if(keep.predictions==TRUE) {
    FitALL<-vector("list",M)
    names(FitALL)<-names.meth # cela ne peut plus marche
    
    for(me in 1:M){
      FitALL[[me]]<-.model[[me]]$predictions
      
    }
    FitALL$sl<-surv.SL
    temp.predictions <- FitALL}
  else {
    temp.predictions <- surv.SL
    temp.predictions<-as.data.frame(temp.predictions)
  }
  
  res<-list(times=time.pred,#
            predictions=temp.predictions,
            data=data.frame(times=data[,times], failures=data[,failures], data[, !(dimnames(data)[[2]] %in% c(times, failures))]),
            outcomes=list(times=times, failures=failures),
            predictors=list(group=group, cov.quanti=cov.quanti, cov.quali=cov.quali), #
            ROC.precision=ROC.precision, #
            cv=cv, #
            pro.time=pro.time, #
            methods=methods,
            models=.model,
            weights=list(coefficients=estim$par, values=w.sl), #
            metric=list(metric=metric, value=estim$value),
            param.tune=list(optimal=.tune.optimal, results=.tune.results))
  
  class(res) <- "sl.time"
  
  if(verbose==T){print("SL is fitted")}
  return(res)
}



