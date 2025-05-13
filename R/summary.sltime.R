summary.sltime <- function(object, newdata=NULL, method="sl",
                           ROC.precision=seq(.01,.99,.01), digits=7,pro.time=NULL, ...) {


  formula<-object$formula
  variables_formula <- all.vars(formula)
  times <- variables_formula[1]
  failures <- variables_formula[2]

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



  if(is.null(pro.time)) {
    pro.time <- median((object$data)[[times]][-1])
  }

  time.pred <- unique(sort(c(0,pro.time,object$times)))



  if(is.null(newdata)){


    time.pred <- unique(sort(c(0,pro.time,object$times)))

    survivals.matrix <- predict(object, newdata=object$data, newtimes=time.pred)$predictions[[method]]

    if(!(method %in% c("sl","LIB_PLANN","LIB_SNN","LIB_RSF","LIB_COXall","LIB_COXaic","LIB_COXlasso","LIB_COXen","LIB_COXridge"))){
      .flex<-object$models[[method]]$model
      .data<-object$models[[method]]$data
      .time<-time.pred
      .hazlist<-predict(
        .flex,
        newdata=.data,
        type = "haz",
        times = .time
      )
      hazards.matrix <- t(sapply(.hazlist$.pred, function(x) x[[2]]))
    }else{
      hazards.matrix<-t(apply(survivals.matrix[,-1],1,haz_function,times=time.pred[-1]))
      }

    return(
      list(metrics=round(  data.frame(
        p_ci = metrics(metric="p_ci", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        uno_ci = metrics(metric="uno_ci", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        auc = metrics(metric="auc", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bs = metrics(metric="bs", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                     hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibs = metrics(metric="ibs", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribs = metrics(metric="ribs", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bll = metrics(metric="bll", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibll = metrics(metric="ibll", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribll = metrics(metric="ribll", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                        hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ll= metrics(metric="ll", times=times, failures=failures, data=object$data, survivals.matrix=survivals.matrix,
                    hazards.matrix=hazards.matrix,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision)), digits = digits ),
        method=method,
        pro.time=pro.time,
        ROC.precision=ROC.precision) )




  }else{

    time<-unique(newdata[[times]])

    time.pred<-sort(c(time,time.pred))

    survivals.matrix <- predict(object, newdata=newdata, newtimes=time.pred)$predictions[[method]]

    if(!(method %in% c("sl","LIB_PLANN","LIB_SNN","LIB_RSF","LIB_COXall","LIB_COXaic","LIB_COXlasso","LIB_COXen","LIB_COXridge"))){
      .flex<-object$models[[method]]$model
      .time<-time.pred
      .hazlist<-predict(
        .flex,
        newdata=newdata,
        type = "haz",
        times = .time
      )
      hazards.matrix <- t(sapply(.hazlist$.pred, function(x) x[[2]]))
      }else{
      hazards.matrix<-t(apply(survivals.matrix[,-1],1,haz_function,times=time.pred[-1]))

    }


    return(
      list(metrics=round(  data.frame(
        p_ci = metrics(metric="p_ci", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        uno_ci = metrics(metric="uno_ci", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        auc = metrics(metric="auc", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bs = metrics(metric="bs", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                     hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibs = metrics(metric="ibs", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribs = metrics(metric="ribs", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bll = metrics(metric="bll", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibll = metrics(metric="ibll", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribll = metrics(metric="ribll", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                        hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ll= metrics(metric="ll", times=times, failures=failures, data=newdata, survivals.matrix=survivals.matrix,
                    hazards.matrix=hazards.matrix,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision)), digits = digits ),
        method=method,
        pro.time=pro.time,
        ROC.precision=ROC.precision) )




  }



}





