
summary.sl.time <- function(object, ..., digits=7, method="sl", pro.time=NULL, newdata=NULL, times=NULL, failures=NULL)
{
  if(is.null(pro.time)) {pro.time <- median(object$data$times)}
  
  if(is.null(newdata))
  {
  return(
    round(  data.frame(
    auc = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                 metric="auc", pro.time=pro.time),
    bs = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                metric="bs", pro.time=pro.time),
    ibs = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                 metric="ibs", pro.time=pro.time),
    ribs = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                  metric="ribs", pro.time=pro.time),
    bll = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                 metric="bll", pro.time=pro.time),
    ibll = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                  metric="ibll", pro.time=pro.time),
    ribll = metric(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions[[method]], prediction.times=object$times,
                   metric="ribll", pro.time=pro.time) ), digits = digits ) )
  }
  
  else
  {
    .pred <- predict(object, newdata=newdata)
    if(is.null(times)) {times <- object$outcomes$times; failures <- object$outcomes$failures}
    
    return(
      round(  data.frame(
      auc = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="auc", pro.time=pro.time),
      bs = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times,  metric="bs", pro.time=pro.time),
      ibs = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="ibs", pro.time=pro.time),
      ribs = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="ribs", pro.time=pro.time),
      bll = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="bll", pro.time=pro.time),
      ibll = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="ibll", pro.time=pro.time),
      ribll = metric(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions[[method]], prediction.times=object$times, metric="ribll", pro.time=pro.time) ), digits = digits ) )
  }
}
