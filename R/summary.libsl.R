
summary.libsl <- function(object, ..., digits=7,
                      pro.time=NULL, newdata=NULL, times=NULL, failures=NULL, ROC.precision=seq(.01,.99,.01))
{
  if(is.null(pro.time)) {pro.time <- median(object$data$times)}

  if(is.null(newdata))
  {
  return(
    round(  data.frame(
    ci = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                prediction.times=object$times, metric="ci", pro.time=pro.time),
    auc = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                 prediction.times=object$times, metric="auc", pro.time=pro.time, ROC.precision=ROC.precision),
    bs = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                prediction.times=object$times,  metric="bs", pro.time=pro.time),
    ibs = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                 prediction.times=object$times,  metric="ibs", pro.time=pro.time),
    ribs = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                  prediction.times=object$times, metric="ribs", pro.time=pro.time),
    bll = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                 prediction.times=object$times, metric="bll", pro.time=pro.time),
    ibll = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                  prediction.times=object$times, metric="ibll", pro.time=pro.time),
    ribll = metrics(times="times", failures="failures", data=object$data, prediction.matrix=object$predictions,
                   prediction.times=object$times, metric="ribll", pro.time=pro.time) ), digits = digits ) )
  }

  else {
    .pred <- predict(object, newdata=newdata)
    if(is.null(times)) {times <- object$outcomes$times; failures <- object$outcomes$failures}

    return(
      round(  data.frame(
      ci = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                  prediction.times=object$times, metric="ci", pro.time=pro.time),
      auc = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                   prediction.times=object$times, metric="auc", pro.time=pro.time, ROC.precision=ROC.precision),
      bs = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                  prediction.times=object$times, metric="bs", pro.time=pro.time),
      ibs = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                   prediction.times=object$times, metric="ibs", pro.time=pro.time),
      ribs = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                    prediction.times=object$times, metric="ribs", pro.time=pro.time),
      bll = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                   prediction.times=object$times, metric="bll", pro.time=pro.time),
      ibll = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                    prediction.times=object$times, metric="ibll", pro.time=pro.time),
      ribll = metrics(times=times, failures=failures, data=newdata, prediction.matrix=.pred$predictions,
                     prediction.times=object$times, metric="ribll", pro.time=pro.time) ), digits = digits ) )
  }
}

