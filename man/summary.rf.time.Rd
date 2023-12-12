\name{summary.rf.time}
\alias{summary.rf.time}

\title{
Summaries of Survival Random Forest
}
\description{
Return goodness-of-fit indicators of Survival Random Forest Tree obtained by the function \code{rf.time}.
}

\usage{
\method{summary}{rf.time}(object, ..., digits, pro.time, newdata, times, failures)
}

\arguments{
  \item{object}{An object returned by the function \code{rf.time}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{pro.time}{The prognostic time up to which the time-dependent indicators are estimated.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the calibration plot is performed from the same subjects of the training sample.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times in \code{newdata}. If \code{NULL}, the default is the name used in the learning sample.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event) in \code{newdata}. If \code{NULL}, the default is the name used in the learning sample.}
}


\value{No return value for this S3 method.}


\details{
The following metrics are returned: "brier" for the Brier score at the prognostic time \code{pro.time}, "ibs" for the Integrated Brier score up to the last observed time of event, "ibll" for the Integrated Binomial Log-likelihood up to the last observed time of event, "bll" for the binomial Log-likelihood, "ribs" for the restricted Integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted Integrated Binomial Log-likelihood Log-likelihood up to the last observed time of event, "bll" for the binomial Log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}

\seealso{
  \code{\link{rf.time}}.
}

\examples{
data(dataDIVAT2)

# The training of the random forest with the first 150 patients
model<-rf.time(times="times", failures="failures", data=dataDIVAT2[1:150,],
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), nodesize=10,
  mtry=2, ntree=100)
  
# The prognostic capacities from the same training sample
# (up to 4 years forseveral indicators)
summary(model, pro.time=4)

# The prognostic capacities from a validation of the next 150 patients
# (up to 4 years forseveral indicators)
summary(model, pro.time=4, newdata=dataDIVAT2[151:300,], times="times",
failures="failures")
}

\keyword{Summarizing}
