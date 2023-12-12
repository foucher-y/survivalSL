\name{summary.sl.time}
\alias{summary.sl.time}

\title{
Summaries of a Super Learner
}
\description{
Return goodness-of-fit indicators of a Super Learner obtained by the function \code{sl.time}.
}

\usage{
\method{summary}{sl.time}(object, ..., digits, method, pro.time, newdata, times, failures)
}

\arguments{
  \item{object}{An object returned by the function \code{sl.time}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{method}{A character string with the name of the algorithm included in the SL for wich the calibration plot is performed. The defaut is "sl" for the Super Learner.}
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
  \code{\link{sl.time}}.
}

\examples{
data(dataDIVAT2)
dataDIVAT2$train <- 1*rbinom(n=dim(dataDIVAT2)[1], size = 1, prob=1/2)

# The training of the super learner with 2 algorithms from the
   # first 100 patients of the training sample
sl<-sl.time(method=c("aft.gamma", "ph.gompertz"),  metric="ibs",
  data=dataDIVAT2[dataDIVAT2$train==1,][1:100,],  times="times", failures="failures",
  pro.time = 12,  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
  cv=3)

# The prognostic capacities from the same training sample
summary(sl)
}

\keyword{Summarizing}
