\name{summary.libsl}
\alias{summary.libsl}

\title{
Summaries of a Parametric Survival Model
}
\description{
Return predictive performences of a model or algorithm obtained by a library  of the class \code{libsl}.
}

\usage{
\method{summary}{libsl}(object, ..., digits=7, pro.time=NULL, newdata=NULL,
times=NULL, failures=NULL, ROC.precision=seq(.01,.99,.01))
}

\arguments{
  \item{object}{An object returned by a library of the class \code{libsl}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{pro.time}{The prognostic time up to which the time-dependent indicators are estimated.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the calibration plot is performed from the same subjects of the training sample.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times in \code{newdata}. If \code{NULL}, the default is the name used in the learning sample.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event) in \code{newdata}. If \code{NULL}, the default is the name used in the learning sample.}
  \item{ROC.precision}{An optional argument with the percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. 0 (min) and 1 (max) are not allowed. By default, the precision is \code{seq(.01,.99,.01)}.}
}


\value{No return value for this S3 method.}

\details{
The following metrics are returned: "brier" for the Brier score at the prognostic time \code{pro.time}, "ibs" for the Integrated Brier score up to the last observed time of event, "ibll" for the Integrated Binomial Log-likelihood up to the last observed time of event, "bll" for the binomial Log-likelihood, "ribs" for the restricted Integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted Integrated Binomial Log-likelihood Log-likelihood up to the last observed time of event, "bll" for the binomial Log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}

\seealso{
  \code{\link{LIB_AFTgamma}}, \code{\link{LIB_AFTggamma}}, \code{\link{LIB_AFTllogis}}, \code{\link{LIB_AFTweibull}}, \code{\link{LIB_PHexponential}}, \code{\link{LIB_PHgompertz}}.
}

\examples{
data(dataDIVAT2)

# The training of the Weibull model with the first 400 patients
model <- LIB_PHgompertz(times="times", failures="failures", data=dataDIVAT2[1:400,],
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"))

# The prognostic capacities from the same training sample
# (up to 4 years forseveral indicators)
summary(model, pro.time=4)

# The prognostic capacities from a validation of the next 150 patients
# (up to 4 years forseveral indicators)
 summary(model, pro.time=4, newdata=dataDIVAT2[401:550,], times="times",
 failures="failures")
}

\keyword{Summarizing}