\name{plot.rsf}
\alias{plot.rsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration Plot for a Survival Random Forest
}
\description{
A calibration plot of an object of the class \code{rsf}.
}

\usage{
\method{plot}{rsf}(x, ..., col, lty, lwd, type, pch, ylab, xlab, ylim, xlim,
   cex, cex.lab, cex.axis, cex.main, n.groups, pro.time, newdata, times, failures)
}

\arguments{
  \item{x}{An object returned by the function \code{RSF}.}
  \item{...}{Additional arguments affecting the plot.}
  \item{col}{A numeric value with the color of the survival curves. The default is 1 for black.}
  \item{lty}{A numeric value with the type of the survival curves. The default is 1.}
  \item{lwd}{A numeric value with the type of the survival curves. The default is 1.}
  \item{type}{A character string giving the type of plot desired. The default is "b" for both the lines and points}
  \item{pch}{A numeric value for plotting characters or symbols. The default is 16.}
  \item{xlab}{The label for the x axis. The default is "Predicted survival".}
  \item{ylab}{The label for the y axis. The default is "Observed survival".}
  \item{cex}{A number indicating the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50\% larger, etc.}
  \item{cex.lab}{The size of the axis label text with a numeric value of length 1.}
  \item{cex.axis}{The size of the tick label numbers/text with a numeric value of length 1.}
  \item{cex.main}{The size of the title text with a numeric value of length 1.}
  \item{xlim}{A vector with the x limits of the plot.}
  \item{ylim}{A vector with the y limits of the plot.}
  \item{n.groups}{A numeric value with the number of groups by their class probabilities. The default is 5.}
  \item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the calibration plot is performed from the same subjects of the training sample.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times in \code{newdata} (optional argument only necessary when newdata is not \code{NULL}).}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event) in \code{newdata} (optional argument only necessary when newdata is not \code{NULL}).}
}

\details{
The plot represents the observed survival and the related 95\% confidence intervals, which are respectively estimated by the Kaplan and Meier estimator and the Greenwood formula, against the mean of the predictive values for individuals stratified into groups of the same size according to the percentiles. The identity line is usually included for reference.
}

\value{No return value for this S3 method.}

\seealso{
  \code{\link{RSF}}.
}

\examples{
data(dataDIVAT2)
# The estimation of the model from the first 200 lignes
model<-RSF(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), nodesize=10,
  mtry=2, ntree=100)

# The calibration plot from the validation sample of 150 patients
plot(model, n.groups=5, pro.time=12, col=2,
     xlab="Predicted 12-year survival", ylab="Observed 12-year survival",
     newdata=dataDIVAT2[151:300,], times="times", failures="failures")
}

\keyword{Plotting}
