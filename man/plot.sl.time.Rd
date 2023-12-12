\name{plot.sl.time}
\alias{plot.sl.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration Plot for Super Learner
}
\description{
A calibration plot of a Super Learner obtained by the function \code{sl.time}.
}

\usage{
\method{plot}{sl.time}(x, ..., col, lty, lwd, type, pch, ylab, xlab, ylim, xlim, cex, cex.lab, 
   cex.axis, cex.main, method, n.groups, pro.time, newdata, times, failures)
}

\arguments{
  \item{x}{An object returned by the function \code{sl.time}.}
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
  \item{method}{A character string with the name of the algorithm included in the SL for wich the calibration plot is performed. The defaut is "sl" for the Super Learner.}
  \item{n.groups}{A numeric value with the number of groups by their class probabilities. The default is 5.}
  \item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the calibration plot is performed from the same subjects of the training sample.}
  \item{times}{The name of the variable related the numeric vector with the follow-up times in \code{newdata} (optional argument only necessary when newdata is not \code{NULL}).}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event) in \code{newdata} (optional argument only necessary when newdata is not \code{NULL}).}
}


\value{No return value for this S3 method.}

\examples{
data(dataDIVAT2)

#The outcome model base on a Super Learner from the first 150 individuals of the data base
sl1<-sl.time( methods=c("aft.gamma", "ph.gompertz"),  metric="ibs",
  data=dataDIVAT2[1:150,],  times="times", failures="failures", group="ecd",
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant"), cv=3)

# The calibration plot from the validation sample of 150 patients
plot(sl1, method="sl", n.groups=5, pro.time=12, col=2,
     xlab="Predicted 12-year survival", ylab="Observed 12-year survival",
     newdata=dataDIVAT2[151:300,], times="times", failures="failures")
}

\keyword{Plotting}
