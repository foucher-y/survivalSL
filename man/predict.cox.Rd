\name{predict.cox}
\alias{predict.cox}

\title{
Prediction from a Penalized Cox Regression
}
\description{
Predict the survival of new observations based on a penalized Cox regression estimated by using a model of the class \code{cox}.
}

\usage{
\method{predict}{cox}(object, ..., newdata, newtimes)
}

\arguments{
  \item{object}{An object returned by one of the following functions: \code{COXlasso}, \code{COXridge}, or \code{COXen}.}
  \item{...}{Further arguments passed.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in \code{cov.quanti} and \code{cov.quali} included in the training sample. The default value is \code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the observed times in the training data frame.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\seealso{
\code{\link{COXall}}, \code{\link{COXen}}, \code{\link{COXlasso}}, \code{\link{COXridge}}.
}

\examples{
data(dataDIVAT2)

# The estimation of the training model
model<-COXlasso(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), lambda=.01)

# Predicted survival from the validation sample
pred <- predict(model,
  newdata=data.frame(age=c(52,52), hla=c(0,1), retransplant=c(1,1), ecd=c(0,1)))

plot(y=pred$predictions[1,], x=pred$times, xlab="Time (years)", ylab="Predicted survival",
     col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

lines(y=pred$predictions[2,], x=pred$times, col=2, type="l", lty=1, lwd=2)

legend("bottomright", col=c(1,2), lty=1, lwd=2, c("Subject #1", "Subject #2"))
}

\keyword{Predicting}
