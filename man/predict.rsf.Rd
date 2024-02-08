\name{predict.rsf}
\alias{predict.rsf}

\title{
Prediction from a Survival Random Survival Forest
}
\description{
Predict the survival of new observations based on a random forest tree  estimated by using \code{RSF}.
}

\usage{
\method{predict}{rsf}(object, ..., newdata, newtimes)
}

\arguments{
  \item{object}{An object returned by the function \code{RSF}.}
  \item{...}{Further arguments passed.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in \code{cov.quanti} and \code{cov.quali} included in the training sample. The default value is \code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the observed times in the training data frame.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}


\seealso{
  \code{\link{RSF}}.
}

\examples{
data(dataDIVAT2)

# The estimation of the training model
model <- RSF(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), nodesize=10,
  mtry=2, ntree=100)

# Predicted survival from the validation sample
pred <- predict(model,
  newdata=data.frame(age=c(52,52), hla=c(0,1), retransplant=c(1,1), ecd=c(0,1)))

plot(y=pred$predictions[1,], x=pred$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

lines(y=pred$predictions[2,], x=pred$times, col=2, type="l", lty=1, lwd=2)

legend("bottomright", col=c(1,2), lty=1, lwd=2, c("Subject #1", "Subject #2"))
}

\keyword{Predicting}
