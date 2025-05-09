\name{predict.sltime}
\alias{predict.sltime}

\title{
Prediction from a Super Learner for Censored Outcomes
}
\description{
Predict the survival of new observations based on an SL by using the \code{survivalSL} function.
}

\usage{
\method{predict}{sltime}(object, newdata, newtimes, ...)
}

\arguments{
  \item{object}{An object returned by the function \code{survivalSL}.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the observed times in the training data frame.}
  \item{...}{For future methods.}
}

\value{
\item{predictions}{A list of matrix with the predictions of survivals of each subject (lines) for each observed time (columns) for each model used for the superlearner construction and the superlearner itself.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
}


\seealso{
  \code{\link{survivalSL}}.
}

\examples{
data("dataDIVAT2")

# The training of the super learner from the first 150 individuals of the data base


formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
sl1 <- survivalSL(formula, method=c("LIB_COXridge", "LIB_AFTggamma"),
metric="auc", data=dataDIVAT2[1:150,], pro.time = 12, cv=3)

# Individual prediction for 2 new subjects
pred <- predict(sl1,
  newdata=data.frame(age=c(52,52),
  hla=c(0,1), retransplant=c(1,1), ecd=c(0,1)))

plot(y=pred$predictions$sl[1,], x=pred$times,
xlab="Time (years)",
ylab="Predicted survival",
col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

lines(y=pred$predictions$sl[2,], x=pred$times,
col=2, type="l", lty=1, lwd=2)

legend("bottomright", col=c(1,2), lty=1, lwd=2,
c("Subject #1", "Subject #2"))
}

\keyword{Predicting}

