\name{LIB_COXaic}
\alias{LIB_COXaic}


\title{
Library of the Super Learner for a Cox Model with Selected Covariates
}

\description{
Fit a  Cox regression for a selection of covariate.
}


\usage{
LIB_COXaic(formula, data, penalty=NULL)
}

\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time, the event,and the other covariables included in the formula.}
  \item{penalty}{A numerical vector that allows the integration of covariates into the final model after selection. We give the value 0 if we want to force the covariate in the model otherwise 1. If \code{NULL}, all covariates undergo the selection process. }
 }



\value{
\item{formula}{The formula object used for model construction.}
\item{model}{The estimated model.}
\item{data}{The data frame used for learning.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}


\references{
Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, https://www.jstatsoft.org/v39/i05/
}

\examples{
data("dataDIVAT2")


# The estimation of the model from the first 200 lines

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_COXaic(formula=formula, data=dataDIVAT2[1:200,])

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}

