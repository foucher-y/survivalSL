\name{LIB_COXen}
\alias{LIB_COXen}

\title{
Library of the Super Learner for Elastic Net Cox Regression
}

\description{
Fit an elastic net Cox regression for fixed values of the regularization parameters.
}

\usage{
LIB_COXen(formula,
data, penalty=NULL, alpha, lambda)
}


\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time, the event,and the other covariables included in the formula.}
 \item{penalty}{A numerical vector that allows the covariates not to be penalized. We give the value 0 if we do not want the covariate to be penalized otherwise 1. If \code{NULL}, all covariates are penalized.}
  \item{alpha}{The value of the regularization parameter alpha for penalizing the partial likelihood.}
  \item{lambda}{The value of the regularization parameter lambda for penalizing the partial likelihood.}
 }



\details{
The elastic net Cox regression is obtained by using the \code{glmnet} package.
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


# The estimation of the model from the first 200 lignes

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_COXen(formula=formula, data=dataDIVAT2[1:200,], lambda=.1, alpha=.1)

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}

