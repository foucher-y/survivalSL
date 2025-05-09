\name{tuneCOXridge}
\alias{tuneCOXridge}

\title{
Tune Ridge Cox Regression
}

\description{
This function finds the optimal lambda parameter for a ridge Cox regression.
}

\usage{
tuneCOXridge(formula, data, penalty=NULL,
 cv=10, parallel=FALSE, lambda, seed=NULL)
}

\arguments{
   \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model with the same covariates as in the formula.}
  \item{penalty}{A numerical vector that allows the covariates not to be penalized. We give the value 0 if we do not want the covariate to be penalized otherwise 1. If \code{NULL}, all covariates are penalized.}
  \item{cv}{The value of the number of folds. The default value is 10.}
  \item{parallel}{If \code{TRUE}, use parallel \code{foreach} to fit each fold. The default is \code{FALSE}.}
  \item{lambda}{The values of the regularization parameter lambda optimized over.}
   \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
 }

\details{
The function runs the \code{cv.glmnet} function of the \code{glmnet} package.
}



\value{
\item{optimal}{The value of lambda that gives the minimum cross-validated deviance.}
\item{results}{The data frame with the cross-validated deviance for each lambda value.}
}

\references{
Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, https://www.jstatsoft.org/v39/i05/
}

\examples{
data("dataDIVAT2")


formula<-Surv(times,failures) ~ age + hla + retransplant + ecd


tune.model <- tuneCOXridge(formula=formula, data=dataDIVAT2,
  cv=5, lambda=seq(0, 10, by=.1))

tune.model$optimal$lambda # the estimated lambda value

# The estimation of the training modelwith the corresponding lambda value
model <- LIB_COXridge(formula, data=dataDIVAT2,
  lambda=tune.model$optimal$lambda)

# The resulted predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Tunning}
