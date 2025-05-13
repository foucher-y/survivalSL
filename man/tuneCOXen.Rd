\name{tuneCOXen}
\alias{tuneCOXen}

\title{
Tune Elastic Net Cox Regression
}

\description{
This function finds the optimal lambda and alpha parameters for an elastic net Cox regression.
}

\usage{
tuneCOXen(formula, data, penalty=NULL,
 cv=10, parallel=FALSE, alpha=seq(.1,.9,.1), lambda=NULL, seed=NULL)
}

\arguments{
   \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model with the same covariates as in the formula.}
  \item{penalty}{A numerical vector that allows the covariates not to be penalized. We give the value 0 if we do not want the covariate to be penalized otherwise 1. If \code{NULL}, all covariates are penalized.}
  \item{cv}{The value of the number of folds. The default value is 10.}
  \item{parallel}{If \code{TRUE}, use parallel \code{foreach} to fit each fold. The default is \code{FALSE}.}
  \item{alpha}{The values of the regularization parameter alpha optimized over.}
  \item{lambda}{The values of the regularization parameter lambda optimized over.}
   \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
 }

\details{
The function runs the \code{cv.glmnet} function of the \code{glmnet} package.
}

\value{
\item{optimal}{The values of lambda and alpha that gives the minimum cross-validated deviance.}
\item{results}{The data frame with the cross-validated deviance for each lambda and alpha values.}
}

\references{
Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, https://www.jstatsoft.org/v39/i05/
}

\examples{
data("dataDIVAT2")

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
tune.model <- tuneCOXen(formula=formula, data=dataDIVAT2, cv=5,
  alpha=seq(.1, 1, by=.1), lambda=seq(.1, 1, by=.1))

tune.model$optimal$lambda # the estimated lambda value

# The estimation of the training modelwith the corresponding lambda value
model <- LIB_COXen(formula, data=dataDIVAT2,
  alpha=tune.model$optimal$alpha,
  lambda=tune.model$optimal$lambda)

# The resulted predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Tunning}
