\name{LIB_COXridge}
\alias{LIB_COXridge}



\title{
Library of the Super Learner for Ridge Cox Regression
}

\description{
Fit a ridge Cox regression for a fixed value of the regularization parameter.
}



\usage{
LIB_COXridge(formula, data, penalty=NULL, lambda)
}


\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame whose columns correspond to the variables present in the formula.}
  \item{penalty}{A numerical vector with a length equals to the number of predictors. It allows the integration of covariates into the final model, i.e. with no selection: the value 0 to force the covariate in the model, 1 otherwise. If \code{NULL}, all covariates undergo the selection process.}
  \item{lambda}{The value of the regularization parameter lambda for penalizing the partial likelihood.}
 }


\details{
The ridge Cox regression is obtained by using the \code{glmnet} package.
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
model <- LIB_COXridge(formula=formula, data=dataDIVAT2[1:200,], lambda=.1)

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}

