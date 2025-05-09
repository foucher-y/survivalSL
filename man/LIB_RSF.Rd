\name{LIB_RSF}
\alias{LIB_RSF}

\title{
Library of the Super Learner for Survival Random Survival Forest
}

\description{
Fit survival random forest tree for given values of the regularization parameters.
}

\usage{
LIB_RSF(formula,
data, nodesize, mtry, ntree, seed=NULL)
}


\arguments{
 \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame whose columns correspond to the variables present in the formula.}
  \item{nodesize}{The value of the node size.}
  \item{mtry}{The number of variables randomly sampled as candidates at each split.}
  \item{ntree}{The number of trees.}
  \item{seed}{A random seed to ensure reproducibility during bootstrap sampling. If \code{NULL}, a seed is randomly assigned. }
 }

\details{
The survival random forest tree is obtained by using the \code{randomForestSRC} package.
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


# The estimation of the model

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_RSF(formula, data=dataDIVAT2, nodesize=10,
  mtry=2, ntree=100, seed=NULL)

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}
