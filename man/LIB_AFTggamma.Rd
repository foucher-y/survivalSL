\name{LIB_AFTggamma}
\alias{LIB_AFTggamma}

\title{
Library of the Super Learner for an Accelerated Failure Time (AFT) Model with a Generalized Gamma Distribution
}

\description{
Fit an AFT parametric model with a generalized gamma distribution.
}

\usage{
LIB_AFTggamma(formula, data)
}

\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame whose columns correspond to the variables present in the formula.}
 }

\details{
The model is obtained by using the \code{dist="gengamma"} in the \code{flexsurvreg} package.
}

\value{
\item{formula}{The formula object used for model construction.}
\item{model}{The estimated model.}
\item{data}{The data frame used for learning.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}


\references{
Jackson, C. (2016). flexsurv: A Platform for Parametric Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33. doi:10.18637/jss.v070.i08
}

\examples{
data("dataDIVAT2")


# The estimation of the model from the first 200 lines

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_AFTggamma(formula=formula, data=dataDIVAT2[1:200,])

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}
