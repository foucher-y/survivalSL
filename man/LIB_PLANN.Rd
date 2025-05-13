\name{LIB_PLANN}
\alias{LIB_PLANN}

\title{
Library of the Super Learner for Survival Neural Network Based on the PLANN Method
}

\description{
Fit a  neural network based on the partial logistic regression.
}

\usage{
LIB_PLANN(formula, data, inter, size, decay,
          maxit, MaxNWts, maxtime=NULL)
}

\arguments{
 \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame whose columns correspond to the variables present in the formula.}
  \item{inter}{The length of the intervals.}
  \item{size}{The number of units in the hidden layer.}
  \item{decay}{The parameter for weight decay.}
  \item{maxit}{The maximum number of iterations.}
  \item{MaxNWts}{The maximum allowable number of weights.}
  \item{maxtime}{A numeric value with the maximum prognostic time. If \code{NULL}, the maximum prognostic time is the maximum value of database times + 1.}
 }


\details{
This function is based is based on the \code{survivalPLANN} from the related package.
}

\value{
\item{formula}{The formula object used for model construction.}
\item{model}{The estimated model.}
\item{data}{The data frame used for learning.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}


\references{
Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.
}


\examples{
data("dataDIVAT2")

# The neural network based from the first 300 individuals of the data base

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_PLANN(formula, data=dataDIVAT2[1:300,],
  inter=0.5, size=32, decay=0.01, maxit=100, MaxNWts=10000, maxtime=NULL)

# The predicted survival of the first subject of the training sample

plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}
