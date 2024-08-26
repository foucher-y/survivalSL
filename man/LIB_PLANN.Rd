\name{LIB_PLANN}
\alias{LIB_PLANN}

\title{
Library of the Super Learner for Survival Neural Network Based on the PLANN Method
}

\description{
Fit a  neural network based on the partial logistic regression.
}

\usage{
LIB_PLANN(times, failures, group=NULL, cov.quanti=NULL, cov.quali=NULL,
data, inter, size, decay, maxit, MaxNWts)
}

\arguments{
  \item{times}{The name of the variable related the numeric vector with the follow-up times.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event).}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones. The default value is NULL: no specific exposure/treatment is considered. When a specific exposure/treatment is considered, it will be forced in the algorithm or related interactions will be tested when possible.}
  \item{cov.quanti}{The name(s) of the variable(s) related to the possible quantitative covariates. These variables must be numeric.}
  \item{cov.quali}{The name(s) of the variable(s) related to the possible qualitative covariates. These variables must be numeric with two levels: 0 and 1. A complete disjunctive form must be used for covariates with more levels.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time (\code{times}), the event (\code{failures}), the optional treatment/exposure (\code{group}) and the covariables included in the previous model (\code{cov.quanti} and \code{cov.quali}).}
  \item{inter}{The length of the intervals.}
  \item{size}{The number of units in the hidden layer.}
  \item{decay}{The parameter for weight decay.}
  \item{maxit}{The maximum number of iterations.}
  \item{MaxNWts}{The maximum allowable number of weights.}
 }


\details{
This function is based is based on the \code{survivalPLANN} from the  related package.
}

\value{
\item{model}{The estimated model.}
\item{group}{The name of the variable related to the exposure/treatment.}
\item{cov.quanti}{The name(s) of the variable(s) related to the possible quantitative covariates.}
\item{cov.quali}{The name(s) of the variable(s) related to the possible qualitative covariates.}
\item{data}{The data frame used for learning. The first column is entitled \code{times} and corresponds to the observed follow-up times. The second column is entitled \code{failures} and corresponds to the event indicators. The other columns correspond to the predictors.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\references{
Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.
}


\examples{
data(dataDIVAT2)

# The neural network based from the first 300 individuals of the data base

model <- LIB_PLANN(times="times", failures="failures", data=dataDIVAT2[1:300,],
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
  inter=0.5, size=32, decay=0.01, maxit=100, MaxNWts=10000)

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}
