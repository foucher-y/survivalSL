\name{tuneSNN}
\alias{tuneSNN}

\title{
Tune a 1-Layer Survival Neural Network
}

\description{
This function finds the optimal n.nodes, decay, batch.size, and epochs parameters for a survival neural network.
}

\usage{
tuneSNN(times, failures, group, cov.quanti,
cov.quali, data, cv, n.nodes, decay, batch.size, epochs)
}

\arguments{
  \item{times}{The name of the variable related the numeric vector with the follow-up times.}
  \item{failures}{The name of the variable related the numeric vector with the event indicators (0=right censored, 1=event).}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones. The default value is \code{NULL}: no specific exposure/treatment is considered. When a specific exposure/treatment is considered, it will be forced in the algorithm or related interactions will be tested when possible.}
  \item{cov.quanti}{The name(s) of the variable(s) related to the possible quantitative covariates. These variables must be numeric.}
  \item{cov.quali}{The name(s) of the variable(s) related to the possible qualitative covariates. These variables must be numeric with two levels: 0 and 1. A complete disjunctive form must be used for covariates with more levels.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time (\code{times}), the event (\code{failures}), the optional treatment/exposure (\code{group}) and the covariables included in the previous model (\code{cov.quanti} and \code{cov.quali}).}
    \item{cv}{The value of the number of folds. The default value is 10.}
  \item{n.nodes}{The number of hidden nodes optimized over.}
  \item{decay}{The value of the weight decay optimized over.}
  \item{batch.size}{The value of batch size}
  \item{epochs}{The value of epochs}
 }

\details{
The function runs the \code{deepsurv} function of the \code{survivalmodels} package. The metric used in the cross-validation is the C-index. Its use requires \code{Python}, but the \code{survivalSL} package can work without if this learner is not included in the library.
}

\value{
\item{optimal}{The value of \code{n.nodes}, \code{decay}, \code{batch.size}, and \code{epochs} that gives the maximum mean cross-validated C-index.}
\item{results}{The data frame with the mean cross-validated C-index according to \code{n.nodes}, \code{decay}, \code{batch.size}, and \code{epochs}.}
}

\references{
Katzman et al.  DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC Medical Research Methodology, 18(1), 24. 1018.
https://doi.org/10.1186/s12874-018-0482-1
}


\keyword{Tunning}
