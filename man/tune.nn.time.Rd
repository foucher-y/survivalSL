\name{tune.nn.time}
\alias{tune.nn.time}

\title{
Tune a 1-Layer Survival Neural Network
}

\description{
This function finds the optimal n.nodes, decay, batch.size, and epochs parameters for a survival neural network.
}

\usage{
tune.nn.time(times, failures, group, cov.quanti,
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
The function runs the \code{deepsurv} function of the \code{survivalmodels} package.
}

\value{
\item{optimal}{The value of lambda that gives the minimum mean cross-validated error.}
\item{results}{The data frame with the mean cross-validated errors for each lambda values.}
}

\references{
Katzman et al.  DeepSurv: personalized treatment recommender system using a Cox proportional hazards deep neural network. BMC Medical Research Methodology, 18(1), 24. 1018.
https://doi.org/10.1186/s12874-018-0482-1
}

%\examples{
%data(dataDIVAT2)

%# The estimation of the hyperparameters
%\donttest{
%tune.model<-tune.nn.time(times="times", failures="failures", data=dataDIVAT2[1:500,],
%    cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
%    cv=3, n.nodes=c(1, 2), decay=0, batch.size=256L, epochs=1L)}

%# the estimated nodesize value
%\donttest{
% tune.model$optimal
% tune.model$results}

%# The estimation of the training modelwith the corresponding lambda value
%\donttest{
% model<-nn.time(times="times", failures="failures", data=dataDIVAT2,
%   cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
%   n.nodes=tune.model$optimal$n.nodes, decay=0, batch.size=256L, epochs=1L)}

%# The resulted predicted survival of the first subject of the training sample
%\donttest{
%plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
% ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))}
%}

\keyword{Tunning}
