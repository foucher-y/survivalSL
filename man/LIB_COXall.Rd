\name{LIB_COXall}
\alias{LIB_COXall}

\title{
Library of the Super Learner for Cox Regression
}

\description{
Fit a Cox regression for all covariates to be used in the super learner.
}

\usage{
LIB_COXall(formula, data)
}

\arguments{
\item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model in which to look for the variables related to the status of the follow-up time, the event,and the other covariables included in the formula.}
 }

\details{
The Cox regression is obtained by using the \code{survival} package.
}

\value{
\item{formula}{The formula object used for model construction.}
\item{model}{The estimated model.}
\item{data}{The data frame used for learning.}
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\references{
Terry M. Therneau (2021). A Package for Survival Analysis in R. R package version 3.2-13, https://CRAN.R-project.org/package=survival.
}

\examples{
data("dataDIVAT2")


# The estimation of the model from the first 200 lines

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
data<-dataDIVAT2[1:200,]
model <- LIB_COXall(formula=formula, data=data)

# The predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Super Learner Library}
