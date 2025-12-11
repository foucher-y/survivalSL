\name{tuneRSF}
\alias{tuneRSF}

\title{
Tune a Survival Random Forest
}

\description{
This function finds the optimal nodesize, mtry, and ntree parameters for a survival random forest tree.
}

\usage{
tuneRSF(formula, data, nodesize=c(2, 4, 6, 10, 20, 30, 50, 100),
   mtry, ntree=500, seed=NULL)
}

\arguments{
 \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model with the same covariates as in the formula.}
    \item{nodesize}{The values of the node size optimized over.}
  \item{mtry}{The numbers of variables randomly sampled as candidates at each split optimized over.}
  \item{ntree}{The numbers of trees optimized over.}
   \item{seed}{A random seed to ensure reproducibility during the bootstrapping process. If \code{NULL}, a seed is randomly assigned.}
 }

\details{
The function runs the \code{rfsrc} function of the \code{randomForestSRC} package.
}

\value{
\item{optimal}{The value of lambda that gives the minimum mean cross-validated error.}
\item{results}{The data frame with the mean cross-validated errors for each lambda values.}
}

\references{
Ishwaran H. and Kogalur U.B. (2007). Random survival forests for R, Rnews, 7(2):25-31.
}

\examples{
data("dataDIVAT2")

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd

tune.model <- tuneRSF(formula, data=dataDIVAT2,
  nodesize=c(100, 250, 500), mtry=1, ntree=100)

tune.model$optimal # the estimated nodesize value

# The estimation of the training modelwith the corresponding lambda value
model <- LIB_RSF(formula, data=dataDIVAT2,
  nodesize=tune.model$optimal$nodesize, mtry=1, ntree=100)

# The resulted predicted survival of the first subject of the training sample
plot(y=model$predictions[1,], x=model$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))
}

\keyword{Tunning}
