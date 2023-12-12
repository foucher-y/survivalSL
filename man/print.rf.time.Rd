\name{print.rf.time}
\alias{print.rf.time}

\title{
S3 Method for Printing an 'rf.time' Object
}
\description{
Print the random survival forest.
}

\usage{
\method{print}{rf.time}(x,  ...)
}

\arguments{
  \item{x}{An object returned by the function \code{rf.time}.}
  \item{...}{Further arguments passed.}
}

\value{No return value for this S3 method.}

\seealso{
  \code{\link{rf.time}}.
}

\examples{
data(dataDIVAT2)

model<-rf.time(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
  nodesize=10, mtry=2, ntree=100)

print(model)
}

\keyword{Printing}