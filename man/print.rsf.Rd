\name{print.rsf}
\alias{print.rsf}

\title{
S3 Method for Printing an 'rsf' Object
}
\description{
Print the random survival forest.
}

\usage{
\method{print}{rsf}(x,  ...)
}

\arguments{
  \item{x}{An object returned by the function \code{RSF}.}
  \item{...}{Further arguments passed.}
}

\value{No return value for this S3 method.}

\seealso{
  \code{\link{RSF}}.
}

\examples{
data(dataDIVAT2)

model <- RSF(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"),
  nodesize=10, mtry=2, ntree=100)

print(model)
}

\keyword{Printing}
