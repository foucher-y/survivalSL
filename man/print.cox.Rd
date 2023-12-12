\name{print.cox}
\alias{print.cox}

\title{
S3 Method for Printing an 'cox' Object
}
\description{
Print the Cox-like model.
}

\usage{
\method{print}{cox}(x,  ...)
}

\arguments{
  \item{x}{An object returned by the function \code{cox}.}
  \item{...}{Further arguments passed.}
}

\value{No return value for this S3 method.}

\seealso{
\code{\link{cox.all}}, \code{\link{cox.en}}, \code{\link{cox.lasso}}, \code{\link{cox.ridge}}.
}

\examples{
data(dataDIVAT2)

model<-cox.all(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"))

print(model)
}

\keyword{Printing}