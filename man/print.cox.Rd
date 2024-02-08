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
\code{\link{COXall}}, \code{\link{COXen}}, \code{\link{COXlasso}}, \code{\link{COXridge}}.
}

\examples{
data(dataDIVAT2)

model <- COXall(times="times", failures="failures", data=dataDIVAT2,
  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"))

print(model)
}

\keyword{Printing}
