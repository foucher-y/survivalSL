\name{print.flexsurv}
\alias{print.flexsurv}

\title{
S3 Method for Printing an 'flexsurv' Object
}
\description{
Print the model.
}

\usage{
\method{print}{flexsurv}(x,  ...)
}

\arguments{
  \item{x}{An object returned by the function \code{flexsurv}.}
  \item{...}{Further arguments passed.}
}

\value{No return value for this S3 method.}

\seealso{
  \code{\link{AFTgamma}}, \code{\link{AFTggamma}}, \code{\link{AFTllogis}}, \code{\link{AFTweibull}}, \code{\link{PHexponential}}, \code{\link{PHgompertz}}.
}

\examples{
data(dataDIVAT2)

model<-AFTgamma(times="times", failures="failures",  data=dataDIVAT2[1:100,],
        cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"))

print(model)
}

\keyword{Printing}
