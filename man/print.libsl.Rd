\name{print.libsl}
\alias{print.libsl}

\title{
S3 Method for Printing an 'libsl' Object
}
\description{
Print the model or algorithm.
}

\usage{
\method{print}{libsl}(x, ...)
}

\arguments{
  \item{x}{A 'libsl' object.}
  \item{...}{For future methods.}
}

\value{No return value for this S3 method.}



\examples{
data("dataDIVAT2")


formula<-Surv(times,failures) ~ age + hla + retransplant + ecd
model <- LIB_AFTgamma(formula,  data=dataDIVAT2[1:100,])

print(model)
}

\keyword{Printing}
