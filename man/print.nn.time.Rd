\name{print.nn.time}
\alias{print.nn.time}

\title{
S3 Method for Printing an 'nn.time' Object
}
\description{
Print the survival neural network.
}

\usage{
\method{print}{nn.time}(x,  ...)
}

\arguments{
  \item{x}{An object returned by the function \code{nn.time}.}
  \item{...}{Further arguments passed.}
}

\seealso{
  \code{\link{nn.time}}.
}

\value{No return value for this S3 method.}

%\examples{
%\donttest{
%data(dataDIVAT2)

%model<-nn.time(times="times", failures="failures", data=dataDIVAT2[1:500,],
%  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), n.nodes=5,
%  decay=.01, batch.size=256L, epochs = 1L)
  
%print(model)}
%}

\keyword{Printing}