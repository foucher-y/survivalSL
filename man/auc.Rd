\name{auc}
\alias{auc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Area Under ROC Curve from Sensitivities and Specificities.
}
\description{
This function computes the area under ROC curve by using the trapezoidal rule.
}
\usage{
auc(sens, spec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sens}{A numeric vector with the sensitivities}
  \item{spec}{A numeric vector with the specificities}
  }
  
\details{
This function computes the area under ROC curve using the trapezoidal rule from two vectors of sensitivities and specificities. The value of the area is directly returned.
}

\value{
The estimation of the area under ROC curve.
}

\examples{
se.temp <- c(0, 0.5, 0.5, 1)
sp.temp <- c(1, 0.5, 0.5, 0)
auc(se.temp, sp.temp)
}


\keyword{Metrics}
