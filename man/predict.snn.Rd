\name{predict.snn}
\alias{predict.snn}

\title{
Prediction from a Survival Neural Network
}
\description{
Predict the survival of new observations based on a neural network estimated by using \code{SNN}.
}

\usage{
\method{predict}{snn}(object, ..., newdata, newtimes)
}

\arguments{
  \item{object}{An object returned by the function \code{SNN}.}
  \item{...}{Further arguments passed.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in \code{cov.quanti} and \code{cov.quali} included in the training sample. The default value is \code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the observed times in the training data frame.}
}


\details{
It requires \code{Python}, but the \code{survivalSL} package can work without if this learner is not included in the library.
}


\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\seealso{
  \code{\link{SNN}}.
}

\keyword{Predicting}
