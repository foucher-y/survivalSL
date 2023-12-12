\name{differentiation}
\alias{differentiation}

\title{Numerical Differentiation with Finite Differences.}

\description{
This function allows to approximate the derivative of a function by the finite differences method.
}

\usage{differentiation(x, fx)}


\arguments{
  \item{x}{A numeric vector with the \eqn{x} values.}
  \item{fx}{A numeric vector with the \eqn{f(x)} values.}
 }
 
\details{
Numerical differentiation is a method of approximating the derivative of a function \eqn{f} at particular value \eqn{x}.
}

\examples{
x <- c(0.0, 0.2, 0.4)
fx <- c(0.00000, 0.74140, 1.3718)
differentiation(x, fx)
}
