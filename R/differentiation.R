differentiation <- function(x, fx) {
  n <- length(x)
  fdx0 <- vector(length = n)
  for (i in 2:n) {
    fdx0[i-1] <- (fx[i-1] - fx[i]) / (x[i-1] - x[i])   }
    fdx0[n] <- (fx[n] - fx[n - 1]) / (x[n] - x[n - 1])
  return(fdx0)  }
