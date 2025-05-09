print.sltime <- function (x, digits=7, ...)
{
  cat("The contribution of the learners: ", "\n\n", sep = "")

  print(data.frame(leaners = x$methods, weights = round(x$weights$values, digits)))

  cat("\n", "Minimum of the ", x$cv, "-fold CV of the metric ",
      x$metric$metric, ":", round(x$metric$value, digits), ".", sep = "")
}
