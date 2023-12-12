\name{lines.rocrisca}
\alias{lines.rocrisca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Lines to a ROC Plot
}
\description{
Used to add an additionnal ROC curve to ROC plot generated with \code{plot.rocrisca}.
}
\usage{
\method{lines}{rocrisca}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{rocrisca}, returned by the functions \code{roc.binary},
  \code{roc.net}, \code{roc.summary}, and \code{roc.time}.}
  \item{...}{Additional arguments affecting the plot line.}
}

\value{No return value for this S3 method.}

\examples{
# import and attach the data example
data(dataDIVAT3)

# A subgroup analysis to reduce the time needed for this exemple

dataDIVAT3 <- dataDIVAT3[1:400,]

# The standardized and weighted time-dependent ROC curve to evaluate the
# capacities of the recipient age for the prognosis of post kidney
# transplant mortality up to 2000 days by taking into account the
# donor age and the recipient gender.

# 0. Compute the raw sensitivity and specificity
roc1 <- roc.time(times="death.time", failures="death", variable="ageR",
confounders=~1, data=dataDIVAT3, pro.time=2000,
precision=seq(0.1,0.9, by=0.2))

# 1. Standardize the marker according to the covariates among the controls
lm1 <- lm(ageR ~ ageD + sexeR, data=dataDIVAT3[dataDIVAT3$death.time >= 2500,])
dataDIVAT3$ageR_std <- (dataDIVAT3$ageR - (lm1$coef[1] + lm1$coef[2] * dataDIVAT3$ageD +
 lm1$coef[3] * dataDIVAT3$sexeR)) / sd(lm1$residuals)

# 2. Compute the sensitivity and specificity from the proposed IPW estimators
roc2 <- roc.time(times="death.time", failures="death", variable="ageR_std",
confounders=~bs(ageD, df=3) + sexeR, data=dataDIVAT3, pro.time=2000,
precision=seq(0.1,0.9, by=0.2))

# The corresponding AUCs
c(roc1$auc, roc2$auc)

# The corresponding ROC graph

plot(roc1, type="b", col=1, pch=2, lty=2, xlab="1-specificity", ylab="sensibility")

lines(roc2, type="b", col=2, pch=1)

legend("bottomright", lty=1:2, lwd=1, pch=1:2, col=1:2,
 c(paste("Crude estimation, (AUC=", round(roc1$auc, 2), ")", sep=""),
 paste("Adjusted estimation, (AUC=", round(roc2$auc, 2), ")", sep="") ) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plotting}
