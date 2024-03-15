survivalSL: an R Package for Predicting Survival by a Super Learner
================

## Description

The R package ‘survivalSL’ contains a variety of functions to construct
a super learner in the presence of censored times-to-event and to
evaluate its prognostic capacities. Several learners are proposed:
proportional hazard (PH) regressions, penalized PH semi-parametric
models, accelerated failure times (AFT) models, neural networks, random
survival forests, etc.). We proposed also a variety of loss functions
for the estimation of the weights (concordance index, Brier score, area
under the time-dependent ROC curve, negative binomial log-likelihood,
etc.). S3 methods are included to evaluate the predictive capacities, as
well as predicting survival curves from new observations.

## Basic Usage

``` r
# Simulate a training and validation samples
n.valid <- 500 # sample size for validation
n.learn <- 200 # sample size for training
n <- n.valid + n.learn # overall sample size

max.time <- 50 # maximum follow-up time

mean.x <- 0; sd.x <- 1 # normal distribution of the quantitative predictors
proba.x <- .5 # proportion of the binary predictors

a <- 2; b <- .05 # Weibull baseline distribution of the PH model
beta <- c(log(1.8), log(1.8), log(1.3), 0, 0, 0) # regression coefficients

# simulation of  the training and validation samples
x1 <- rnorm(n, mean.x, sd.x)
x2 <- rbinom(n, 1, proba.x)
x3 <- rbinom(n, 1, proba.x)
x4 <- rnorm(n, mean.x, sd.x)
x5 <- rbinom(n, 1, proba.x)
x6 <- rbinom(n, 1, proba.x)
x <- cbind(x1, x2, x3, x4, x5, x6) # matrix of the potential predictors
  
times <- 1/b*((-exp(-1*(x %*% beta))*(log(1-runif(n, 0, 1))))**(1/a)) # time to event
censoring <- runif(n, min=0, max=max.time)

status <- ifelse(times <= censoring, 1, 0) # event status
obs.times <- ifelse(times <= censoring, times, censoring) # follow-up times

data <- cbind(obs.times, status, as.data.frame(x))
  
data.simul <- list(data[1:n.valid,], data[(n.valid+1):n,])

# model estimation with default parameters and three learners
slres <- survivalSL(
  methods=c("LIB_COXen", "LIB_AFTgamma", "LIB_PHexponential"),
  metric="ci",  data=data.simul[[1]],  times="obs.times",
  failures="status", cov.quanti=c("x1","x4"),
  cov.quali=c("x2","x3","x5","x6"), progress = FALSE)
#> Warning in min(group %in% colnames(data)): no non-missing arguments to min;
#> returning Inf

# prognostic capacities from training sample
summary(slres, digits=3) 
#>      ci   auc    bs   ibs  ribs   bll  ibll ribll
#> 1 0.671 0.723 0.197 0.094 0.092 0.581 0.303 0.306

# prognostic capacities from validation sample
summary(slres, newdata=data.simul[[2]], digits=3) 
#>      ci  auc    bs   ibs ribs   bll  ibll ribll
#> 1 0.677 0.75 0.187 0.092 0.09 0.557 0.297 0.299
```

## Installation

To install the latest release from CRAN:

``` r
install.packages("survivalSL")
```

To install the development version from GitHub:

``` r
remotes::install_github("foucher-y/survivalSL")
```

## Reporting bugs

You can report any issues at this
[link](https://github.com/foucher-y/survivalSL/issues).
