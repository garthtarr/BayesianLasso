# BayesianLasso

``` r
library(BayesianLasso)
```

## Introduction

The package provides efficient sampling algorithms for Bayesian Lasso
regression, modified Hans and PC samplers. The modified Hans sampler is
based on a newly defined Lasso distribution.

This vignette introduces the main functionality of the package and
demonstrates how to apply the samplers to real and simulated data.

## Simulated Example

``` r
# Simulate data
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), nrow = n)
beta <- c(rep(2, 3), rep(0, p - 3))
y <- X %*% beta + rnorm(n)

# Run the modified Hans sampler
output_Hans <- Modified_Hans_Gibbs(
               X = X, y = y, a1=2, b1=1, u1=2, v1=1, 
               nsamples=1000, beta_init= rep(1,10), lambda_init=1, sigma2_init=1, verbose=100)
#> iter: 0 lambda2: 4.91771 sigma2: 0.960406
#> iter: 100 lambda2: 3.58811 sigma2: 0.997586
#> iter: 200 lambda2: 0.933625 sigma2: 0.901476
#> iter: 300 lambda2: 0.794242 sigma2: 0.991402
#> iter: 400 lambda2: 3.81367 sigma2: 1.01224
#> iter: 500 lambda2: 2.58749 sigma2: 0.906112
#> iter: 600 lambda2: 1.17227 sigma2: 1.2331
#> iter: 700 lambda2: 3.38756 sigma2: 1.27917
#> iter: 800 lambda2: 0.676862 sigma2: 0.881326
#> iter: 900 lambda2: 2.87214 sigma2: 1.21127
colMeans(output_Hans$mBeta)
#>  [1]  2.05825310  1.96658242  1.85249926  0.12747215  0.06827941 -0.04384691
#>  [7]  0.07039560  0.10960198 -0.04258196  0.12511862

# Run the modified PC sampler
output_PC <- Modified_PC_Gibbs(
               X = X, y = y, a1=2, b1=1, u1=2, v1=1, 
               nsamples=1000, lambda_init=1, sigma2_init=1, verbose=100)
#> iter: 0
#> iter: 100
#> iter: 200
#> iter: 300
#> iter: 400
#> iter: 500
#> iter: 600
#> iter: 700
#> iter: 800
#> iter: 900
colMeans(output_PC$mBeta)
#>  [1]  2.05197868  1.96502229  1.84986538  0.12355534  0.06721269 -0.04147922
#>  [7]  0.06621934  0.11292249 -0.03955249  0.12324313
```

## Data Sources

This vignette demonstrates the usage of the `BayesianLasso` package on
several well-known datasets. To avoid redistribution issues and keep the
package lightweight, none of these datasets are included directly in the
package. Instead, we provide guidance on how to access them.

- **Diabetes**: Available in the package and loaded using
  `data(diabetes, package = "lars")`.

- **Kakadu**: Available in the package and loaded using
  `data(Kakadu, package = "Ecdat")`.

- **Communities and Crime**: This dataset is available from the [UCI
  Machine Learning
  Repository](https://archive.ics.uci.edu/dataset/183/communities+and+crime).  
  Due to licensing and redistribution considerations, the *Communities
  and Crime* dataset is not included in this package. Users can download
  the dataset manually by visiting the repository and following the
  preprocessing steps outlined in this vignette.

## Reproducing Table 1: Performance Comparison

This section demonstrates how to reproduce the performance comparison
table from the manuscript, comparing the mixing percentages, sampling
efficiencies, and runtimes for the Diabetes, Kakadu, and Communities and
Crime datasets.

``` r
effective_sample_size <- function(samples) 
{
  # Using the spectral method
  N <- length(samples)
  
  # Compute spectral density at zero using spectrum0.ar (uses AR smoothing)
  spectral_density_zero <- spectrum0.ar(samples)$spec
  
  # Compute effective sample size using variance ratio
  sample_variance <- var(samples)
  tau <- spectral_density_zero / sample_variance  # Integrated autocorrelation time
  ESS <- N / tau
  if (ESS > N) {
    ESS = N
  }
  
  return(ESS)
}
```

``` r
# Load libraries
library(BayesianLasso)
# This code requires bayesreg

if (requireNamespace("monomvn", quietly = TRUE)) {
  monomvn::monomvn(...)  # use the function safely
} else {
  message("monomvn package not installed; skipping example.")
}
if (requireNamespace("bayeslm", quietly = TRUE)) {
  bayeslm::bayeslm(...)  # use the function safely
} else {
  message("bayeslm package not installed; skipping example.")
}
if (requireNamespace("rstan", quietly = TRUE)) {
  rstan::rstan(...)  # use the function safely
} else {
  message("rstan package not installed; skipping example.")
}
if (requireNamespace("bayesreg", quietly = TRUE)) {
  bayesreg::bayesreg(...)  # use the function safely
} else {
  message("bayesreg package not installed; skipping example.")
}


# Example: dataset_name <- "diabetes2"

if (dataset_name=="diabetes2") 
{
  if (!requireNamespace("lars", quietly = TRUE)) {
    install.packages("lars")
  }
  data(diabetes, package = "lars")

  y = diabetes$y
  x = diabetes$x
  inds = 1:ncol(x)
  
  # Normalizing and scaling the dataset by function normalize()
  norm = normalize(y,x, scale = TRUE)
  x = norm$mX
  x <- model.matrix(~.^2, data=data.frame(x=x))[,-1]
  y <- norm$vy
  
}

if (dataset_name=="Kakadu2") 
{
  if (!requireNamespace("Ecdat", quietly = TRUE)) {
    install.packages("Ecdat")
  }
  dat <- data("Kakadu", package = "Ecdat")
  

  # Get y vector and X matrix
  y <- as.vector(dat$income)
  x <- dat[,c(2:21,23)]  
  
  x <- model.matrix(~.^2,data=x)[,-1]
}

# Make sure Crime.csv (or comData.Rdata) is downloaded manually before running this
if (dataset_name=="Crime") 
{
  rdata_path <- system.file("extdata", "comData.Rdata", package = "BayesianLasso")
  load(rdata_path)
  # Show only if file exists
  crime_path <- "Crime.csv"
  if (file.exists(crime_path)) {
    crime_data <- read.csv(crime_path)
  } else {
    message("Please download 'communities.data' from the UCI Machine Learning Repository and convert to Crime.csv.")
  }
  

  mX <- t(t(X) %>% na.omit())  
  mX <- mX[,-which(colnames(X)%in%c("ownHousQrange","rentUpperQ"))]
  varnames <- colnames(mX)
  vy <- Y[,"murders"]
  x <- mX
  y <- vy
  inds = 1:ncol(x)
}


# Set prior hyperparameter constants
a1 = 1.0E-2  # Prior shape for sigma2
b1 = 1.0E-2  # Prior scale for sigma2
u1 = 1.0E-2  # Prior shape for lambda2
v1 = 1.0E-2  # Prior scale for lambda2

# Initial values for lambda2 and sigma2
lambda2_init  = 10
lambda_init = sqrt(lambda2_init)
sigma2_init = 1

# Number of samples to run the MCMC
nburn = 1000
nsamples = 5000
inds_use = (nburn + 1):nsamples
N = length(inds_use)

# To store elapsed time and results of the PC sampler
vtime_val_PC = c()
results_PC <- NULL  # to be initialized after the first run

# Running the modified PC sampler 5 times and taking the average of the results across runs.
for (i in 1:5) {
  time_val <- system.time({
    res_PC <- Modified_PC_Gibbs(
      X = x, y = y, a1, b1, u1, v1,
      nsamples,
      lambda_init = lambda_init,
      sigma2_init = sigma2_init,
      verbose = 1000
    )
  })[3]

  vtime_val_PC[i] <- time_val

  # Initialize accumulators after first run
  if (is.null(results_PC)) {
    results_PC <- list(
      mBeta = res_PC$mBeta,
      vsigma2 = res_PC$vsigma2,
      vlambda2 = res_PC$vlambda2
    )
  } else {
    results_PC$mBeta    <- results_PC$mBeta + res_PC$mBeta
    results_PC$vsigma2  <- results_PC$vsigma2 + res_PC$vsigma2
    results_PC$vlambda2 <- results_PC$vlambda2 + res_PC$vlambda2
  }
}

# Take averages
mBeta = (res_PC$mBeta)/5
vsigma2 = (res_PC$vsigma2)/5
vlambda2 = (res_PC$vlambda2)/5
time_val_PC = mean(vtime_val_PC)

# Compute effective sample sizes  
  vESS <- c()
  for(j in 1:p) {
    vESS[j] <- effective_sample_size(mBeta[inds_use,j])
  }
  Ef_PC = median(vESS)/time_val_PC
  
  ESS_sigma2_PC  = effective_sample_size(as.vector(vsigma2[inds_use]))
  Ef_sigma2_PC = ESS_sigma2_PC/time_val_PC
  
  ESS_lambda2_PC  = effective_sample_size(as.vector(vlambda2[inds_use]))
  Ef_lambda2_PC = ESS_lambda2_PC/time_val_PC
  
  stat_vec_PC = c(
    100*median(vESS)/N,
    Ef_PC,
    100*ESS_sigma2_PC/N,
    Ef_sigma2_PC,
    100*ESS_lambda2_PC/N,
    Ef_lambda2_PC,
    time_val_PC)
  
name_vals = c("mix_beta", "eff_beta", "mix_sigma2", "eff_sigma2", "mix_lambda2", "eff_lambda2", "time")
names(stat_vec_PC) = name_vals 


# To store elapsed time and results of the Hans sampler
vtime_val_Hans = c()
results_Hans <- NULL  # to be initialized after the first run

# Running the modified Hans sampler 5 times and taking the average of the results across runs.
for (i in 1:5) {
  time_val <- system.time({
    res_Hans <- Modified_Hans_Gibbs(
      X = x, y = y, a1, b1, u1, v1,
      nsamples,
      beta_init = as.vector(colMeans(mBeta)),
      lambda_init = lambda_init,
      sigma2_init = sigma2_init,
      verbose = 1000
    )
  })[3]

  vtime_val_Hans[i] <- time_val

  # Initialize accumulators after first run
  if (is.null(results_Hans)) {
    results_Hans <- list(
      mBeta = res_Hans$mBeta,
      vsigma2 = res_Hans$vsigma2,
      vlambda2 = res_Hans$vlambda2
    )
  } else {
    results_Hans$mBeta    <- results_Hans$mBeta + res_Hans$mBeta
    results_Hans$vsigma2  <- results_Hans$vsigma2 + res_Hans$vsigma2
    results_Hans$vlambda2 <- results_Hans$vlambda2 + res_Hans$vlambda2
  }
}

# Take averages
mBeta = (results_Hans$mBeta)/5
vsigma2 = (results_Hans$vsigma2)/5
vlambda2 = (results_Hans$vlambda2)/5
time_val_Hans = mean(vtime_val_Hans)

# Compute effective sample sizes
  vESS <- c()
  for(j in 1:p) {
    vESS[j] <- effective_sample_size(mBeta[inds_use,j])
  }
  Ef_Hans = median(vESS)/time_val_Hans
  
  ESS_sigma2_Hans  = effective_sample_size(as.vector(vsigma2[inds_use]))
  Ef_sigma2_Hans = ESS_sigma2_Hans/time_val_Hans
  
  ESS_lambda2_Hans  = effective_sample_size(as.vector(vlambda2[inds_use]))
  Ef_lambda2_Hans = ESS_lambda2_Hans/time_val_Hans
  
  stat_vec_Hans = c(
    100*median(vESS)/N,
    Ef_Hans,
    100*ESS_sigma2_Hans/N,
    Ef_sigma2_Hans,
    100*ESS_lambda2_Hans/N,
    Ef_lambda2_Hans,
    time_val_Hans)
  
name_vals = c("mix_beta", "eff_beta", "mix_sigma2", "eff_sigma2", "mix_lambda2", "eff_lambda2", "time")
names(stat_vec_Hans) = name_vals 
```
