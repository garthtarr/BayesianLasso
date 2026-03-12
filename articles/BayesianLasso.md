# BayesianLasso

``` r
library(BayesianLasso)

helpers_path <- system.file("helpers", package = "BayesianLasso")

helper_files <- list.files(helpers_path, pattern = "\\.R$", full.names = TRUE)

invisible(lapply(helper_files, source))
```

## Introduction

The package provides efficient sampling algorithms for Bayesian Lasso
regression, modified Hans and PC samplers. The modified Hans sampler is
based on a newly defined Lasso distribution.

This vignette introduces the main functionality of the package and
demonstrates how to apply the samplers to real and simulated data.

``` r
effective_sample_size <- function(samples) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    message(
      "Package 'posterior' is required to compute effective sample size (ESS).\n",
      "Install it with install.packages('posterior') to enable ESS diagnostics."
    )
    return(NA_real_)
  }
  
  as.numeric(posterior::ess_bulk(as.numeric(samples)))
}
```

## Dataset selection

Select a dataset name for dataset_name.

``` r

datasets_ngtp <- c("simulated", "diabetes2", "Kakadu2", "Crime")

datasets_pgtn <- c("cookie", "eyedata")

dataset_name <- datasets_ngtp[1]
```

## Simulated Example

``` r
if(dataset_name == datasets_ngtp[1]){

  # Simulate data
  set.seed(123)
  Ns <- 2000
  ns <- 100
  ps <- 10
  X <- matrix(rnorm(ns * ps), nrow = ns)
  beta <- c(rep(2, 3), rep(0, ps - 3))
  y <- X %*% beta + rnorm(ns)
  
  vtime_val_Hans = c()
  results_Hans <- NULL 

  # Run the modified Hans sampler
  for(i in 1:5){
    time_val <- system.time({
      res_Hans <- Modified_Hans_Gibbs(
      X = X, y = y, beta_init= rep(1,10), a1=2, b1=1, u1=2, v1=1,
      nsamples=Ns, lambda_init=1, sigma2_init=1, verbose=100, tune_lambda2 = TRUE,
      rao_blackwellization = FALSE)
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
  

  vESS <- c()
  for(j in 1:ps) {
    vESS[j] <- effective_sample_size(results_Hans$mBeta[200:Ns,j])
  }
  Ef_Hans = median(vESS)/time_val_Hans
  
  ESS_sigma2_Hans  = effective_sample_size(as.vector(results_Hans$vsigma2[200:Ns]))
  Ef_sigma2_Hans = ESS_sigma2_Hans/time_val_Hans
  
  ESS_lambda2_Hans  = effective_sample_size(as.vector(results_Hans$vlambda2[200:Ns]))
  Ef_lambda2_Hans = ESS_lambda2_Hans/time_val_Hans
  
  stat_vec_Hans = c(
    100*median(vESS)/Ns,
    Ef_Hans,
    100*ESS_sigma2_Hans/Ns,
    Ef_sigma2_Hans,
    100*ESS_lambda2_Hans/Ns,
    Ef_lambda2_Hans,
    time_val_Hans)
  
  name_vals = c("mix_beta", "eff_beta", "mix_sigma2", "eff_sigma2", "mix_lambda2", "eff_lambda2", "time")
  names(stat_vec_Hans) = name_vals 

  
  # ======================== PC sampler =======================================
  
  # Run the modified PC sampler
  
  vtime_val_PC = c()
  results_PC <- NULL 
  
 for(i in 1:5){
   time_val <- system.time({
     res_PC <- Modified_PC_Gibbs(
      X = X, y = y, a1=2, b1=1, u1=2, v1=1,
      nsamples=Ns, lambda_init=1, sigma2_init=1, verbose=100)
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
  mBeta = (results_PC$mBeta)/5
  vsigma2 = (results_PC$vsigma2)/5
  vlambda2 = (results_PC$vlambda2)/5
  time_val_PC = mean(vtime_val_PC)
  
  # colMeans(results_PC$mBeta)
  
  vESS <- c()
  for(j in 1:ps) {
    vESS[j] <- effective_sample_size(results_PC$mBeta[200:Ns,j])
  }
  Ef_PC = median(vESS)/time_val_PC
  
  ESS_sigma2_PC  = effective_sample_size(as.vector(results_PC$vsigma2[200:Ns]))
  Ef_sigma2_PC = ESS_sigma2_PC/time_val_PC
  
  ESS_lambda2_PC  = effective_sample_size(as.vector(results_PC$vlambda2[200:Ns]))
  Ef_lambda2_PC = ESS_lambda2_PC/time_val_PC
  
  stat_vec_PC = c(
    100*median(vESS)/Ns,
    Ef_PC,
    100*ESS_sigma2_PC/Ns,
    Ef_sigma2_PC,
    100*ESS_lambda2_PC/Ns,
    Ef_lambda2_PC,
    time_val_PC)
  
  name_vals = c("mix_beta", "eff_beta", "mix_sigma2", "eff_sigma2", "mix_lambda2", "eff_lambda2", "time")
  names(stat_vec_PC) = name_vals 
}
#> iter: 0 lambda2: 0.858958 sigma2: 1.00941
#> iter: 100 lambda2: 1.4063 sigma2: 0.958329
#> iter: 200 lambda2: 2.75352 sigma2: 0.864837
#> iter: 300 lambda2: 1.29269 sigma2: 1.04531
#> iter: 400 lambda2: 2.95129 sigma2: 0.953137
#> iter: 500 lambda2: 1.58173 sigma2: 0.90526
#> iter: 600 lambda2: 1.89708 sigma2: 0.962534
#> iter: 700 lambda2: 1.74399 sigma2: 1.25343
#> iter: 800 lambda2: 2.71632 sigma2: 0.992663
#> iter: 900 lambda2: 2.31175 sigma2: 0.915144
#> iter: 1000 lambda2: 1.87945 sigma2: 1.12389
#> iter: 1100 lambda2: 0.883053 sigma2: 1.23857
#> iter: 1200 lambda2: 2.21911 sigma2: 1.07196
#> iter: 1300 lambda2: 1.45428 sigma2: 1.0355
#> iter: 1400 lambda2: 2.10524 sigma2: 1.01913
#> iter: 1500 lambda2: 4.05418 sigma2: 1.07186
#> iter: 1600 lambda2: 2.45761 sigma2: 0.90785
#> iter: 1700 lambda2: 3.45027 sigma2: 1.34169
#> iter: 1800 lambda2: 3.95832 sigma2: 0.881618
#> iter: 1900 lambda2: 1.88884 sigma2: 1.24943
#> iter: 0 lambda2: 3.08123 sigma2: 1.43085
#> iter: 100 lambda2: 1.64498 sigma2: 1.07332
#> iter: 200 lambda2: 1.25335 sigma2: 0.980203
#> iter: 300 lambda2: 2.31797 sigma2: 1.28925
#> iter: 400 lambda2: 1.86186 sigma2: 0.977749
#> iter: 500 lambda2: 3.49106 sigma2: 1.06686
#> iter: 600 lambda2: 1.54199 sigma2: 1.09826
#> iter: 700 lambda2: 2.45581 sigma2: 1.2903
#> iter: 800 lambda2: 1.65776 sigma2: 1.18259
#> iter: 900 lambda2: 1.85546 sigma2: 1.03839
#> iter: 1000 lambda2: 2.24132 sigma2: 1.18037
#> iter: 1100 lambda2: 3.54543 sigma2: 1.2783
#> iter: 1200 lambda2: 1.37582 sigma2: 1.21287
#> iter: 1300 lambda2: 1.87356 sigma2: 1.12916
#> iter: 1400 lambda2: 1.33734 sigma2: 1.12834
#> iter: 1500 lambda2: 1.80435 sigma2: 1.16175
#> iter: 1600 lambda2: 1.29339 sigma2: 1.34696
#> iter: 1700 lambda2: 2.56371 sigma2: 0.954944
#> iter: 1800 lambda2: 2.45902 sigma2: 1.22607
#> iter: 1900 lambda2: 3.40984 sigma2: 1.19406
#> iter: 0 lambda2: 1.25219 sigma2: 0.960657
#> iter: 100 lambda2: 1.20597 sigma2: 1.01437
#> iter: 200 lambda2: 3.48412 sigma2: 1.20471
#> iter: 300 lambda2: 1.80046 sigma2: 1.19462
#> iter: 400 lambda2: 1.05257 sigma2: 1.0919
#> iter: 500 lambda2: 1.57488 sigma2: 0.96958
#> iter: 600 lambda2: 0.780784 sigma2: 0.961727
#> iter: 700 lambda2: 0.89589 sigma2: 0.859694
#> iter: 800 lambda2: 2.21023 sigma2: 0.895006
#> iter: 900 lambda2: 5.1711 sigma2: 1.04788
#> iter: 1000 lambda2: 4.29003 sigma2: 1.34692
#> iter: 1100 lambda2: 2.18293 sigma2: 0.981105
#> iter: 1200 lambda2: 4.07203 sigma2: 1.18968
#> iter: 1300 lambda2: 1.88576 sigma2: 1.06622
#> iter: 1400 lambda2: 2.09095 sigma2: 1.20465
#> iter: 1500 lambda2: 3.93494 sigma2: 1.24238
#> iter: 1600 lambda2: 1.4093 sigma2: 1.04512
#> iter: 1700 lambda2: 1.23756 sigma2: 0.891891
#> iter: 1800 lambda2: 1.25073 sigma2: 1.05552
#> iter: 1900 lambda2: 4.1184 sigma2: 1.01891
#> iter: 0 lambda2: 2.40241 sigma2: 0.969217
#> iter: 100 lambda2: 2.32971 sigma2: 0.983008
#> iter: 200 lambda2: 2.88361 sigma2: 1.09809
#> iter: 300 lambda2: 2.26611 sigma2: 1.03886
#> iter: 400 lambda2: 1.79334 sigma2: 0.994797
#> iter: 500 lambda2: 3.25743 sigma2: 1.10414
#> iter: 600 lambda2: 3.10534 sigma2: 0.995841
#> iter: 700 lambda2: 2.01161 sigma2: 0.925017
#> iter: 800 lambda2: 2.56433 sigma2: 1.45843
#> iter: 900 lambda2: 2.2818 sigma2: 1.05038
#> iter: 1000 lambda2: 3.98816 sigma2: 1.26292
#> iter: 1100 lambda2: 2.94681 sigma2: 1.03657
#> iter: 1200 lambda2: 2.28902 sigma2: 1.03245
#> iter: 1300 lambda2: 1.82286 sigma2: 1.10949
#> iter: 1400 lambda2: 1.46105 sigma2: 1.25766
#> iter: 1500 lambda2: 0.787921 sigma2: 0.911566
#> iter: 1600 lambda2: 2.7443 sigma2: 0.988332
#> iter: 1700 lambda2: 3.13911 sigma2: 1.32975
#> iter: 1800 lambda2: 2.80028 sigma2: 1.29629
#> iter: 1900 lambda2: 2.21642 sigma2: 0.93613
#> iter: 0 lambda2: 4.47944 sigma2: 1.01204
#> iter: 100 lambda2: 1.4653 sigma2: 1.02341
#> iter: 200 lambda2: 1.19634 sigma2: 1.10096
#> iter: 300 lambda2: 2.2262 sigma2: 0.999395
#> iter: 400 lambda2: 2.9787 sigma2: 0.922988
#> iter: 500 lambda2: 4.71374 sigma2: 0.85405
#> iter: 600 lambda2: 1.83211 sigma2: 0.96172
#> iter: 700 lambda2: 3.05674 sigma2: 1.16421
#> iter: 800 lambda2: 1.98883 sigma2: 1.01803
#> iter: 900 lambda2: 1.16042 sigma2: 1.18473
#> iter: 1000 lambda2: 0.936166 sigma2: 1.0562
#> iter: 1100 lambda2: 0.864847 sigma2: 1.18938
#> iter: 1200 lambda2: 3.82901 sigma2: 0.846567
#> iter: 1300 lambda2: 1.27783 sigma2: 0.878591
#> iter: 1400 lambda2: 1.11838 sigma2: 0.972412
#> iter: 1500 lambda2: 2.4752 sigma2: 0.87348
#> iter: 1600 lambda2: 1.46618 sigma2: 0.732426
#> iter: 1700 lambda2: 1.74941 sigma2: 1.0051
#> iter: 1800 lambda2: 2.86923 sigma2: 1.00961
#> iter: 1900 lambda2: 1.03685 sigma2: 0.887191
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
#> iter: 1000
#> iter: 1100
#> iter: 1200
#> iter: 1300
#> iter: 1400
#> iter: 1500
#> iter: 1600
#> iter: 1700
#> iter: 1800
#> iter: 1900
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
#> iter: 1000
#> iter: 1100
#> iter: 1200
#> iter: 1300
#> iter: 1400
#> iter: 1500
#> iter: 1600
#> iter: 1700
#> iter: 1800
#> iter: 1900
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
#> iter: 1000
#> iter: 1100
#> iter: 1200
#> iter: 1300
#> iter: 1400
#> iter: 1500
#> iter: 1600
#> iter: 1700
#> iter: 1800
#> iter: 1900
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
#> iter: 1000
#> iter: 1100
#> iter: 1200
#> iter: 1300
#> iter: 1400
#> iter: 1500
#> iter: 1600
#> iter: 1700
#> iter: 1800
#> iter: 1900
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
#> iter: 1000
#> iter: 1100
#> iter: 1200
#> iter: 1300
#> iter: 1400
#> iter: 1500
#> iter: 1600
#> iter: 1700
#> iter: 1800
#> iter: 1900
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
# Load libraries
library(BayesianLasso)


# Example: dataset_name <- "diabetes2"

if (dataset_name == datasets_ngtp[2]) {
  if (!requireNamespace("lars", quietly = TRUE)) {
    message("Package 'lars' is required for the diabetes example in this vignette.")
  }else{
    
    data("diabetes", package = "lars", envir = environment())
    
    y <- diabetes$y
    x <- diabetes$x
    inds <- seq_len(ncol(x))
    
    # Normalizing and scaling the dataset by function normalize()
    norm <- normalize(y, x, scale = TRUE)
    x <- norm$mX
    x <- model.matrix(~ .^2, data = data.frame(x = x))[ , -1]
    y <- norm$vy
  }
}

if (dataset_name == datasets_ngtp[3]) {
    if (!requireNamespace("Ecdat", quietly = TRUE)) {
        message("Package 'Ecdat' is required for the Kakadu example in this vignette.")
    }else{
      
      data("Kakadu", package = "Ecdat", envir = environment())
      
      # Get y vector and X matrix
      y <- as.vector(Kakadu$income)
      x <- Kakadu[, c(1:20, 22)]
      
      x <- model.matrix(~ .^2, data = x)[ , -1]
    }
}

# Make sure Crime.csv (or comData.Rdata) is downloaded manually before running this
if (dataset_name == datasets_ngtp[4]) 
{
  rdata_path <- system.file("extdata", "comData.Rdata", package = "BayesianLasso") 
  if (nzchar(rdata_path) && file.exists(rdata_path)) { load(rdata_path) 
  }else { stop("File 'comData.Rdata' not found in package extdata.") 
  }
  
  
  # ---- Remove rows with NA while keeping X and Y aligned ----
  datXY <- na.omit(cbind(as.data.frame(X), as.data.frame(Y)))
  X2 <- as.matrix(datXY[, colnames(X), drop = FALSE])
  Y2 <- as.matrix(datXY[, colnames(Y), drop = FALSE])
  
  # ---- Drop unwanted columns safely ----
  drop_cols <- c("ownHousQrange", "rentUpperQ")
  X2 <- X2[, !colnames(X2) %in% drop_cols, drop = FALSE]
  
  # ---- Define regression inputs ----
  x <- X2
  y <- as.vector(Y2[, "murders"])
  varnames <- colnames(x)
  inds <- seq_len(ncol(x))
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



if(dataset_name != datasets_ngtp[1]){
  
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
  mBeta = (results_PC$mBeta)/5
  vsigma2 = (results_PC$vsigma2)/5
  vlambda2 = (results_PC$vlambda2)/5
  time_val_PC = mean(vtime_val_PC)
  
  # Compute effective sample sizes  
  vESS <- c()
  p <- ncol(x)
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
        X = x, y = y, beta_init = as.vector(colMeans(mBeta)),
        a1, b1, u1, v1,
        nsamples,
        lambda_init = lambda_init,
        sigma2_init = sigma2_init,
        verbose = 1000, tune_lambda2 = TRUE, rao_blackwellization = FALSE
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
  p <- ncol(x)
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
  
  
  if(dataset_name != datasets_ngtp[3]){
    res_monomvn <- NULL
    
    if (requireNamespace("monomvn", quietly = TRUE)) {
      
      
      res_monomvn = benchmark_blasso_monomvn(y, x,
                                             lambda_init, sigma2_init, a1, b1, u1, v1, 
                                             nburn, nsamples, trials=5, 
                                             beta_inds=NA) 
      
      stat_vec_monomvn <- colMeans(res_monomvn$mStat)
      
      
    }else {
      message("Skipping monomvn benchmark: package 'monomvn' is not installed.")
    }
  }
  
  
  # ======================= Run bayeslm ============================
  
  res_bayeslm <- NULL
  
  if (requireNamespace("bayeslm", quietly = TRUE)) {
    
    res_bayeslm = benchmark_blasso_bayeslm(y, x,
                                           lambda_init, sigma2_init, a1, b1, u1, v1, 
                                           nburn, nsamples, trials=5, 
                                           beta_inds=NA) 
    
    stat_vec_bayeslm <- colMeans(res_bayeslm$mStat)
    
    
  }else {
    message("Skipping bayeslm benchmark: package 'bayeslm' is not installed.")
  }
  
  
  # ======================= Run rstan ============================
  
  res_rstan <- NULL
  
  if (requireNamespace("rstan", quietly = TRUE)) {
    
    res_rstan = benchmark_blasso_rstan(y, x,
                                       lambda_init, sigma2_init, a1, b1, u1, v1, 
                                       nburn, nsamples, trials=5, 
                                       beta_inds=NA)
    
    stat_vec_rstan <- colMeans(res_rstan$mStat)
    
    
  }else {
    message("Skipping rstan benchmark: package 'rstan' is not installed.")
  }
  
  
  # ======================= Run bayesreg ============================
  
  
  res_bayesreg <- NULL
  
  if (requireNamespace("bayesreg", quietly = TRUE)) {
    
    res_bayesreg = benchmark_blasso_bayesreg(y, x,
                                             lambda_init, sigma2_init, a1, b1, u1, v1, 
                                             nburn, nsamples, trials=5, 
                                             beta_inds=NA) 
    
    stat_vec_bayesreg <- colMeans(res_bayesreg$mStat)
    
    
  }else {
    message("Skipping bayesreg benchmark: package 'bayesreg' is not installed.")
  }
  
}
```

## Convergence diagnostics

Assessing convergence of the MCMC chains is an important step in
Bayesian analysis. The package provides helper functions such as
`mcmc_stats()` and `mcmc_diagnostics()` to summarize mixing behavior and
visualize the behavior of the chains. The function `mcmc_stats()`
computes numerical summaries of the MCMC output, while
`mcmc_diagnostics()` produces basic diagnostic plots such as trace and
autocorrelation plots for selected parameters. These tools allow users
to assess mixing and stability of the chains in practice. Posterior
draws can also be passed to packages such as `coda`, `posterior`, or
`bayesplot` for additional diagnostics.
