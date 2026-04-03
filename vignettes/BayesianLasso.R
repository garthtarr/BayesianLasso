## ----include = FALSE---------------------------------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----helper-functions, include=FALSE-----------------------------------------------------------------
mcmc_stats = function(mBeta, vsigma2, vlambda2, time_val, inds_use)
{
  mBeta <- as.matrix(mBeta)
  N = length(inds_use)
  p <- ncol(mBeta)

  vESS <- numeric(p)
  for(j in 1:p) {
    vESS[j] <- effective_sample_size(mBeta[inds_use,j])
  }
  Ef = stats::median(vESS)/time_val
  ESS_sigma2  = effective_sample_size(as.vector(vsigma2[inds_use]))
  Ef_sigma2 = ESS_sigma2/time_val
  ESS_lambda2  = effective_sample_size(as.vector(vlambda2[inds_use]))
  Ef_lambda2 = ESS_lambda2/time_val

  stat_vec = c(
    100*stats::median(vESS)/N,
    Ef,
    100*ESS_sigma2/N,
    Ef_sigma2,
    100*ESS_lambda2/N,
    Ef_lambda2,
    time_val)

  name_vals = c("mix_beta", "eff_beta", "mix_sigma2", "eff_sigma2", "mix_lambda2", "eff_lambda2", "time")
  names(stat_vec) = name_vals

  return(stat_vec)
}

################################################################################

mcmc_diagnostics <- function(mBeta, vsigma2, vlambda2, beta_inds, mStat, doplots=TRUE)
{
  mBeta <- as.matrix(mBeta)

  if (doplots) {
    # Plot the acf for sigma2 and lambda2
    stats::acf(vsigma2)
    stats::acf(vlambda2)


    # Trace plots for sigma2 and lambda2
    graphics::plot(vsigma2, type = "l")
    graphics::plot(vlambda2, type = "l")
  }

  # Compute R-hat using posterior if available
  if (!requireNamespace("posterior", quietly = TRUE)) {
    message(
      "Package 'posterior' is required to compute R-hat diagnostics.\n",
      "Install it with install.packages('posterior') to enable rhat()."
    )
    rhat_sigma2  <- NA_real_
    rhat_lambda2 <- NA_real_
  } else {
    # posterior::rhat expects draws with chain dimension.
    # If we only have a single chain vector, rhat is not defined.
    # So we try rhat() and fall back to NA with a message if it errors.
    rhat_sigma2 <- tryCatch(
      as.numeric(posterior::rhat(vsigma2)),
      error = function(e) {
        message("posterior::rhat() failed for sigma2 (likely only 1 chain). Returning NA.")
        NA_real_
      }
    )
    rhat_lambda2 <- tryCatch(
      as.numeric(posterior::rhat(vlambda2)),
      error = function(e) {
        message("posterior::rhat() failed for lambda2 (likely only 1 chain). Returning NA.")
        NA_real_
      }
    )
  }

  # Choose beta indices to report densities for
  if (any(is.na(beta_inds))) {
    vbeta_hat <- colMeans(mBeta)
    beta_inds <- order(abs(vbeta_hat), decreasing = TRUE)[1:min(10, ncol(mBeta))]
  } else {
    beta_inds <- beta_inds[beta_inds >= 1 & beta_inds <= ncol(mBeta)]
    if (length(beta_inds) == 0) {
      beta_inds <- order(abs(colMeans(mBeta)), decreasing = TRUE)[1:min(10, ncol(mBeta))]
    }
  }

  ldens_beta <- lapply(beta_inds, function(j) {
    stats::density(mBeta[, j])
  })





  dens_sigma2  = stats::density(vsigma2)
  dens_lambda2 = stats::density(vlambda2)

  return(list(
    ldens_beta = ldens_beta,
    dens_sigma2 = dens_sigma2,
    dens_lambda2 = dens_lambda2,
    mStat = mStat,
    beta_inds = beta_inds,
    rhat_sigma2 = rhat_sigma2,
    rhat_lambda2 = rhat_lambda2
  ))
}


## ----ess-function, include=FALSE---------------------------------------------------------------------
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


## ----simulated-example, message=FALSE----------------------------------------------------------------

library(BayesianLasso)
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
      nsamples=Ns, lambda_init=1, sigma2_init=1, verbose=0, tune_lambda2 = TRUE,
      rao_blackwellization = FALSE)
  })[3]
  
  vtime_val_Hans[i] <- time_val
  
  # Initialize accumulators after first run
  if (is.null(results_Hans)) {
    results_Hans <- list(
      mBeta_Hans = res_Hans$mBeta,
      vsigma2_Hans = res_Hans$vsigma2,
      vlambda2_Hans = res_Hans$vlambda2
    )
  } else {
    results_Hans$mBeta_Hans    <- results_Hans$mBeta_Hans + res_Hans$mBeta
    results_Hans$vsigma2_Hans  <- results_Hans$vsigma2_Hans + res_Hans$vsigma2
    results_Hans$vlambda2_Hans <- results_Hans$vlambda2_Hans + res_Hans$vlambda2
  }
}

# Take averages
mBeta_Hans    <- (results_Hans$mBeta_Hans) / 5
vsigma2_Hans  <- (results_Hans$vsigma2_Hans) / 5
vlambda2_Hans <- (results_Hans$vlambda2_Hans) / 5
time_val_Hans = mean(vtime_val_Hans)




# ======================== PC sampler =======================================

# Run the modified PC sampler

vtime_val_PC = c()
results_PC <- NULL 

for(i in 1:5){
  time_val <- system.time({
    res_PC <- Modified_PC_Gibbs(
      X = X, y = y, a1=2, b1=1, u1=2, v1=1,
      nsamples=Ns, lambda_init=1, sigma2_init=1, verbose=0)
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




## ----convergence-diagnostics-------------------------------------------------------------------------

stats_Hans <- mcmc_stats(
  mBeta_Hans, vsigma2_Hans, vlambda2_Hans,
  time_val = time_val_Hans, inds_use = 200:Ns
)

print(stats_Hans)

mcmc_diagnostics(
  mBeta_Hans, vsigma2_Hans, vlambda2_Hans,
  beta_inds = 1:3, mStat = stats_Hans, doplots = TRUE
)


