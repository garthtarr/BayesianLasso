
#' @export
benchmark_blasso_bayesreg <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{
  if (!requireNamespace("bayesreg", quietly = TRUE)) {
    stop(
      "The 'bayesreg' package is required for this benchmark.\n",
      "Install it with install.packages('bayesreg').",
      call. = FALSE
    )
  }

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }


  # Run the Gibbs sampler trials times
  dat <- data.frame(vy = vy, mX)

  mStat = c()
  for (i in seq_len(trials)) {
    # Initial fit

    time_val = system.time({
      res_mcmc <- bayesreg::bayesreg(
        vy~ .,
        data=dat,
        model = "normal",
        prior = "lasso",
        n.samples = (nsamples-nburn),
        burnin = nburn,
        thin = 1,
        t.dof = 25,
        n.cores = 1
      )
    })[3]
    # print(time_val)

    inds_use = 1:(nsamples - nburn)

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(t(res_mcmc$beta), res_mcmc$sigma2, res_mcmc$tau2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  colname_vals = c("eff_beta", "mix_beta", "eff_sigma2", "mix_sigma2", "eff_lambda2", "mix_lambda2", "time")
  colnames(mStat) <- colname_vals

  return(mcmc_diagnostics(t(res_mcmc$beta), res_mcmc$sigma2, res_mcmc$tau2, beta_inds, mStat, doplots = FALSE))
}
