
############################################################
# Reproducibility Script for BayesianLasso
#
# IMPORTANT NOTE:
#
# Running all R Markdown files in this script with the full
# MCMC settings (e.g., 11000 iterations per chain) may take
# several days on a standard desktop machine.
#
# The long runtime is due to:
#   - High-dimensional MCMC sampling
#   - Large benchmark datasets
#   - Multiple competing methods fitted for comparison
#
############################################################


results_dir <- "results"
stopifnot(dir.exists(results_dir))

datasets_ngtp <- c("diabetes2", "Kakadu2", "Crime")



for (d in datasets_ngtp) {
  print(d)
  rmarkdown::render(
    "benchmark_horseshoe_ngtp.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}

fun = median


#####################
# ---- Running models with Lasso prior and datasets n>p -----
#####################

for (d in datasets_ngtp) {
  print(d)
  rmarkdown::render(
    "benchmark_lasso_ngtp.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}




load(file.path(results_dir,"diabetes2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)


stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_diabetes2 = rbind(
  stat_hans,
  stat_4BG,
  
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"Crime_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)


stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_crime = rbind(
  stat_hans,
  stat_4BG,
  
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"Kakadu2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)

stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_Kakadu2 = rbind(
  stat_hans,
  stat_4BG,

  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


reff_crime = table_crime[,6]/table_crime[6,6]
reff_crime = round(reff_crime,2)

reff_diabetes2 = table_diabetes2[,6]/table_diabetes2[6,6]
reff_diabetes2 = round(reff_diabetes2,2)



reff_Kakadu2 = table_Kakadu2[,6]/table_Kakadu2[6,6]
reff_Kakadu2 = round(reff_Kakadu2,2)



table1  = cbind(reff_crime, reff_diabetes2, reff_Kakadu2)

table1


table2 = rbind(
  table_crime,
  table_diabetes2,
  table_Kakadu2
)


table2[table2[,1]>100,1] = 100
table2[table2[,3]>100,3] = 100
table2[table2[,5]>100,5] = 100

round(table2,1)




