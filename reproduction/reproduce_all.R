
############################################################
# Reproducibility Script for BayesianLasso

#
############################################################

# if (!requireNamespace("BayesianLasso", quietly = TRUE)) {
#   devtools::load_all(".")
# }
# reproduce_all.R



repro_dir <- dirname(normalizePath(sys.frame(1)$ofile))
setwd(repro_dir)

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

datasets_ngtp <- c("simulated", "diabetes2", "Kakadu2", "Crime")

time_val_total_ngtp <- system.time({
  for (d in datasets_ngtp) {
    dataset_name <- d
    source("Benchmarks_ngtp.R", local = FALSE)
  }
})[3]

print(time_val_total_ngtp)

datasets_pgtn <- c("cookie", "eyedata")



for (d in datasets_pgtn) {
  print(d)
  rmarkdown::render(
    "benchmarks_pgtn.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}


load(file.path(results_dir,"eyedata_results_pgtn.Rdata"))

print(colMeans(res_PC$mStat))
print(colMeans(res_hans$mStat))
print(colMeans(res_bayesreg$mStat))
print(time_val_total_pgtn)



