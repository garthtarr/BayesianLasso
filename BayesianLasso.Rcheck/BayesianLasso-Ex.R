pkgname <- "BayesianLasso"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "BayesianLasso-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BayesianLasso')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LassoDistribution")
### * LassoDistribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LassoDistribution
### Title: The Lasso Distribution
### Aliases: LassoDistribution zlasso dlasso plasso qlasso rlasso elasso
###   vlasso mlasso MillsRatio Modified_Hans_Gibbs Modified_PC_Gibbs

### ** Examples

a <- 2; b <- 1; c <- 3
x <- seq(-3, 3, length.out = 1000)
plot(x, dlasso(x, a, b, c, logarithm = FALSE), type = 'l')

r <- rlasso(1000, a, b, c)
hist(r, breaks = 50, probability = TRUE, col = "grey", border = "white")
lines(x, dlasso(x, a, b, c, logarithm = FALSE), col = "blue")

plasso(0, a, b, c)
qlasso(0.25, a, b, c)
elasso(a, b, c)
vlasso(a, b, c)
mlasso(a, b, c)
MillsRatio(2)




# The Modified_Hans_Gibbs() function uses the Lasso distribution to draw 
# samples from the full conditional distribution of the regression coefficients.

y <- 1:20
X <- matrix(c(1:20,12:31,7:26),20,3,byrow = TRUE)

a1 <- b1 <- u1 <- v1 <- 0.01
sigma2_init <- 1
lambda_init <- 0.1
beta_init <- rep(1, ncol(X))
nsamples <- 1000
verbose <- 100
tune_lambda2 <- TRUE
rao_blackwellization <- FALSE

Output_Hans <- Modified_Hans_Gibbs(
                X, y, beta_init, a1, b1, u1, v1,
                nsamples, lambda_init, sigma2_init, 
                verbose, tune_lambda2, rao_blackwellization
)

colMeans(Output_Hans$mBeta)
mean(Output_Hans$vlambda2)


Output_PC <- Modified_PC_Gibbs(
               X, y, a1, b1, u1, v1, 
               nsamples, lambda_init, sigma2_init, verbose)

colMeans(Output_PC$mBeta)
mean(Output_PC$vlambda2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LassoDistribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalize")
### * normalize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalize
### Title: Normalize Response and Covariates
### Aliases: normalize

### ** Examples

set.seed(1)
X <- matrix(rnorm(100 * 10), 100, 10)
beta <- c(2, -3, rep(0, 8))
y <- as.vector(X %*% beta + rnorm(100))
norm_result <- normalize(y, X)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
