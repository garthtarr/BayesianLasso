
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_samplers.h"
#include "lasso_distribution.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List Modified_Hans_Gibbs(const arma::vec& y,
                                     const arma::mat& X,
                                     const arma::vec& beta_init,
                                     double a1 = 1,
                                     double b1 = 1,
                                     double u1 = 1,
                                     double v1 = 1,
                                     int nsamples = 100,
                                     double lambda_init = 0.5,
                                     double sigma2_init = 1,
                                     int thin = 1,
                                     int verbose = 100,
                                     bool tune_lambda2 = true,
                                     bool rao_blackwellization = false)
{
  // Initialize dimensions of the problem
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  
  if (thin <= 0)
    Rcpp::stop("thin must be >= 1");

  const int maxiter = nsamples / thin;
  if (maxiter <= 0) stop("thin is larger than nsamples: no draws would be stored.");

  // Initialize storage of samples for MCMC
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  if (beta_init.size() != p) stop("beta_init must have length ncol(X).");

  // Declare variables for Rao-Blackwellization
  arma::mat mA;
  arma::mat mB;
  arma::mat mC;

  // Initialize variables for Rao-Blackwellization
  if (rao_blackwellization) {
    mA = zeros(maxiter,p);
    mB = zeros(maxiter,p);
    mC = zeros(maxiter,p);
  }

  // Initialize problem constants
  const arma::vec vone_n = arma::ones(n);
  const arma::vec vone_p = arma::ones(p);
  arma::vec dgXTX  = (X%X).t() * vone_n; // Diagonal values of XTX
  arma::mat XTX;
  if (n>p) {
    XTX = X.t() * X;
  }
  arma::vec XTy = X.t() * y;
  double yTy = dot(y, y);

  // Set the current values of the parameters
  arma::vec vd = arma::ones(p);
  arma::vec vnu = arma::ones(p);

  // Assign initial values
  arma::vec vbeta = beta_init;
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda*lambda;

  // Constant values
  const double a_til = a1 + 0.5*(n + p);
  const double u_til = u1 + 0.5*p;
  
  int keep_idx = 0;

  for (int i = 0; i < maxiter; ++i)
  {
    arma::vec va_vals = dgXTX/sigma2;
    arma::vec vb_vals = zeros(p);
    arma::vec vc_vals = vone_p*sqrt(lambda2/sigma2);

    arma::vec XTy_hat;
    arma::vec vy_hat;
    arma::vec vy_hat_mj;
    arma::vec vx_j;

    arma::vec vu = randu(p);
    if (n > p)
    {
      XTy_hat = XTX * vbeta;
      for (arma::uword j=0; j < p; ++j) {
        vx_j = XTX.col(j);
        XTy_hat = XTy_hat - vx_j*vbeta[j]; // This might not be exactly right
        vb_vals[j] = (XTy[j] - XTy_hat[j])/sigma2;
        vbeta[j] =  qlasso_fast_c_v2(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        XTy_hat = XTy_hat + vx_j*vbeta[j];
      }
    }
    else
    {
      vy_hat = X*vbeta;
      for (arma::uword j=0; j < p; ++j) {
        vx_j = X.col(j);
        vy_hat_mj = vy_hat - vx_j*vbeta[j];
        vb_vals[j] = (XTy[j] -  dot(vx_j, vy_hat_mj))/sigma2;
        vbeta[j] = qlasso_fast_c_v2(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        vy_hat = vy_hat_mj + vx_j*vbeta[j];
      }
    }

    // Store sufficient statistics for Rao-Blackwellization
    if (rao_blackwellization) {
      if (i % thin == 0) {
      mA.row(keep_idx) = va_vals.as_row();
      mB.row(keep_idx) = vb_vals.as_row();
      mC.row(keep_idx) = vc_vals.as_row();
      
      keep_idx++;
      }
    }

    ////////////////////////////////////////////////////////////////////////////

    // Slice from lambda2|rest
    double sum_abs_vbeta = sum(abs(vbeta));
    double RSS;
    if (n>p) {
      RSS = yTy - 2*dot(vbeta, XTy) + dot(vbeta, XTX*vbeta); // O(p^2)
    } else {
      RSS = sum(pow(y - vy_hat,2.0)); // O(n)
    }

    // Slice sample from sigma2|rest
    lambda = sqrt(lambda2);
    double a_val = (a_til-1);
    double b_val = b1 + 0.5*RSS;
    double c_val = lambda*sum_abs_vbeta;

    // Slice sampler for tau=1/sigma2 and then invert.
    double tau = 1/sigma2;
    tau = slice_sampler_precision(tau, a_val, b_val, c_val);
    sigma2 = 1/tau;
    sigma = sqrt(sigma2);

    ////////////////////////////////////////////////////////////////////////////

    // Slice sampler for lambda2
    if (tune_lambda2)
    {
      a_val = u_til - 1;
      b_val = v1;
      c_val = sum_abs_vbeta/sigma;

      lambda2 = slice_sampler_precision(lambda2, a_val, b_val, c_val);
      lambda  = sqrt(lambda2);
    }
    ////////////////////////////////////////////////////////////////////////////

    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    if (i % thin == 0) {
      
      mBeta.row(keep_idx) = vbeta.as_row();
      vsigma2[keep_idx]  = sigma2;
      vlambda2[keep_idx] = lambda2;
      
      keep_idx++;
      
    }

    // Store MCMC samples
   //  mBeta.row(i) = vbeta.as_row();
   // vsigma2[i] = sigma2;
   //  vlambda2[i] = lambda2;
  }

  return List::create(_["mBeta"] = mBeta,
                      _["vsigma2"] = vsigma2,
                      _["vlambda2"] = vlambda2,
                      _["mA"] = mA,
                      _["mB"] = mB,
                      _["mC"] = mC);
}



////////////////////////////////////////////////////////////////////////////////
