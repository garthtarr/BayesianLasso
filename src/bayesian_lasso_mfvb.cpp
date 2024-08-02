//' @title bayesian lasso mfvb
//' @name bayesian lasso mfvb
//' @description TBA
//' @usage TBA
//' @details TBA
//' @returns TBA
//' @references Where did this come from?
//' @examples
//' rnorm(10)
//'

#include <cmath>

#include <RcppArmadillo.h>
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

typedef Map<MatrixXd> MapMat;
typedef Map<VectorXd> MapVec;

// bayesian_lasso_mfvb_c
// MFVB for Bayesian lasso regression
// Assumes that lambda is fixed
// Assumes that n>p

// [[Rcpp::export]]
List bayesian_lasso_mfvb_c(arma::mat mX, arma::vec vy, double lambda, double sigma2_hat, 
                           double a, double b, int maxiter, double tol, bool verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  // Initialisation
  arma::vec vmu_til = zeros(p); // vector of zeros
  arma::mat mSigma_til = eye(p, p); // Identity matrix
  arma::vec vsigma2_til;
   
  arma::vec va_til = ones(p); // expected value of auxiliary variables under q
  
  arma::mat mD;
  arma::mat mQ;
  arma::mat mQ_inv;
  
  double lambda2 = lambda*lambda;
  double a_til = a + 0.5*(n + p);
  double b_til = sigma2_hat*a_til;
  double sigma2inv_til = a_til/b_til;
  
  // Used to monitor convergence
  arma::vec vtheta = vmu_til;
  arma::vec vtheta_old;
    
  // Main MFVB loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Store old values
    vtheta_old = vmu_til;
    
    // Update q(beta)
    mD = diagmat(lambda2*va_til);
    mQ_inv = XTX + mD;
    mQ = inv(mQ_inv);
    vmu_til = mQ*XTy;
    mSigma_til = mQ/sigma2inv_til;
    vsigma2_til = diagvec(mSigma_til);

    arma::vec vmu2_til = pow(vmu_til,2.0);
    double RSS = sum(pow(vy-mX*vmu_til,2.0));
    double sum_Ea_vmu2 = sum(va_til%vmu2_til); 
        
    // Update q(sigma2)
    b_til = b + 0.5*(RSS + lambda2*sum_Ea_vmu2 + p/sigma2inv_til);
    sigma2inv_til = a_til/b_til;
   
    // Update q(va) 
    va_til = sqrt(1/(sigma2inv_til*lambda2*(vmu2_til + vsigma2_til)));
    
    // Check for convergence
    vtheta = vmu_til;
    double delta = norm(vtheta - vtheta_old);
    
//if (remainder(iter,100)==0) {
      if (verbose) {
        Rcout << "iter: " << i << " err:" << delta << "\n";
      }
      
//    }
    
    if (delta < tol) {
      break;
    }
  }
  
  return List::create(_["vmu_til"] = vmu_til, _["mSigma_til"] = mSigma_til, _["a_til"] = a_til, _["b_til"] = b_til);
}

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List bayesian_lasso_mfvb_tune_c(arma::mat mX, arma::vec vy, double lambda, double sigma2_hat, 
                           double a, double b, double u, double v,
                           int maxiter, double tol, bool verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  // Initialisation
  arma::vec vmu_til = zeros(p); // vector of zeros
  arma::mat mSigma_til = eye(p, p); // Identity matrix
  arma::vec vsigma2_til;
  
  arma::vec va_til = ones(p); // expected value of auxiliary variables under q
  
  arma::mat mD;
  arma::mat mQ;
  arma::mat mQ_inv;
  
  double a_til = a + 0.5*(n + p);
  double b_til = sigma2_hat*a_til;
  double sigma2inv_til = a_til/b_til;
  
  double u_til = u + 0.5*p;
  double v_til = u_til/(lambda*lambda);
  double lambda2_til = u_til/v_til;
  
  // Used to monitor convergence
  arma::vec vtheta = vmu_til;
  arma::vec vtheta_old;
  
  // Main MFVB loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Store old values
    vtheta_old = vmu_til;
    
    // Update q(beta)
    mD = diagmat(lambda2_til*va_til);
    mQ_inv = XTX + mD;
    mQ = inv(mQ_inv);
    vmu_til = mQ*XTy;
    mSigma_til = mQ/sigma2inv_til;
    vsigma2_til = diagvec(mSigma_til);
    
    arma::vec vmu2_til = pow(vmu_til,2.0);
    double RSS = sum(pow(vy-mX*vmu_til,2.0));
    double sum_Ea_vmu2 = sum(va_til%vmu2_til);
    
    // Update q(sigma2)
    b_til = b + 0.5*(RSS + lambda2_til*sum_Ea_vmu2 + p/sigma2inv_til);
    sigma2inv_til = a_til/b_til;
    
    // Update q(lambda2)
    v_til = v + 0.5*sigma2inv_til*(sum_Ea_vmu2 + sum(va_til%vsigma2_til));
    lambda2_til = u_til/v_til;
    
    // Update q(va) 
    va_til = sqrt(1/(sigma2inv_til*lambda2_til*(vmu2_til + vsigma2_til)));
    
    // Check for convergence
    vtheta = vmu_til;
    double delta = norm(vtheta - vtheta_old);
    
    //if (remainder(iter,100)==0) {
    if (verbose) {
      Rcout << "iter: " << i << " err:" << delta << "\n";
    }
    //    }
    
    if (delta < tol) {
      break;
    }
  }
  
  return List::create(_["vmu_til"] = vmu_til, 
                      _["mSigma_til"] = mSigma_til, 
                      _["a_til"] = a_til, 
                      _["b_til"] = b_til, 
                      _["u_til"] = u_til, 
                      _["v_til"] = v_til,
                      _["sigma2inv_til"] = sigma2inv_til, 
                      _["lambda2_til"] = lambda2_til);
}

////////////////////////////////////////////////////////////////////////////////

// looping through each column and element wise multiplication
// [[Rcpp::export]]
arma::mat matTimesVec(arma::mat mat, arma::vec v) {
  for(int i; i < mat.n_cols; i++){
    mat.col(i)  %=  v;
  }
  return mat;
}

// form a diagonal matrix with the vector and then use matrix multiplication
// [[Rcpp::export]]
arma::mat matTimesVec2(arma::mat mat, arma::vec v) {
  return arma::diagmat(v) * mat;
}

// use the functionality described at http://arma.sourceforge.net/docs.html#each_colrow 
// to "Apply a vector operation to each column or row of a matrix "
// [[Rcpp::export]]
arma::mat matTimesVec3(arma::mat mat, arma::vec v) {
  mat.each_col() %= v;
  return mat; 
}

////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List bayesian_lasso_mfvb_tune_pgtn_c(
    arma::mat mX, arma::vec vy, double lambda, double sigma2_hat, double a, double b, 
    double u, double v, int maxiter, double tol, bool verbose) {
  
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  //arma::mat XTX = mX.t()*mX;
  //arma::vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  // Initialisation
  arma::vec vmu_til = zeros(p); // vector of zeros
  arma::vec vmu2_til;
  arma::vec vsigma2_til;
  
  arma::mat mI_n = eye(n, n); // Identity matrix
  arma::vec va_til = ones(p); // expected value of auxiliary variables under q
  arma::vec vd_inv;
  arma::vec vv;
  arma::mat mQ;
  arma::mat XD;
  arma::mat mQ_inv;
  arma::mat QXD;

  double a_til = a + 0.5*(n + p);
  double b_til = sigma2_hat*a_til;
  double sigma2inv_til = a_til/b_til;
  
  double u_til = u + 0.5*p;
  double v_til = u_til/(lambda*lambda);
  double lambda2_til = u_til/v_til;
  
  // Used to monitor convergence
  arma::vec vtheta = vmu_til;
  arma::vec vtheta_old;
  
  // Main MFVB loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Store old values
    vtheta_old = vmu_til;
    
    /////////////////////////////////////
    // via Woodbury identity
    
    vd_inv = 1/(lambda2_til*va_til);
    
    //mQ = inv_sympd(mX*matTimesVec3(mX.t(), vd_inv) + mI_n, inv_opts::allow_approx);
    //vmu_til = vd_inv%(mX.t()*(mQ*vy));
    //vv = sum(mX % (mQ*mX), 0).t();
    //vsigma2_til = (vd_inv - pow(vd_inv,2.0) % vv)/Eq_sigma2inv;
    
    XD = mX.each_row() % vd_inv.t();
    mQ_inv = XD*mX.t() + mI_n;
    mQ = inv_sympd(mQ_inv, inv_opts::allow_approx);
    //QXD = solve(mQ_inv, XD, solve_opts::fast);
    QXD = mQ*XD; 
    vmu_til = QXD.t()*vy;
    vv = (XD % QXD).t()*ones(n);
    vsigma2_til = (vd_inv - vv)/sigma2inv_til;
    
    ////////////////////////////////////
    
    vmu2_til = pow(vmu_til,2.0);
    
    double sum_Ea_vmu2 = sum(va_til%vmu2_til);
    
    // Update q(sigma2)
    b_til = b + 0.5*(sum(pow(vy-mX*vmu_til,2.0)) + lambda2_til*sum_Ea_vmu2 + p/sigma2inv_til);
    sigma2inv_til = a_til/b_til;
    
    // Update q(lambda2)
    v_til = v + 0.5*sigma2inv_til*(sum_Ea_vmu2 + sum(va_til % vsigma2_til));
    lambda2_til = u_til/v_til;
    
    // Update q(va) 
    va_til = sqrt(1/(sigma2inv_til*lambda2_til*(vmu2_til + vsigma2_til)));
    
    // Check for convergence
    vtheta = vmu_til;
    double delta = norm(vtheta - vtheta_old);
    
    //if (remainder(iter,100)==0) {
    if (verbose) {
      Rcout << "iter: " << i << " err:" << delta << "\n";
    }
    //    }
    
    if (delta < tol) {
      break;
    }
  }
  
  return List::create(_["vmu_til"] = vmu_til, 
                      _["vsigma_til"] = vsigma2_til, 
                      _["a_til"] = a_til, 
                      _["b_til"] = b_til, 
                      _["u_til"] = u_til, 
                      _["v_til"] = v_til,
                      _["lambda2_til"] = lambda2_til, 
                      _["sigma2inv_til"] = sigma2inv_til,
                      _["va_til"] = va_til);
}
