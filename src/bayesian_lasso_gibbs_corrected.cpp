
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"

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

////////////////////////////////////////////////////////////////////////////////

// bayesian_lasso_gibbs_c
// Gibbs sampling for Bayesian lasso regression
// Assumes that lambda is fixed
// Assumes that n>p

// [[Rcpp::export]]
List bayesian_lasso_gibbs_c(arma::mat mX, arma::vec vy, double lambda, double sigma2_hat, 
                            double a, double b, double u, double v, int nsamples, 
                            int verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  double yTy = sum(vy%vy);
  //double yTy = vy.t()*vy;
  
  int maxiter = nsamples;
  
  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mV(maxiter,p);
  arma::vec va_til(maxiter);
  arma::vec vb_til(maxiter);
  
  // Initialisation
  arma::vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_hat;
  double lambda2 = lambda*lambda;
  double a_til = a + 0.5*(n + p);
  double b_til;
  double u_til = u + 0.5*p;
  double v_til;
  
  arma::mat vmu_til;
  arma::mat mSigma_til;
  arma::vec vbeta;
  arma::mat mQ;
  
  // Parameters of va|rest which is inverse Gaussian
  arma::vec vmu;
  arma::vec vlambda = ones(p);
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    mQ = inv(XTX + diagmat(lambda2*va));
    vmu_til = mQ*XTy;
    mSigma_til = sigma2*mQ;
    vbeta = mvnrnd(vmu_til, mSigma_til, 1);
    
    // Sample from sigma2|rest
    arma::vec XTX_b = XTX * vbeta;
    //Rcout << "iter: " << i << "\n";
    double bT_XTX_b = sum(vbeta%XTX_b);
    double bT_XTy = sum(vbeta%XTy);
    double RSS  =  yTy - 2*bT_XTy + bT_XTX_b;//sum(pow(vy-mX*vbeta,2.0));
    
    
    b_til = b + 0.5*(RSS + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(va%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
   
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rinvgaussian_c(vmu, vlambda);
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << "\n";
      }
    }
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = diagvec(mSigma_til).as_row();
    va_til[i] = a_til;
    vb_til[i] = b_til;
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["va_til"] = va_til, 
                      _["vb_til"] = vb_til);
}

////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List Modified_PC_Gibbs(arma::mat mX, arma::vec vy,double a1, double b1, double u1, double v1, 
                       int nsamples, double lambda_init, double sigma2_init, 
                       int verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  int maxiter = nsamples;
  
  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mV(maxiter,p);
  arma::vec va_til(maxiter);
  arma::vec vb_til(maxiter);
  arma::vec vu_til(maxiter);
  arma::vec vv_til(maxiter);
  
  // Initialisation
  arma::vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_init;
  double lambda2 = lambda_init*lambda_init;
  double a_til = a1 + 0.5*(n + p);
  double b_til;
  double u_til = u1 + 0.5*p;
  double v_til;
  
  arma::mat vmu_til;
  arma::mat mSigma_til;
  arma::vec vbeta;
  arma::mat mQ;
  
  // Parameters of va|rest which is inverse Gaussian
  arma::vec vmu;
  arma::vec vlambda = ones(p);
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    mQ = inv(XTX + diagmat(lambda2*va));
    vmu_til = mQ*XTy;
    mSigma_til = sigma2*mQ;
    vbeta = mvnrnd(vmu_til, mSigma_til, 1);
    
    // Sample from sigma2|rest
    b_til = b1 + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v1 + 0.5*sum(va%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rinvgaussian_c(vmu, vlambda);
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << "\n";
      }
    }
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = diagvec(mSigma_til).as_row();
    va_til[i] = a_til;
    vb_til[i] = b_til;
    vu_til[i] = u_til;
    vv_til[i] = v_til;
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["va_til"] = va_til, 
                      _["vb_til"] = vb_til,
                      _["vu_til"] = vu_til, 
                      _["vv_til"] = vv_til);
}


// looping through each column and element wise multiplication
// this is included in another file I think
// // [[Rcpp::export]]
// arma::mat matTimesVec(arma::mat mat, arma::vec v) {
//   for(int i; i < mat.n_cols; i++){
//     mat.col(i)  %=  v;
//   }
//   return mat;
// }

// form a diagonal matrix with the vector and then use matrix multiplication
// this is included in another file I think
// // [[Rcpp::export]]
// arma::mat matTimesVec2(arma::mat mat, arma::vec v) {
//   return arma::diagmat(v) * mat;
// }

// use the functionality described at http://arma.sourceforge.net/docs.html#each_colrow 
// to "Apply a vector operation to each column or row of a matrix "
// this is included in another file I think
// // [[Rcpp::export]]
// arma::mat matTimesVec3(arma::mat mat, arma::vec v) {
//   mat.each_col() %= v;
//   return mat; 
// }

arma::mat crossprod_cpp(arma::mat mX, arma::vec vw, arma::mat mD, double trunc) 
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  arma::mat out = mD;
  for (int i = 0; i<p; ++i) {
    if (vw[i]>trunc) {
      arma::vec vx = mX.col(i);
      out = out + (vw[i]*vx) * vx.t();
    }
  }
  return out;
}




// [[Rcpp::export]]
List bayesian_lasso_gibbs_tune_pgtn_c(const arma::mat& mX, const arma::vec& vy, 
                                      arma::vec vbeta_init, double lambda_init, double sigma2_init, 
                                      const double a, const double b, const double u, const double v,  
                                      const int nsamples, const int verbose, const double trunc) {
  const int n = mX.n_rows;
  const int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  const arma::mat& mI_n = eye(n, n);
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mV(maxiter,p);
  arma::vec vb_til(maxiter);
  arma::vec vv_til(maxiter);
  
  const arma::vec one_n = ones(n);
  
  // Initialisation
  arma::vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_init;
  double lambda2 = lambda_init*lambda_init;
  const double a_til = a + 0.5*(n + p);
  double b_til;
  const double u_til = u + 0.5*p;
  double v_til;
  
  arma::mat vmu_til;
  arma::mat vsigma2_til;
  arma::vec vbeta;
  arma::mat mQ;
  arma::mat mQ_inv;
  arma::mat mC;
  arma::mat XDh;
  //arma::mat XYV = join_horiz(join_horiz(mX,vy),vy);
  arma::mat QiXYV;
  arma::mat XTQ;
  arma::mat QX;
  arma::mat XD;
  arma::mat QXD;
  arma::mat mU;
  arma::mat mG;
  
  arma::vec vd_inv;
  arma::vec vd_invsqrt;
  arma::vec vu;
  arma::vec vv;
  arma::vec vw;
  
  //Rcpp::Clock clock;
 
  // Parameters of va|rest which is inverse Gaussian
  arma::vec vmu;
  const arma::vec vlambda = ones(p);
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    
    /**
    //clock.tick("A");
    vd_inv = 1/(lambda2*va);
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = -(mX*vu + randn(n));

    XD = mX.each_row() % vd_inv.t();
    
    if (trunc>0) {
      uvec ids = find( vd_inv>trunc );
      mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
    } else {
      mQ_inv = XD*mX.t() + mI_n;
    }
    
    mQ = inv_sympd(mQ_inv, inv_opts::allow_approx );
    QXD = mQ*XD;
    vmu_til = QXD.t()*vy;
    vw = (sum(XD % QXD,0)).t();
    vsigma2_til = sigma2*(vd_inv - vw);
    vw = QXD.t()*vv;
    vbeta = vmu_til + sqrt(sigma2)*(vu + vw);
    **/
    
    vd_inv = 1/(lambda2*va);
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = mX*vu + randn(n);
    
    XD = mX.each_row() % vd_inv.t();
    if (trunc>0) {
      uvec ids = find( vd_inv>trunc );
      mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
    } else {
      mQ_inv = XD*mX.t() + mI_n;
    }
    
    // Force symmetry
    // Sometimes for numerical reasons it isn't symmetric.
    mQ_inv = 0.5*(mQ_inv + mQ_inv.t());
    
    mU = trimatu(chol(mQ_inv));
    arma::mat UT = trimatl(mU.t()); 
        
    mG = solve(UT, mX);
    vw = (mG % mG).t() * one_n; 
    vsigma2_til = sigma2*(vd_inv - vd_inv%vd_inv%vw);
        
    vw = solve(mU, solve(mU.t(), -vv));
    vmu_til = mX.t() * solve(mU, solve(UT, vy));
    vmu_til = vd_inv % vmu_til;
    vbeta = vmu_til + sqrt(sigma2)*(vu + vd_inv%(mX.t()*vw));

    
    // Sample from sigma2|rest
    b_til = b + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(va%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rinvgaussian_c(vmu, vlambda);
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    // Storing samples
    //Rcout << "A: \n";
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    //Rcout << "B: \n";
    mM.row(i) = vmu_til.as_row();
    //Rcout << "C: \n";
    mV.row(i) = vsigma2_til.as_row();
    //Rcout << "D: \n";
    vb_til[i] = b_til;
    vv_til[i] = v_til;
  }
  
  //clock.stop("Rcpp_times");
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["a_til"] = a_til, 
                      _["vb_til"] = vb_til,
                      _["u_til"] = u_til, 
                      _["vv_til"] = vv_til);
}




// [[Rcpp::export]]
List bayesian_lasso_gibbs_default_c(const arma::mat& mX, const arma::vec& vy, 
                                    arma::vec vbeta_init, double lambda_init, double sigma2_init, 
                                    const double a, const double b, const double u, const double v, 
                                    const int nsamples, const int verbose, const int inv_gauss_type) {
  const int n = mX.n_rows;
  const int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  const arma::mat& mI_n = eye(n, n);
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mV(maxiter,p);
  arma::vec vb_til(maxiter);
  arma::vec vv_til(maxiter);
  
  const arma::vec one_n = ones(n);
  const arma::vec one_p = ones(p);
  
  // Initialisation
  arma::vec vb = ones(p); 
  double sigma2 = sigma2_init;
  double sigma2_hat;
  double lambda2 = lambda_init*lambda_init;
  const double a_til = a + 0.5*n;
  double b_til;
  const double u_til = u + p;
  double v_til;
  
  arma::vec XTy = mX.t() * vy;
  double yTy = sum(vy%vy);
  
  arma::mat vmu_til;
  arma::mat vsigma2_til;
  arma::vec vbeta = vbeta_init;
  arma::mat mQ;
  arma::mat mQ_inv;
  arma::mat XD;
  arma::mat mU;
  arma::mat mR;
  arma::mat mG;
  
  arma::vec vd_inv;
  arma::vec vd_invsqrt;
  arma::vec vu;
  arma::vec vv;
  arma::vec vw;
  
 
  // Parameters of vb|rest which is inverse Gaussian
  arma::vec vmu;
  arma::vec vlambda;
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    //Rcout << "Sample from b|rest \n";
    // Sample from b|rest
    vmu = sqrt(lambda2*sigma2)/abs(vbeta);
    vlambda = lambda2*one_p;
    
    // To avoid numerical problems when calculating mQ_inv
    // we clamp large values of vmu
    vmu.clamp(0.0, 1.0E8);
    
    if (inv_gauss_type==0) {
      vb = rrinvgauss(vmu, vlambda);
    }
    if (inv_gauss_type==1) {
      vb  = rinvgaussian_c(vmu, vlambda);
    }
    if (inv_gauss_type==2) {
      vb = rinvgaussian_slice_c(vb, vmu, vlambda);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    //Rcout << "Sample from lambda2|rest \n";
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(1/vb);
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    ////////////////////////////////////////////////////////////////////////////

    vd_inv = 1/vb;
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = mX*vu + randn(n);
    
    XD = mX.each_row() % vd_inv.t();
    //if (trunc>0) {
    //  uvec ids = find( vd_inv>trunc );
    //  mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
    //} else {
      mQ_inv = XD*mX.t() + mI_n;
    //}
    
    if (!mQ_inv.is_sympd()) {
      // Something has gone horribly wrong.
      vd_inv.print();
      Rcout << "Catastrophic error!!!!  mQ_inv is not positive definite \n";
      return List::create(
                          _["vb"] = vb,
                          _["vmu"] = vmu,
                          _["vlambda"] = vlambda,      
                          _["vd_inv"] = vd_inv,
                          _["mQ_inv"] = mQ_inv);
    }
 
    // Force symmetry
    // Sometimes for numerical reasons it isn't symmetric.
    mQ_inv = 0.5*(mQ_inv + mQ_inv.t());

    mU = trimatu(chol(mQ_inv));
    arma::mat mL = trimatl(mU.t()); 
    mG = solve(mL, mX);
    arma::vec vq = (mG % mG).t() * one_n; 
    vw = solve(mU, solve(mL, -vv));
    vmu_til = XD.t() * solve(mU, solve(mL, vy));

 
    //Rcout << "Sample from sigma2|rest\n";
    // Sample from sigma2|rest
    sigma2_hat = (yTy - sum(XTy%vmu_til))/n;
    b_til = b + 0.5*n*sigma2_hat;
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Calculate diagonal elements of mSigma_til 
    // for Rao-Blackwellization
    vsigma2_til = sigma2*(vd_inv - vd_inv%vd_inv%vq);
    
    //Rcout << "Sample from beta|rest\n";
    // Sample from beta|rest
    vbeta = vmu_til + sqrt(sigma2)*(vu + XD.t()*vw);
    
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = vsigma2_til.as_row();
    vb_til[i] = b_til;
    vv_til[i] = v_til;
  }
  
  //clock.stop("Rcpp_times");
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["a_til"] = a_til, 
                      _["vb_til"] = vb_til,
                      _["u_til"] = u_til, 
                      _["vv_til"] = vv_til);
}





// [[Rcpp::export]]
List bayesian_lasso_gibbs_approx_c(const arma::mat& mX, const arma::vec& vy, 
                                    arma::vec vbeta_init, double lambda_init, double sigma2_init, 
                                    const double a, const double b, const double u, const double v, 
                                    const int nsamples, const int verbose, const int inv_gauss_type,
                                    const double trunc) {
  const int n = mX.n_rows;
  const int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  const arma::mat& mI_n = eye(n, n);
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mV(maxiter,p);
  arma::vec vb_til(maxiter);
  arma::vec vv_til(maxiter);
  
  const arma::vec one_n = ones(n);
  const arma::vec one_p = ones(p);
  
  // Initialisation
  arma::vec vb = ones(p); 
  double sigma2 = sigma2_init;
  double sigma2_hat;
  double lambda2 = lambda_init*lambda_init;
  const double a_til = a + 0.5*n;
  double b_til;
  const double u_til = u + 0.5*p;
  double v_til;
  
  arma::vec XTy = mX.t() * vy;
  double yTy = sum(vy%vy);
  
  arma::mat vmu_til;
  arma::mat vsigma2_til;
  arma::vec vbeta = vbeta_init;
  arma::mat mQ;
  arma::mat mQ_inv;
  arma::mat XD;
  arma::mat mU;
  arma::mat mR;
  arma::mat mG;
  
  arma::vec vd_inv;
  arma::vec vd_invsqrt;
  arma::vec vu;
  arma::vec vv;
  arma::vec vw;
  
  
  // Parameters of vb|rest which is inverse Gaussian
  arma::vec vmu;
  arma::vec vlambda;
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    //Rcout << "Sample from b|rest \n";
    // Sample from b|rest
    //vmu = sqrt(lambda2*sigma2)/abs(vbeta);
    //vlambda = lambda2*one_p;
    
    vmu = sqrt(sigma2/lambda2)/abs(vbeta);
    vlambda = lambda2*one_p;
    
    // To avoid numerical problems when calculating mQ_inv
    // we clamp large values of vmu
    vmu.clamp(0.0, 1.0E8);
    
    if (inv_gauss_type==0) {
      vb = rrinvgauss(vmu, vlambda);
    }
    if (inv_gauss_type==1) {
      vb  = rinvgaussian_c(vmu, vlambda);
    }
    if (inv_gauss_type==2) {
      vb = rinvgaussian_slice_c(vb, vmu, vlambda);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    //Rcout << "Sample from lambda2|rest \n";
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(pow(vbeta,2.0)/vb)/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    ////////////////////////////////////////////////////////////////////////////
    
    vd_inv = 1/(lambda2*vb);
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = mX*vu + randn(n);
    
    XD = mX.each_row() % vd_inv.t();
    if (trunc>0) {
      uvec ids = find( vd_inv>trunc );
      mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
      mQ_inv = (mQ_inv + mQ_inv)/2;
    } else {
      mQ_inv = XD*mX.t() + mI_n;
      mQ_inv = (mQ_inv + mQ_inv)/2;
    }
    
    if (!mQ_inv.is_sympd()) {
      // Something has gone horribly wrong.
      vd_inv.print();
      Rcout << "Catastrophic error!!!!  mQ_inv is not positive definite \n";
      return List::create(
        _["vb"] = vb,
        _["vmu"] = vmu,
        _["vlambda"] = vlambda,      
        _["vd_inv"] = vd_inv,
        _["mQ_inv"] = mQ_inv);
    }
    
    // Force symmetry
    // Sometimes for numerical reasons it isn't symmetric.
    mQ_inv = 0.5*(mQ_inv + mQ_inv.t());
    
    mU = trimatu(chol(mQ_inv));
    arma::mat mL = trimatl(mU.t()); 
    mG = solve(mL, mX);
    arma::vec vq = (mG % mG).t() * one_n; 
    vw = solve(mU, solve(mL, -vv));
    vmu_til = XD.t() * solve(mU, solve(mL, vy));
    
    //Rcout << "Sample from sigma2|rest\n";
    // Sample from sigma2|rest
    sigma2_hat = (yTy - sum(XTy%vmu_til))/n;
    b_til = b + 0.5*n*sigma2_hat;
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Calculate diagonal elements of mSigma_til 
    // for Rao-Blackwellization
    vsigma2_til = sigma2*(vd_inv - vd_inv%vd_inv%vq);
    
    //Rcout << "Sample from beta|rest\n";
    // Sample from beta|rest
    vbeta = vmu_til + sqrt(sigma2)*(vu + XD.t()*vw);
    
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = vsigma2_til.as_row();
    vb_til[i] = b_til;
    vv_til[i] = v_til;
  }
  
  //clock.stop("Rcpp_times");
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["a_til"] = a_til, 
                      _["vb_til"] = vb_til,
                      _["u_til"] = u_til, 
                      _["vv_til"] = vv_til);
}



