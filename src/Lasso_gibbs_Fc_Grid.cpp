
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "lasso_distribution.h"

// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

/////////////////


// [[Rcpp::export]]
List svd_cpp(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, X, "both","dc");
  return List::create(_["mU"] = U,
                      _["vs"] = S,
                      _["mV"] = V);
}

// [[Rcpp::export]]
List svd_left_cpp(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, X, "left","dc");
  return List::create(_["mU"] = U,
                      _["vs"] = S,
                      _["mV"] = V);
}

// [[Rcpp::export]]
List svd_base(const NumericMatrix & X) {
  Environment base("package:base");
  Function svd = base["svd"];
  List svdResult = svd(X);
  return svdResult;
}

///////////////////////////////////

// [[Rcpp::export]]
List demmler_reinsch(arma::mat mX, arma::vec vy, arma::vec va, double sigma2, arma::vec vlambda2) 
{
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  double yTy = sum(vy%vy);
  
  arma::vec sqrt_va = sqrt(va);
  arma::vec sqrt_va_inv = 1.0/sqrt_va;
  arma::mat mX1_til = mX.each_row() % sqrt_va_inv.as_row();
  
  
  arma::mat mU, mV;
  arma::vec vs;
  svd(mU, vs, mV, mX1_til, "dc");
  
  
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  
  Rcout << "n = " << vy.n_elem << "\n";
  Rcout << "p = " << mX.n_cols << "\n";
  
  Rcout << "mU.n_rows = " << mU.n_rows << "\n";
  Rcout << "mU.n_cols = " << mU.n_cols << "\n";
  Rcout << "vs.n_elem = " << vs.n_elem << "\n";
  Rcout << "mV.n_rows = " << mV.n_rows << "\n";
  Rcout << "mV.n_cols = " << mV.n_cols << "\n";
  
  arma::mat mS(n,p);
  mS.diag() = vs;
  
  
  arma::mat mX1_til_hat = mU * mS * mV.t();
  
  arma::mat mD = mX1_til_hat - mX1_til;  
  
  Rcout << max(abs(mD)) << "\n";
  
  
  arma::vec UTy = mU.t()*vy;
  arma::vec vc = vs%UTy;
  arma::vec vc2 = vc%vc;
  arma::vec vs2 = vs%vs;
  
  int N = vlambda2.n_elem;
  
  
  arma::vec vsigma2_hat(N);   
  arma::vec vlogdet_mSigma(N); 
  
  arma::mat mMu(N,p);
  arma::mat mSigma(N,p);
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  
  for (int k = 0; k < N; ++k) 
  {
    arma::vec vw(p);
    for (int i = 0; k < n; ++k) {
      vw[i] = 1.0/(vs2[i] + vlambda2[k]);
    }
    for (int i = 0; i < (p-n); ++k) {
      vw[i+n] = 1.0/vlambda2[k];
    }
    

    arma::mat mB = diagmat(sqrt_va_inv);
    arma::mat mW = diagmat(vw);
    
 
    arma::mat mQ = mB * mV * mW * mV.t() * mB;
    arma::vec vsigma2 = diagvec(mQ);
    
    arma::mat mSigma_k = inv(XTX + diagmat(vlambda2[k]*va));
    arma::vec vsigma2_k = diagvec(mSigma_k);
    
    Rcout << "Demmler-Reinsch: \n";
    vsigma2.as_row().print();
    Rcout << "Naive approach: \n";
    
    vsigma2_k.as_row().print();
    stop("blah");
    
   }  
  
  return List::create(_["vsigma2_hat"] = vsigma2_hat,
                      _["vlogdet_mSigma"] = vlogdet_mSigma,
                      _["mMu"] = mMu,
                      _["mSigma"] = mSigma);  
}

///////////////////////////////////

// [[Rcpp::export]]
List calc_penreg_stats(arma::mat mX, arma::vec vy, arma::vec va, double sigma2, arma::vec vlambda2)  
{
  int N = vlambda2.n_elem;
  int n = vy.n_elem;
  int p = mX.n_cols; 
  
  arma::mat XTX = mX.t()*mX;
  arma::vec XTy = mX.t()*vy;
  double yTy = sum(vy%vy);
  
  arma::vec vlogdet_mSigma(N);
  arma::vec vsigma2_hat(N);      
  arma::mat mMu(N,p);       // vmu matrix for each vlambda2 value
  arma::mat mSigma(N,p);
  
  double val; // for log_det call
  double sign; // for log_det call
  bool ok;
  
  for (int k = 0; k < N; ++k) {
    arma::vec vmu = inv(XTX + diagmat(vlambda2[k]*va))*XTy;
    ok = log_det(val, sign, XTX + diagmat(vlambda2[k]*va));
    vsigma2_hat[k] = (yTy - sum(vmu % XTy))/n;
    vlogdet_mSigma[k] = p*log(sigma2) - val;
    
    arma::mat mSigma_k = sigma2*inv(XTX + diagmat(vlambda2[k]*va));
    arma::vec vsigma2 = diagvec(mSigma_k);
    
    mMu.row(k) = vmu.as_row();
    mSigma.row(k) = vsigma2.as_row();
  }  
  
  return List::create(_["vsigma2_hat"] = vsigma2_hat,
                      _["vlogdet_mSigma"] = vlogdet_mSigma,
                      _["mMu"] = mMu,
                      _["mSigma"] = mSigma); 
}

///////////////////////////////////

double f_fun(double lambda2, arma::vec vs2, arma::vec vc2, double a_til, double b, double u_til, double v, double yTy) 
{
  arma::vec vw = 1/(vs2 + lambda2);
  double val = u_til*log(lambda2) - v*lambda2 + 0.5*sum(log(vw)) - a_til*log(b + 0.5*(yTy - sum(vc2%vw)));
  return val;
}

double g_fun(double lambda2, arma::vec vs2, arma::vec vc2, double a_til, double b, double u_til, double v, double yTy) 
{
  arma::vec vw = 1/(vs2 + lambda2);
  arma::vec vw2 = vw%vw;
  double val = u_til/lambda2 - v - 0.5*sum(vw)  - a_til*sum(vc2%vw2)/(2*b + yTy - sum(vc2%vw));
  return val;
}

List slice_sampler_svd(double x, arma::vec vs2, arma::vec vc2, double a_til, double b, double u_til, double v, double yTy, double lower, double upper) 
{
  double tol = 1.0E-3;
  int  maxiter = 100;
  
  
  double a_val = lower;
  double b_val = upper;
  double c_val = (a_val + b_val)/2.0;
  double ga = g_fun(a_val, vs2, vc2, a_til, b, u_til, v, yTy); 
  double gb = g_fun(b_val, vs2, vc2, a_til, b, u_til, v, yTy);
  double gc = g_fun(c_val, vs2, vc2, a_til, b, u_til, v, yTy);
  
  if (sign(gb) == sign(ga)) {
    Rcout << "Error in slice_sampler_svd: bisection search signs are equal \n";
    Rcout << " a_val: " << a_val << " b_val: " << b_val << " c_val: " << c_val << " ga: " << ga << " gb: " << gb << " gc: " << gc << "\n";
    return List::create(
      _["error"] = true,
      _["x"] = x,
      _["vs2"] = vs2,
      _["vc2"] = vc2,
      _["b"] = b,
      _["u_til"] = u_til,
      _["v"] = v,
      _["yTy"] = yTy,
      _["lower"] = lower,
      _["upper"] = upper);
  }
  
  double x_star;
  
  for (int i = 0; i < maxiter; ++i) 
  {
    c_val = (a_val + b_val)/2.0;
    gc = g_fun(c_val, vs2, vc2, a_til, b, u_til, v, yTy);

    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      x_star = c_val;
      break;
    }
    
    if (i==(maxiter-1)) {
      Rcout << "Error in slice_sampler_svd: maximum iterations reached \n";
      Rcout << " a_val: " << a_val << " b_val: " << b_val << " c_val: " << c_val << " ga: " << ga << " gb: " << gb << " gc: " << gc << "\n";
      return List::create(
        _["error"] = true,
        _["x"] = x,
        _["vs2"] = vs2,
        _["vc2"] = vc2,
        _["b"] = b,
        _["u_til"] = u_til,
        _["v"] = v,
        _["yTy"] = yTy,
        _["lower"] = lower,
        _["upper"] = upper);
    }
    
    
    if (sign(gc) == sign(ga)) {
      a_val = c_val;
      ga = gc;
    } else {
      b_val = c_val;
      gb = gc;
    }
  }
  
  double f_star = f_fun(x_star, vs2, vc2, a_til, b, u_til, v, yTy);
  double p_star = exp(f_star);
  
 
  
  double x0 = x;
  double f0 = f_fun(x0, vs2, vc2, a_til, b, u_til, v, yTy)  - f_star;
  double p0 = exp(f0);
  
  if (p0<tol) {
    // Something has gone horribly wrong because the 
    //  initial x is far outside the effective domain
    x0 = x_star;
    f0 = f_fun(x0, vs2, vc2, a_til, b, u_til, v, yTy)  - f_star;
    p0 = exp(f0);
  }
  
  // Sample uniform vertically
  double u = as_scalar(randu(1, distr_param(0.0,p0)));
  double log_u = log(u);
  
  double xL;
  double xR;
  
  // Bisection method to find left point
  a_val = lower;
  b_val = x_star; 
  ga = f_fun(a_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
  gb = f_fun(b_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
  for (int i = 0; i < maxiter; ++i) 
  {
    c_val = (a_val + b_val)/2;
    gc = f_fun(c_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
    
 
    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      xL = c_val;
      break;
    }
    
    if (i==(maxiter-1)) {
      Rcout << "Error in slice_sampler_svd: maximum iterations reached \n";
      Rcout << " a_val: " << a_val << " b_val: " << b_val << " c_val: " << c_val << " ga: " << ga << " gb: " << gb << " gc: " << gc << "\n";
      return List::create(
        _["error"] = true,
        _["x"] = x,
        _["vs2"] = vs2,
        _["vc2"] = vc2,
        _["b"] = b,
        _["u_til"] = u_til,
        _["v"] = v,
        _["yTy"] = yTy,
        _["lower"] = lower,
        _["upper"] = upper);
    }    
    
    
    if (sign(gc) == sign(ga)) {
      a_val = c_val;
      ga = gc;
    } else {
      b_val = c_val;
      gb = gc;
    }
  }
  
  // Bisection method to find right point
  a_val = x_star;
  b_val = upper;
  ga = f_fun(a_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
  gb = f_fun(b_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
  for (int i = 0; i < maxiter; ++i) 
  {
    c_val = (a_val + b_val)/2;
    gc = f_fun(c_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u;
    
     
    
    
    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      xR = c_val;
      break;
    }
    
    if (i==(maxiter-1)) {
      Rcout << "Error in slice_sampler_svd: maximum iterations reached \n";
      Rcout << " a_val: " << a_val << " b_val: " << b_val << " c_val: " << c_val << " ga: " << ga << " gb: " << gb << " gc: " << gc << "\n";
      return List::create(
        _["error"] = true,
        _["x"] = x,
        _["vs2"] = vs2,
        _["vc2"] = vc2,
        _["b"] = b,
        _["u_til"] = u_til,
        _["v"] = v,
        _["yTy"] = yTy,
        _["lower"] = lower,
        _["upper"] = upper);
    }
    
    
    if (sign(gc) == sign(ga)) {
      a_val = c_val;
      ga = gc;
    } else {
      b_val = c_val;
      gb = gc;
    }
  }
  
  // Sample uniformly between left and right end-points
  double x_new =  as_scalar(randu(1, distr_param(xL,xR)));   
  
  return List::create(
    _["error"] = false,
    _["x_new"] = x_new);
}

///////////////////////////////////
// Sampler F

// [[Rcpp::export]]
List lasso_gibbs_Fc(arma::mat mX, arma::vec vy, double a, double b, double u, double v, 
                    int nsamples, arma::vec vbeta_init, double lambda_init, double sigma2_init,
                    arma::vec va_init, int verbose, double lower, double upper, int N)
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  arma::vec np(2);
  np(0) = n;
  np(1) = p;
  int min_np = (int) np.min();
  int max_np = (int) np.max();
  
  bool debug = false;
  
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples for MCMC
  arma::mat mBeta(maxiter,p);
  arma::vec vsigma2(maxiter);
  arma::vec vlambda2(maxiter);
  // Initialize storage for Rao-Blackwellization
  arma::mat mM(maxiter,p);
  arma::mat mS(maxiter,p);
  arma::vec vb_til(maxiter);
  arma::vec vv_til(maxiter);
  
  arma::vec va = va_init;
  
  // Assign initial values
  arma::vec vbeta = vbeta_init;
  
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda*lambda;
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  

  arma::vec XTy = mX.t() * vy;

  double yTy = sum(vy%vy);
  
  arma::vec vmu_til(p);
  arma::vec vsigma2_til(p);
  
  // Constant values 
  const double a_til = a + 0.5*n;
  double b_til = b;
  const double u_til = u + 0.5*p;
  double v_til = v;
  
  arma::mat mI_n;
  
  if(p>n){
    mI_n = eye(n, n);
  }
  else{
    mI_n = eye(p, p);
  }
  
  
  
  
  for (int i = 0; i < maxiter; ++i) 
  {
    arma::vec sqrt_va_inv = sqrt(1/va);
    arma::mat mX1_til = mX.each_row() % sqrt_va_inv.as_row();
    
    arma::mat mU, mV;
    arma::vec vs;
    svd_econ(mU, vs, mV, mX1_til, "left", "dc");
    
    arma::vec UTy = mU.t()*vy;
    arma::vec vc = vs%UTy;
    arma::vec vc2 = vc%vc;
    arma::vec vs2 = vs%vs;

    arma::vec vw;
    

    double u_star = u + 0.5*min_np - 1.0;
    List res_slice = slice_sampler_svd(lambda2, vs2, vc2, a_til, b, u_star, v, yTy, lower, upper);
    if (res_slice["error"]) {
      Rcout << "Error from  slice_sampler_svd \n";
      return res_slice;
    } else {
      lambda2 = res_slice["x_new"];
    }
    vw = 1.0/(vs2 + lambda2);
    double sigma2_hat = (yTy - sum(vc2%vw))/n;
    b_til = b + 0.5*n*sigma2_hat;

    
    //Sampling from sig^2|D,u,lam^2,a ~ IG()
 
    sigma2 = 1.0/randg(distr_param(a_til,1.0/b_til));
    

    //Sample from vbeta|D,sig^2,lam^2,a ~ N()

     if (false) {
      arma::vec vd_inv = 1/(lambda2*va);

      arma::vec vd_invsqrt = sqrt(vd_inv);
      arma::vec vu = randn(p) % vd_invsqrt;
      arma::vec vv = mX*vu + randn(n);
      
      arma::mat XD = mX.each_row() % vd_inv.t();

      arma::mat mQ_inv = XD*mX.t() + mI_n;

      
      if (!mQ_inv.is_sympd()) {
        // Something has gone horribly wrong.
        vd_inv.print();
        Rcout << "Catastrophic error!!!!  mQ_inv is not positive definite \n";
        return List::create(
          _["va"] = va,
          _["sigma2"] = sigma2,
          _["vlambda"] = vlambda2,      
          _["vd_inv"] = vd_inv,
          _["mQ_inv"] = mQ_inv);
      }
      
      // Force symmetry
      // Sometimes for numerical reasons it isn't symmetric.
      mQ_inv = 0.5*(mQ_inv + mQ_inv.t());
      
      mU = trimatu(chol(mQ_inv));
      arma::mat mL = trimatl(mU.t()); 
      arma::mat mG = solve(mL, mX);
      arma::vec vq = (mG % mG).t() * one_n; 
      arma::vec vw = solve(mU, solve(mL, -vv));
      vmu_til = XD.t() * solve(mU, solve(mL, vy));
      
      // Calculate diagonal elements of mSigma_til 
      // for Rao-Blackwellization
      vsigma2_til = sigma2*(vd_inv - vd_inv%vd_inv%vq);

      // Sample from beta|rest
      vbeta = vmu_til + sqrt(sigma2)*(vu + XD.t()*vw);    
    } else {
      double sigma = sqrt(sigma2);
      double lambda = sqrt(lambda2);

      arma::vec vu = (sigma/lambda)*(sqrt_va_inv%randn(p));

      arma::vec vdelta = randn(n);
      arma::vec vv = (1/sigma)*(mX*vu) + vdelta;
      

      arma::vec vw = mU*((1.0/(vs2/lambda2 + 1.0))%(mU.t()*(vy/sigma - vv)));

      vbeta = vu + (sigma/lambda2)*(1/va)%(mX.t()*vw);

      
    }
    
    
    //Sample from va|beta,sig^2,lam^2 ~ IG()
    arma::vec vnu = sqrt(sigma2/lambda2)/abs(vbeta);
    va = rinvgaussian_c(vnu, one_p);

    
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    // Store MCMC samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage sufficient statistics for Rao-Blackwellization for sI
    mM.row(i) = vmu_til.as_row();
    mS.row(i) = vsigma2_til.as_row();
    vb_til[i] = b_til;
    vv_til[i] = v_til;
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mS,
                      _["a_til"] = a_til, 
                      _["vb_til"] = vb_til,
                      _["u_til"] = u_til, 
                      _["vv_til"] = vv_til);
}

