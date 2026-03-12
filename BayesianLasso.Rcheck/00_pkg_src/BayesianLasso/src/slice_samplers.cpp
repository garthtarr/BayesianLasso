

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

////////////////////////////////////////////////////////////////////////////////

// The log of the target distribution
double f_fun_precision(double x, double a_val, double b_val, double c_val) {
  double val = a_val*log(x) - b_val*x - c_val*sqrt(x);
  return val;
}

// It's gradient
double g_fun_precision(double x, double a_val, double b_val, double c_val) {
  double val = a_val/x - b_val - 0.5*c_val/sqrt(x);
  return val;
}

// It's Hessian
double h_fun_precision(double x, double a_val, double b_val, double c_val) {
  double val = -a_val/(x*x)  + 0.25*c_val/(sqrt(x)*x);
  return val;
}

/*
 * This is a slice sampler for the precision parameter and penalty parameter
 * for the Hans (2009) Gibbs sampler using lambda2 ~ Gamma(u,v) instead of
 * lambda ~ Gamma(u,v). Note we use this prior to make a more appropriate
 * comparison with other Gibbs samplers.
 */


double slice_sampler_precision(double x, double a_val, double b_val, double c_val)
{
  double TOL = 1.0E-8;
  double MAXITER = 100;

  // Calculate the mode
  double z_star = (-0.5*c_val + sqrt(0.25*c_val*c_val + 4*a_val*b_val))/(2*b_val);
  double x_star = z_star*z_star;

  // Target at the mode
  double f_star = f_fun_precision(x_star, a_val, b_val, c_val);
  // double p_star = exp(f_star);

  // Calculate approximate variaiance and sd
  double h_star = h_fun_precision(x_star, a_val, b_val, c_val);
  double sigma2 = -1/h_star;
  // double sigma  = sqrt(sigma2);

  // Calculate target at initial x
  double x0 = x;
  double f0 = f_fun_precision(x0, a_val, b_val, c_val) - f_star;
  double p0 = exp(f0);

  if (p0==0.0) {
    // Something has gone horribly wrong because the
    // initial x is far outside the effective domain
    x0 = x_star;
    f0 = f_fun_precision(x0, a_val, b_val, c_val) - f_star;
    p0 = exp(f0);
  }

  // Sample uniform vertically
  double u = R::runif(0,p0);
  double log_u = log(u);

  // Initial guesses of the left and right end-points based on Laplace approximation
  double xL = x_star - sqrt(2*sigma2*( - log_u));
  double xR = x_star + sqrt(2*sigma2*( - log_u));

  // Set left end-point initial guess to be small if Laplace approximation
  // gives negative end-point
  if (xL<0) {
    xL = 1.0E-12;
  }

  double fL;
  double gL;
  double fR;
  double gR;

  // Use Newton's method to find left end point
  for (int i=0; i<MAXITER; ++i) {
    fL = f_fun_precision(xL, a_val, b_val, c_val) - f_star;
    gL = g_fun_precision(xL, a_val, b_val, c_val);
    xL = xL - (fL - log_u)/gL;
    if (abs(gL)<TOL) { break; }
  }

  // Use Newton's method to find right end-point
  for (int i=0; i<MAXITER; ++i) {
    fR = f_fun_precision(xR, a_val, b_val, c_val) - f_star;
    gR = g_fun_precision(xR, a_val, b_val, c_val);
    xR = xR - (fR - log_u)/gR;
    if (abs(gR)<TOL) { break; }
  }

  // Sample uniformly between left and right end-points
  double x_new = R::runif(xL, xR);

  return x_new;
}

////////////////////////////////////////////////////////////////////////////////

void print_inputs(double x, arma::vec vs2, arma::vec vc2, double a_til,
                        double b, double u_til, double v, double yTy,
                        double lower, double upper, double a_val, double b_val,
                        double c_val, double ga, double gb, double gc)
{
  Rcout << " a_val: " << a_val << " b_val: " << b_val << " c_val: " << c_val << " ga: " << ga << " gb: " << gb << " gc: " << gc << "\n";
  stop("Error!");
   /**
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
   **/
}

////////////////////////////////////////////////////////////////////////////////

double f_fun_orthogonal(double lambda2, arma::vec vs2, arma::vec vc2, double a_til,
             double b, double u_til, double v, double yTy)
{
  arma::vec vw = 1/(vs2 + lambda2);
  double val = u_til*log(lambda2) - v*lambda2 + 0.5*sum(log(vw)) - a_til*log(b + 0.5*(yTy - dot(vc2,vw)));
  return val;
}

double g_fun_orthogonal(double lambda2, arma::vec vs2, arma::vec vc2, double a_til,
             double b, double u_til, double v, double yTy)
{
  arma::vec vw = 1/(vs2 + lambda2);
  double val = u_til/lambda2 - v - 0.5*sum(vw)  - a_til*dot(vc2,square(vw))/(2.0*b + yTy - dot(vc2,vw));
  return val;
}


double slice_orthogonal(double x, arma::vec vs2, arma::vec vc2, double a_til,
                            double b, double u_til, double v, double yTy,
                            double lower, double upper)
{
  double tol = 1.0E-3;
  int  maxiter = 100;
  int num_fun_evals = 0;

  double a_val = lower;
  double b_val = upper;
  double c_val = (a_val + b_val)/2.0;
  double ga = g_fun_orthogonal(a_val, vs2, vc2, a_til, b, u_til, v, yTy); num_fun_evals ++;
  double gb = g_fun_orthogonal(b_val, vs2, vc2, a_til, b, u_til, v, yTy); num_fun_evals ++;
  double gc = g_fun_orthogonal(c_val, vs2, vc2, a_til, b, u_til, v, yTy); num_fun_evals ++;

  if (sign(gb) == sign(ga)) {
    Rcout << "Error in slice_orthogonal: bisection search signs are equal \n";
    print_inputs(x, vs2, vc2, a_til, b, u_til, v, yTy, lower, upper, a_val, b_val, c_val, ga, gb, gc);
  }

  double x_star;

  for (int i = 0; i < maxiter; ++i)
  {
    c_val = (a_val + b_val)/2.0;
    gc = g_fun_orthogonal(c_val, vs2, vc2, a_til, b, u_til, v, yTy); num_fun_evals ++;
    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      x_star = c_val;
      break;
    }

    if (i==(maxiter-1)) {
      Rcout << "Error in slice_orthogonal: maximum iterations reached \n";
      print_inputs(x, vs2, vc2, a_til, b, u_til, v, yTy, lower, upper, a_val, b_val, c_val, ga, gb, gc);
    }

    if (sign(gc) == sign(ga)) {
      a_val = c_val;
      ga = gc;
    } else {
      b_val = c_val;
      gb = gc;
    }
  }

  double f_star = f_fun_orthogonal(x_star, vs2, vc2, a_til, b, u_til, v, yTy); num_fun_evals ++;
  // double p_star = exp(f_star);

  double x0 = x;
  double f0 = f_fun_orthogonal(x0, vs2, vc2, a_til, b, u_til, v, yTy)  - f_star; num_fun_evals ++;
  double p0 = exp(f0);

  if (p0<tol) {
    // Something has gone horribly wrong because the initial x is far outside
    // the effective domain
    x0 = x_star;
    f0 = f_fun_orthogonal(x0, vs2, vc2, a_til, b, u_til, v, yTy)  - f_star; num_fun_evals ++;
    p0 = exp(f0);
  }

  if (x0<=0) {
    // Something has gone horribly wrong
    x0 = x_star;
    f0 = f_fun_orthogonal(x0, vs2, vc2, a_til, b, u_til, v, yTy)  - f_star; num_fun_evals ++;
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
  ga = f_fun_orthogonal(a_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
  gb = f_fun_orthogonal(b_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
  for (int i = 0; i < maxiter; ++i)
  {
    c_val = (a_val + b_val)/2;
    gc = f_fun_orthogonal(c_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      xL = c_val;
      break;
    }

    if (i==(maxiter-1)) {
      Rcout << "Error in slice_orthogonal: maximum iterations reached \n";
      print_inputs(x, vs2, vc2, a_til, b, u_til, v, yTy, lower, upper, a_val, b_val, c_val, ga, gb, gc);
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
  ga = f_fun_orthogonal(a_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
  gb = f_fun_orthogonal(b_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
  for (int i = 0; i < maxiter; ++i)
  {
    c_val = (a_val + b_val)/2;
    gc = f_fun_orthogonal(c_val, vs2, vc2, a_til, b, u_til, v, yTy) - f_star - log_u; num_fun_evals ++;
    if ((abs(gc)<tol)||(abs(b_val - a_val)<tol)) {
      xR = c_val;
      break;
    }

    if (i==(maxiter-1)) {
      Rcout << "Error in slice_orthogonal: maximum iterations reached \n";
      print_inputs(x, vs2, vc2, a_til, b, u_til, v, yTy, lower, upper, a_val, b_val, c_val, ga, gb, gc);
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

  return x_new;

  //List::create(
  //  _["error"] = false,
  //  _["x_new"] = x_new,
  //  _["num_fun_evals"] = num_fun_evals);
}

////////////////////////////////////////////////////////////////////////////////

inline double f_fun_pcg(double x, double a_val, double b_val, double c_val,
                    double d_val, double e_val, double log_kappa = 0) {
  double val = a_val*log(x) - b_val*x - c_val*log(d_val + e_val*x) - log_kappa;
  return val;
}

inline double p_fun_pcg(double x, double a_val, double b_val, double c_val,
                    double d_val, double e_val) {
  double val = exp(f_fun_pcg(x, a_val, b_val, c_val, d_val, e_val));
  return val;
}

inline double g_fun_pcg(double x, double a_val, double b_val, double c_val,
                    double d_val, double e_val) {
  return a_val/x - b_val -  c_val*e_val/(d_val + e_val*x);
}

inline double h_fun_pcg(double x, double a_val, double b_val, double c_val,
                    double d_val, double e_val) {
  double val = e_val/(d_val + e_val*x);
  return -a_val/(x*x) + c_val*(val*val);
}


double solve_quadratic_pcg(double a_til, double b_til, double c_til)
{
  // Compute discriminant
  double D = b_til * b_til - 4 * a_til * c_til;

  if (D < 0) {
    stop("Error in solve_quadratic: discriminant is negative");
  }

  double sqrtD = sqrt(D);
  double x1, x2;

  // Choose stable form for x1
  if (b_til >= 0) {
    x1 = (-b_til - sqrtD) / (2 * a_til);
  } else {
    x1 = (-b_til + sqrtD) / (2 * a_til);
  }

  // Compute x2 using product of roots
  x2 = c_til / (a_til * x1);

  double bigger = 0;
  if (x1 > x2) {
    bigger = x1;
  } else {
    bigger = x2;
  }

  return bigger;
}


double mode_p_fun_pcg(double a_val, double b_val, double c_val, double d_val,
                  double e_val)
{
  double A_val = b_val*e_val;
  double B_val = b_val*d_val + e_val*(c_val - a_val);
  double C_val = -a_val*d_val;
  double val = solve_quadratic_pcg(A_val, B_val, C_val);

  return val;
}


/**
 * Solve f_fun(x, a_val, b_val, c_val, d_val, e_val, kappa) = 0 using Newton
 * steps where possible and the bisection method as a fallback when Newton steps
 * return nonsense.
 *
 * x_init is an initial guess
 * a_val, b_val, c_val, d_val, e_val are parameters of f
 * kappa is a crossing point
 * xL is is a left bracket for the solution
 * xR is is a right bracket for the solution
 * left is a boolean for whether we are looking for the left solution
 * maxiter is the maximum number of iterations (100 is overkill)
 * tol - terminate when the absolute value of the stepsize is smaller than tol
 */


double solve_pcg(double x_init,
                 double a_val,
                 double b_val,
                 double c_val,
                 double d_val,
                 double e_val,
                 double log_kappa,
                 double xL,
                 double xR,
                 bool left = true,
                 int maxiter = 100,
                 double tol = 1.0E-8)
{
  // Looking for positive solutions
  // If the input is nonsense use the midpoint
  if (x_init <= 0) {
    x_init = 0.5 * (xL + xR);
  }

  // Update (x,f)
  double x = x_init;
  double f = f_fun_pcg(x, a_val, b_val, c_val, d_val, e_val, log_kappa);

  // Update left and right brackets
  if (left) {
    if (f < 0) xL = x; else xR = x;
  } else {
    if (f < 0) xR = x; else xL = x;
  }

  for (int iter = 0; iter < maxiter; iter++)
  {
    // Calculate Newton step
    double g = g_fun_pcg(x, a_val, b_val, c_val, d_val, e_val);
    double delta = f / g;

    // Check whether Newton step is positive
    double x_temp = x - delta;
    if (x_temp <= 0) {
      x = 0.5 * (xL + xR);
    } else {
      x = x_temp;
    }

    // Update left and right brackets
    f = f_fun_pcg(x, a_val, b_val, c_val, d_val, e_val, log_kappa);
    if (left) {
      if (f < 0) xL = x; else xR = x;
    } else {
      if (f < 0) xR = x; else xL = x;
    }

    if (std::abs(delta) < tol) break;
  }

  return x;
}

/**
 * Samples from p_fun (defined above) using slice sampling.
 */


double slice_pcg(double x_curr,
                 double a_val,
                 double b_val,
                 double c_val,
                 double d_val,
                 double e_val,
                 double nsd = 20,
                 int maxiter = 100,
                 double tol = 1.0E-8)
{
  // Compute the mode
  double x_mode = mode_p_fun_pcg(a_val, b_val, c_val, d_val, e_val);
  double f_mode = f_fun_pcg(x_mode, a_val, b_val, c_val, d_val, e_val);

  // Get current value at un-normalized
  double f_curr = f_fun_pcg(x_curr, a_val, b_val, c_val, d_val, e_val);

  double p_curr = exp(f_curr - f_mode);
  if (p_curr < 1.0E-8) {
    // If outside the effective domain for some reason reset
    x_curr = x_mode;
    p_curr = 1.0;
  }

  double u_val = R::runif(0, 1);
  double log_u_val = log(u_val);
  double log_kappa = f_curr + log_u_val;

  // Compute the second derivative at the mode
  double h = h_fun_pcg(x_mode, a_val, b_val, c_val, d_val, e_val);

  // Use Laplace approximation as a starting point
  double sigma2 = -1.0 / h;
  double sigma = sqrt(sigma2);

  // Get approximate starting points using normal approximation
  double term = 2 * sigma2 * (f_mode - log_kappa);
  double sqrt_term = sqrt(term);

  // Two possible solutions
  double x1 = x_mode + sqrt_term;
  double x2 = x_mode - sqrt_term;
  double L = std::min(x1, x2);
  double R = std::max(x1, x2);

  // Ensure left starting value is positive
  if (L < 0) L = x_mode / 2;

  // Solve for left and right bounds
  double xL = solve_pcg(L, a_val, b_val, c_val, d_val, e_val, log_kappa, 1.0E-14, x_mode, true, maxiter, tol);
  double xR = solve_pcg(R, a_val, b_val, c_val, d_val, e_val, log_kappa, x_mode, x_mode + nsd * sigma, false, maxiter, tol);

  // Sample a new y uniformly between xL and xR
  double x_new = R::runif(xL, xR);

  return x_new;
}

////////////////////////////////////////////////////////////////////////////////

// Slice sampler for p(z) \propto exp(-alpha*z)/(1+z)
// given the current value z_curr and rate alpha.
// Returns a new sample of z.


double slice_sample_z(double z_curr, double alpha)
{
  // Step 1: Sample u1 ~ Uniform(0, 1/(1 + z_curr))
  double u1 = R::runif(0.0, 1.0 / (1.0 + z_curr));

  // Step 2: Compute the truncation boundary:
  // u1 < 1/(1+z)  <=>  z < (1-u1)/u1.
  double z_max = (1.0 - u1) / u1;

  // Step 3: Sample z from an exponential RV with rate alpha
  //         truncated between (0, z_max).
  //
  // Note: The CDF for Exp(alpha) is
  //    F(z) = 1 - exp(-alpha*z)
  // so the truncated CDF is:
  //    F_trunc(z) = (1 - exp(-alpha*z)) / (1 - exp(-alpha*z_max))
  // Inverting:
  //    z = -log(1 - u2*(1 - exp(-alpha*z_max)))/alpha,
  // where u2 ~ Uniform(0,1).
  double u2 = R::runif(0.0, 1.0);
  double z_new = -log(1 - u2 * (1 - exp(-alpha * z_max))) / alpha;

  return z_new;
}

// This function performs one slice sampling update for eta,
// where the target for eta is:
//    p(eta) \propto exp(-alpha * eta^2)*(eta/(1+eta^2)).
// We transform via z = eta^2.


double sample_eta(double eta_curr, double alpha)
{
  // Transform to z:
  double z_curr = eta_curr * eta_curr;

  // Slice-sample a new z:
  double z_new = slice_sample_z(z_curr, alpha);

  // Transform back:
  double eta_new = sqrt(z_new);
  return eta_new;
}


// This function generates a series of samples from the posterior for eta.
// n: number of samples, init_eta: starting value for eta,
// alpha: parameter (typically mu_j^2/2).


arma::vec sample_eta_vec(int n, double eta_init, double alpha)
{
  arma::vec veta(n);
  double eta_curr = eta_init;
  for (int i = 0; i < n; i++) {
    eta_curr = sample_eta(eta_curr, alpha);
    veta[i] = eta_curr;
  }
  return veta;
}


