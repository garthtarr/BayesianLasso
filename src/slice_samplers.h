#ifndef SLICE_SAMPLERS_H
#define SLICE_SAMPLERS_H

// Slice sampler used in modified Hans Gibbs sampler
double slice_sampler_precision(double x, double a_val, double b_val, double c_val);

// Slice sampler used by horseshoe auxiliary parameters. Original by Polson et al (2014).
double sample_eta(double lambda2, double mu);

// Since sampler used by partially collapsed Gibbs samplers
double slice_pcg(double x_curr, double a_val, double b_val, double c_val, 
                 double d_val, double e_val, double nsd = 20, int maxiter = 100, double tol = 1.0E-8);


// Slice sampler for lambda2
double slice_orthogonal(double x, arma::vec vs2, arma::vec vc2, double a_til, 
                            double b, double u_til, double v, double yTy, 
                            double lower, double upper);

#endif


 



 


 


