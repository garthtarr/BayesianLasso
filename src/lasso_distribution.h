
#ifndef LASSO_DISTRIBUTION_H
#define LASSO_DISTRIBUTION_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////////////////////////////////////

double expit_c(double x);
//Rcpp::List calculate_lasso_dist_stats_c(double a_val, double b_val, double c_val);
//double logSumExp_c(vec vx);
double zlasso_c_v1(double a_val, double b_val, double c_val, bool logarithm);
double zlasso_c_v2(double a_val, double b_val, double c_val, bool logarithm);
arma::vec dlasso_c_v1(arma::vec vx, double a_val, double b_val, double c_val, bool logarithm);
arma::vec dlasso_c_v2(arma::vec vx, double a_val, double b_val, double c_val, bool logarithm);
arma::vec plasso_c_v1(arma::vec vx, double a_val, double b_val, double c_val);
arma::vec plasso_c_v2(arma::vec vx, double a_val, double b_val, double c_val);
double qlasso_fast_c_v1(double u, double a_val, double b_val, double c_val);
double qlasso_fast_c_v2(double u, double a_val, double b_val, double c_val);
double rlasso_fast_c_v1(double a_val, double b_val, double c_val);
double rlasso_fast_c_v2(double a_val, double b_val, double c_val);
double elasso_c_v1(double a_val, double b_val, double c_val);
double elasso_c_v2(double a_val, double b_val, double c_val);
double vlasso_c_v1(double a_val, double b_val, double c_val);
double vlasso_c_v2(double a_val, double b_val, double c_val);
double mlasso_c(double a_val, double b_val, double c_val);

////////////////////////////////////////////////////////////////////////////////

#endif
