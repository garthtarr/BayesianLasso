
#ifndef LASSO_DISTRIBUTION_H
#define LASSO_DISTRIBUTION_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////////////////////////////////////

double expit_c(double x);
//Rcpp::List calculate_lasso_dist_stats_c(double a_val, double b_val, double c_val);
//double logSumExp_c(vec vx);
double zlasso_c_v1(double a_val, double b_val, double c_val, bool logarithm);
double zlasso(double a_val, double b_val, double c_val, bool logarithm);
arma::vec dlasso_c_v1(arma::vec vx, double a_val, double b_val, double c_val, bool logarithm);
arma::vec dlasso(arma::vec vx, double a_val, double b_val, double c_val, bool logarithm);
arma::vec plasso_c_v1(arma::vec vx, double a_val, double b_val, double c_val);
arma::vec plasso(arma::vec vx, double a_val, double b_val, double c_val);
arma::vec qlasso_fast_c_v1(arma::vec u, double a_val, double b_val, double c_val);
arma::vec qlasso(arma::vec u, double a_val, double b_val, double c_val);
arma::vec rlasso_fast_c_v1(double n, double a_val, double b_val, double c_val);
arma::vec rlasso(double n, double a_val, double b_val, double c_val);
double elasso_c_v1(double a_val, double b_val, double c_val);
double elasso(double a_val, double b_val, double c_val);
double vlasso_c_v1(double a_val, double b_val, double c_val);
double vlasso(double a_val, double b_val, double c_val);
arma::vec mlasso(arma::vec a_val, arma::vec b_val, arma::vec c_val);

////////////////////////////////////////////////////////////////////////////////

#endif
