#include <Rcpp.h>
#include "../inst/include/RcppEigenSS.h"
using namespace Rcpp;
using Eigen::MatrixXd;

// [[Rcpp::export]]
List rcpp_hello() {
  Eigen::MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  Rcpp::Rcout << m << std::endl << std::endl;

  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}
