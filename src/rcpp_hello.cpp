#include "../inst/include/RcppEigenSS.h"

using namespace Rcpp;
using namespace Eigen;
//using Eigen::MatrixXd;

// [[Rcpp::export]]
List rcpp_hello(SEXP mat) {
  MatrixXd MAT = as<MatrixXd>(mat);
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  Rcpp::Rcout << m << std::endl << std::endl;
  SparseMatrix<double> Spm = m.sparseView();
  Rcpp::Rcout << Spm << std::endl << std::endl;

  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y,
                                   Rcpp::RcppEigenSS::eigen_wrap_plain_dense(MAT),
                                   Rcpp::RcppEigenSS::eigen_wrap_plain_sparse(Spm));
  return z;
}
