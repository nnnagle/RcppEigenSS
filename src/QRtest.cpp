#include "../inst/include/RcppEigenSS.h"
#include <Eigen/SPQRsupport>

using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
List QRtest( SEXP spmat ) {
  SparseMatrix<double> MAT = as<SparseMatrix<double> >(spmat);
  SPQR<SparseMatrix<double> > QR(MAT);
  // Create an identity matrix to get the Q matrix
  int size = QR.rows();
  SparseMatrix<double> eye(size,size);
  eye.setIdentity();
  MatrixXd eye2(size,size);
  eye2.setIdentity();
  SparseMatrix<double> R = QR.matrixR();
  //SparseMatrix<double> Q = QR.matrixQ() * eye;
  MatrixXd Q = QR.matrixQ() * eye2;
  MatrixXd P = QR.colsPermutation(); // This is an actual matrix
  List z  = List::create(Named("Q")=Rcpp::RcppEigenSS::eigen_wrap_plain_dense(Q),
                         Named("R")=Rcpp::RcppEigenSS::eigen_wrap_plain_sparse(R),
                         Named("P")=Rcpp::RcppEigenSS::eigen_wrap_plain_dense(P));
  return(z);
}
