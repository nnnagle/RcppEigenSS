// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigenSS.h: Rcpp/Eigen glue
//
// Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is part of MyEigen.
//
// RcppEigenSS is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigenSS is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MyEigen.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RcppEigenSS__RcppEigenSS__h
#define RcppEigenSS__RcppEigenSS__h

#include <iterator>
#include <Rcpp.h>
#include <RcppCommon.h>
#include <Rconfig.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/* forward declarations */
namespace Rcpp{
  /* support for wrap */


  namespace RcppEigenSS{
    // helper trait to identify if T is a plain object type
    // TODO: perhaps move this to its own file
    template <typename T> struct is_plain : Rcpp::traits::same_type<T,typename T::PlainObject>{} ;
    // helper trait to identify if the object has dense storage
    template <typename T> struct is_dense : Rcpp::traits::same_type<typename T::StorageKind,Eigen::Dense>{} ;
    // for plain dense objects
    template <typename T>
    SEXP eigen_wrap_plain_dense( const T& obj ){
      typename Eigen::internal::conditional<T::IsRowMajor,
                                            Eigen::Matrix<typename T::Scalar,
                                                          T::RowsAtCompileTime,
                                                          T::ColsAtCompileTime>,
                                                          const T&>::type objCopy(obj);
      int m = obj.rows(), n = obj.cols();
      SEXP ans = PROTECT(::Rcpp::wrap(objCopy.data(), objCopy.data() + m * n));
      if( T::ColsAtCompileTime != 1 ) {
        SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
        int *d = INTEGER(dd);
        d[0] = m;
        d[1] = n;
        ::Rf_setAttrib(ans, R_DimSymbol, dd);
        UNPROTECT(1);
      }
      UNPROTECT(1);
      return ans;
    }

    // for plain sparse objects
    template <typename T>
    SEXP eigen_wrap_plain_sparse( const T& object ){
      typedef typename T::Scalar     Scalar;
      const int  RTYPE = Rcpp::traits::r_sexptype_traits<Scalar>::rtype;
      std::string klass;
      switch(RTYPE) {
      case REALSXP: klass = T::IsRowMajor ? "dgRMatrix" : "dgCMatrix";
        break;
      default:
        throw std::invalid_argument("RTYPE not matched in conversion to sparse matrix");
      }
      Rcpp::S4           ans(klass);
      const int    nnz = object.nonZeros();
      ans.slot("Dim")  = Dimension(object.rows(), object.cols());
      ans.slot(T::IsRowMajor ? "j" : "i") =
        IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
      ans.slot("p")    = IntegerVector(object.outerIndexPtr(),
               object.outerIndexPtr() + object.outerSize() + 1);
      ans.slot("x")    = Vector<RTYPE>(object.valuePtr(), object.valuePtr() + nnz);
      return  ans;
    }
  }

  namespace traits {
    /* support for as */
  } // namespace traits
} // namespace Rcpp
#include <Rcpp.h>
#endif
