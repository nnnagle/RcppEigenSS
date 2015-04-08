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
    template <typename T, typename value_type>
    class MatrixExporterForEigen {
    public:
      typedef value_type r_export_type;

      MatrixExporterForEigen(SEXP x) : object(x){}
      ~MatrixExporterForEigen(){}

      T get() {
        Shield<SEXP> dims( ::Rf_getAttrib( object, R_DimSymbol ) );
        if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
          throw ::Rcpp::not_a_matrix();
        }
        int* dims_ = INTEGER(dims);
        T result( dims_[0], dims_[1] );
        value_type *data = result.data();
        ::Rcpp::internal::export_indexing<value_type*, value_type>( object, data );
        return result ;
      }

    private:
      SEXP object;
    };

    template <typename T>
    class Exporter<Eigen::Matrix<T, Eigen::Dynamic, 1> >
      : public IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T> {
    public:
      Exporter(SEXP x) : IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T >(x){}
    };

    template <typename T>
    class Exporter< Eigen::Matrix<T, 1, Eigen::Dynamic> >
      : public IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T > {
    public:
      Exporter(SEXP x) : IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T >(x){}
    };

    template <typename T>
    class Exporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
      : public MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T > {
    public:
      Exporter(SEXP x) :
      MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T >(x){}
    };

    template<typename T>
    class Exporter<Eigen::SparseMatrix<T> > {
    public:
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
      Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")), xx(d_x.slot("x")) {
        if (!d_x.is("dgCMatrix"))
          throw std::invalid_argument("Need S4 class dgCMatrix for a sparse matrix");
      }
      Eigen::SparseMatrix<T> get() {
        Eigen::SparseMatrix<T>  ans(d_dims[0], d_dims[1]);
        ans.reserve(d_p[d_dims[1]]);
        for(int j = 0; j < d_dims[1]; ++j) {
          ans.startVec(j);
          for (int k = d_p[j]; k < d_p[j + 1]; ++k) ans.insertBack(d_i[k], j) = xx[k];
        }
        ans.finalize();
        return ans;
      }
    protected:
      S4            d_x;
      IntegerVector d_dims, d_i, d_p;
      Vector<RTYPE> xx ;
    };

    template<typename T>
    class Exporter<Eigen::SparseMatrix<T, Eigen::RowMajor> > {
    public:
      const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;

      Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_j(d_x.slot("j")), d_p(d_x.slot("p")), xx(d_x.slot("x")) {
        if (!d_x.is("dgRMatrix"))
          throw std::invalid_argument("Need S4 class dgRMatrix for a sparse matrix");
      }
      Eigen::SparseMatrix<T, Eigen::RowMajor> get() {
        Eigen::SparseMatrix<T, Eigen::RowMajor>  ans(d_dims[0], d_dims[1]);
        ans.reserve(d_p[d_dims[0]]);
        for(int i = 0; i < d_dims[0]; ++i) {
          ans.startVec(i);
          for (int k = d_p[i]; k < d_p[i + 1]; ++k) ans.insertBack(i, d_j[k]) = xx[k];
        }
        ans.finalize();
        return ans;
      }
    protected:
      S4            d_x;
      IntegerVector d_dims, d_j, d_p;
      Vector<RTYPE> xx ;
    };

  } // namespace traits
} // namespace Rcpp
#include <Rcpp.h>
#endif
