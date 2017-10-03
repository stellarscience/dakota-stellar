// $Id: HOPSPACK_LapackWrappers.hpp 208 2012-08-01 22:59:33Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_LapackWrappers.hpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This file is part of HOPSPACK.
//
// HOPSPACK is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see http://www.gnu.org/licenses/.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
//                 or Todd Plantenga (tplante@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*!
  @file HOPSPACK_LapackWrappers.hpp
  @brief Declaration for LAPACK wrapper methods used by HOPSPACK.
*/
#ifndef HOPSPACK_LAPACKWRAPPERS_HPP
#define HOPSPACK_LAPACKWRAPPERS_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Provides a uniform interface to selected LAPACK functions.
/*!
 *  The class provides a static singleton instance that offers selected
 *  LAPACK functions in a simplified, type-safe form.
 *
 *  The implementation coordinates with CMake to support specific calling
 *  conventions used by various LAPACK providers.  To add a new convention,
 *  look at src/ConfigureLinConstr.cmake and HOPSPACK_LapackWrappers.cpp.
 */
//----------------------------------------------------------------------
class LapackWrappers
{
  public:

    //! Applications call this instead of a constructor.
    static LapackWrappers &  getTheInstance (void);

    //! Applications should not call this directly (managed internally).
    ~LapackWrappers (void);


    //! Return the dot product of two vectors.
    /*!
     *  In most libraries this is a BLAS routine.
     *
     *  @param[in] n   Length of vectors x and y.
     *  @param[in] x   Input vector.
     *  @param[in] y   Input vector.
     *  @return        Dot product of dx and dy.
     */
    double ddot (const int             n,
                 const double * const  x,
                 const double * const  y) const;


    //! Replace y with y = alpha*Ax + beta*y (or A-transpose).
    /*!
     *  Compute
     *  \f[
     *    y := \alpha A x + \beta y, \;{\rm or }\;  y := \alpha A^T x + \beta y
     *  \f]
     *  where \f$ \alpha \f$ and \f$ \beta \f$ are scalars,
     *  \f$ x \f$ and \f$ y \f$ are vectors,
     *  and \f$ A \f$ is an \f$ m \times n \f$ matrix.
     *
     *  The dense matrix is assumed to be in column-major order (all elements
     *  of column 1, then all elements of column 2, etc.).
     *
     *  @param[in]     trans  'N' for normal, 'T' means take the transpose of A.
     *  @param[in]     m      Number of rows in A.
     *  @param[in]     n      Number of columns in A.
     *  @param[in]     alpha  Scalar value.
     *  @param[in]     A      Dense matrix in column-major order.
     *  @param[in]     x      Dense vector of length n (or length m if trans='T').
     *  @param[in]     beta   Scalar value.
     *  @param[in,out] y      Dense vector of length m (or length n if trans='T').
     */
    void  dgemv (const char            trans,
                 const int             m,
                 const int             n,
                 const double          alpha,
                 const double * const  A,
                 const double * const  x,
                 const double          beta,
                       double * const  y) const;


    //! Replace C with C = alpha*(A*B) + beta*C (or the transpose of A or B).
    /*!
     *  Compute
     *  \f[
     *    C := \alpha * op( A ) * op( B ) + \beta* C
     *  \f]
     *  where  \f$ op( X ) \f$
     *  is one of \f$ op( X ) = X \f$ or \f$ op( X ) = X^T \f$,
     *  and \f$ \alpha \f$ and \f$ \beta \f$ are scalars.
     *  \f$ A \f$, \f$ B \f$, and \f$C\f$ are matrices,
     *  with \f$ op(A) \f$ an \f$ m \times k \f$ matrix,
     *  \f$ op(B) \f$  a \f$ k \times n \f$ matrix,
     *  and \f$ C \f$ an \f$m \times n \f$ matrix.
     *
     *  The dense matrices are assumed to be in column-major order (all elements
     *  of column 1, then all elements of column 2, etc.).
     *
     *  @param[in]     transA 'N' for normal, 'T' means take the transpose of A.
     *  @param[in]     transB 'N' for normal, 'T' means take the transpose of B.
     *  @param[in]     m      Number of rows in op(A) and C.
     *  @param[in]     n      Number of columns in op(B) and C.
     *  @param[in]     k      Number of columns in op(A) and rows in op(B).
     *  @param[in]     alpha  Scalar value.
     *  @param[in]     A      Dense matrix, m by k.
     *  @param[in]     B      Dense matrix, k by n.
     *  @param[in]     beta   Scalar value.
     *  @param[in,out] C      Dense matrix, m by n.
     */
    void  dgemm (const char            transA,
                 const char            transB,
                 const int             m,
                 const int             n,
                 const int             k,
                 const double          alpha,
                 const double * const  A,
                 const double * const  B,
                 const double          beta,
                       double * const  C) const;

    //! Compute the SVD of rectangular matrix A.
    /*!
     *  Given \f$ m \times n \f$ matrix \f$ A \f$, compute
     *  \f[
     *    A := U \Sigma V^T
     *  \f]
     *  where \f$ \Sigma \f$ is an \f$ m \times n \f$ matrix which is zero
     *  except for its min(m,n) diagonal elements,
     *  \f$ U \f$ is an \f$ m \times m \f$ orthogonal matrix, and
     *  \f$ V \f$ is an \f$ n \times n \f$ orthogonal matrix
     *.  The diagonal elements of \f$ \Sigma \f$ are returned with the
     *  singular values of \f$ A \f$; they are real and nonnegative, and
     *  are returned in descending order.  The first min(m,n) columns of
     *  \f$ U \f$ and \f$ V \f$ are the left and right singular vectors of
     *  \f$ A \f$.
     *  Note that the routine returns \f$ V^T \f$, not \f$ V \f$.
     *
     *  @param[in]     jobU   'A' means return all \f$m\f$ cols of \f$U\f$ in U,
     *                        'S' means return first min(m,n) cols in U,
     *                        'O' means return first min(m,n) cols in A,
     *                        'N' means no left singular vectors are returned.
     *  @param[in]     jobVT  'A' means return all \f$n\f$ rows of \f$V^T\f$ in VT,
     *                        'S' means return first min(m,n) rows in VT,
     *                        'O' means return first min(m,n) rows in A,
     *                        'N' means no righft singular vectors are returned.
     *  @param[in]     m       Number of rows in A.
     *  @param[in]     n       Number of columns in A.
     *  @param[in,out] A       Dense matrix, m by n, overwritten always.
     *  @param[out]    Sigma   Dense vector of length min(m,n) for singular values.
     *  @param[out]    U       Dense matrix, typically m by m.
     *  @param[in]     ldU     Leading dimension of U, typically m.
     *  @param[out]    VT      Dense matrix, typically n by n.
     *  @param[in]     ldVT    Leading dimension of VT, typically n.
     *  @return True if OK, false if there was an error.
     */
    bool  dgesvd (const char            jobU,
                  const char            jobVT,
                  const int             m,
                  const int             n,
                        double * const  A,
                        double * const  Sigma,
                        double * const  U,
                  const int             ldU,
                        double * const  VT,
                  const int             ldVT) const;


    //! Solve a linear least squares problem with equality constraints.
    /*!
     *  Given \f$ m \times n \f$ matrix \f$ A \f$, solve
     *  \f[
     *    \mbox{minimize} \; \| c - Ax \|_2 \quad
     *    \mbox{subject to} \;\; Bx = d,
     *  \f]
     *  where \f$ B \f$ is a \f$ p \times n \f$ matrix,
     *  \f$ c \f$ is a vector of length \f$ m \f$,
     *  and \f$ d \f$ is a vector of length \f$ p \f$.
     *  To ensure the problem has a unique solution, assume that
     *  \f[
     *    p <= n <= m+p, \quad
     *    \mbox{rank}(B) = p, \quad
     *    \mbox{and} \;\;
     *    \mbox{rank} \left( \begin{array}{c} A \\
     *                                        B \end{array} \right)
     *    = n.
     *  \f]
     *
     *  @param[in]     m  Number of rows in A.
     *  @param[in]     n  Number of columns in A.
     *  @param[in]     p  Number of rows in B.
     *  @param[in,out] A  Dense matrix, m by n, overwritten.
     *  @param[in,out] B  Dense matrix, p by n, overwritten.
     *  @param[in,out] c  Dense vector of length m, returns the residual sum
     *                    of squares of elements n-p+1 to m of vector c.
     *  @param[in,out] d  Dense vector of length p, overwritten.
     *  @param[in,out] x  Dense vector of length n, returns the solution.
     *  @return True if OK, false if there was an error.
     */
    bool  dgglse (const int             m,
                  const int             n,
                  const int             p,
                        double * const  A,
                        double * const  B,
                        double * const  c,
                        double * const  d,
                        double * const  x) const;


    //! Solve an overdetermined linear least squares problem without constraints.
    /*!
     *  Given \f$ m \times n \f$ matrix \f$ A \f$ with \f$ m \geq n \f$, solve
     *  \f[
     *    \mbox{minimize} \; \| c - Ax \|_2 \quad
     *  \f]
     *  where \f$ c \f$ is a vector of length \f$ m \f$.
     *  To ensure the problem has a unique solution, assume \f$ A \f$
     *  has full rank.
     *
     *  @param[in]     m  Number of rows in A.
     *  @param[in]     n  Number of columns in A.
     *  @param[in,out] A  Dense matrix, m by n, overwritten.
     *  @param[in]     c  Dense vector of length m.
     *  @param[out]    x  Dense vector of length n, returns the solution.
     *  @return True if OK, false if there was an error.
     */
    bool  dgelss (const int             m,
                  const int             n,
                        double * const  A,
                  const double * const  c,
                        double * const  x) const;


    //! Compute the LQ factorization of rectangular matrix A.
    /*!
     *  Given \f$ m \times n \f$ matrix \f$ A \f$, compute
     *  \f[
     *    A = L * Q
     *  \f]
     *  where \f$ L \f$ is lower triangular and \f$ Q \f$ is orthogonal.
     *
     *  @param[in]     m  Number of rows in A.
     *  @param[in]     n  Number of columns in A.
     *  @param[in,out] A  Dense matrix, m by n, overwritten with L.
     *  @param[out]    T  Dense matrix, m by n, returns scalar factors for Q.
     *  @return True if OK, false if there was an error.
     */
    bool  dgelqf (const int             m,
                  const int             n,
                        double * const  A,
                        double * const  T) const;


  private:

    //! Singleton implementation hides the constructor.
    LapackWrappers (void);
    //! By design, there is no copy constructor.
    LapackWrappers (const LapackWrappers &);
    //! By design, there is no assignment operator.
    LapackWrappers & operator= (const LapackWrappers &);

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_LAPACKWRAPPERS_HPP
