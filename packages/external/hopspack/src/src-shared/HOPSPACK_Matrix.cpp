// $Id: HOPSPACK_Matrix.cpp 217 2013-11-25 21:59:49Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_Matrix.cpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2013 Sandia Corporation
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
  \file HOPSPACK_Matrix.cpp
  \brief Implement HOPSPACK::Matrix.
*/

#include <iomanip>
#include <math.h>                //-- FOR fabs

#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LapackWrappers.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_Print.hpp"

HOPSPACK::Matrix::Matrix() :
  fmatvecSet(false),
  fmatvecTSet(false)
{
}

HOPSPACK::Matrix::Matrix(const HOPSPACK::Matrix& source, TransposeType ttype) :
  fmatvecSet(false),
  fmatvecTSet(false)
{
  if (ttype == TRANSPOSE)
    transpose(source);
  else
    operator=(source);
}

HOPSPACK::Matrix::Matrix(const HOPSPACK::Matrix& source, const Vector& s, TransposeType ttype) :
  fmatvecSet(false),
  fmatvecTSet(false)
{
  if (ttype == TRANSPOSE)
    transpose(source);
  else
    operator=(source);
  scale(s);
}

HOPSPACK::Matrix::Matrix(double ** const Aptr, int nrows, int ncols) :
  fmatvecSet(false),
  fmatvecTSet(false)
{
  resize(nrows, ncols);
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      matrix[i][j] = Aptr[i][j];
}

HOPSPACK::Matrix &  HOPSPACK::Matrix::operator= (const HOPSPACK::Matrix &  m)
{
  matrix = m.matrix;
  matrixChanged();
  return( *this );
}

HOPSPACK::Matrix::~Matrix()
{
}

bool HOPSPACK::Matrix::empty() const
{
  return (matrix.size() == 0);
}

int HOPSPACK::Matrix::getNrows() const
{
    return( (int) matrix.size() );
}

int HOPSPACK::Matrix::getNcols() const
{
  if (matrix.empty())
    return 0;

  return matrix[0].size();
}

const HOPSPACK::Vector& HOPSPACK::Matrix::getRow(int i) const
{
  if ((i < 0) || (i >= getNrows()))
  {
    cerr << "ERROR: Matrix row " << i
         << " out of range  <HOPSPACK::Matrix.getRow()>" << endl;
    throw INTERNAL_ERROR;
  }

  return matrix[i];
}

void HOPSPACK::Matrix::getModifiableRowPointers(vector< double *>& Aptr)
{ 
  for (int i = 0; i < getNrows(); i++)
    Aptr.push_back(&matrix[i][0]);
  matrixChanged();
}

void HOPSPACK::Matrix::clear()
{
  resize(0,0);
  matrixChanged();;
}

void HOPSPACK::Matrix::addRow(const HOPSPACK::Vector& r)
{
  if ((!matrix.empty()) && (r.size() != getNcols()))
  {
    cerr << "ERROR: Matrix row size mismatch  <HOPSPACK::Matrix.addRow()>"
         << endl;
    throw INTERNAL_ERROR;
  }
  
  matrix.push_back(r);
  matrixChanged();
}

void HOPSPACK::Matrix::deleteRow(int i)
{
  if ((i < 0) || (i >= getNrows()))
  {
    cerr << "ERROR: Matrix row " << i
         << " out of range  <HOPSPACK::Matrix.deleteRow()>" << endl;
    throw INTERNAL_ERROR;
  }

  matrix.erase(matrix.begin()+i);
  matrixChanged();
}

void HOPSPACK::Matrix::addRow(const HOPSPACK::Vector& r, double alpha)
{
  addRow(r);
  matrix[matrix.size()-1].scale(alpha);
  matrixChanged();
}

void HOPSPACK::Matrix::addMatrix(const HOPSPACK::Matrix& B)
{
  for (int i = 0; i < B.getNrows(); i ++)
    addRow(B.getRow(i));
  matrixChanged();
}

void HOPSPACK::Matrix::addMatrix(const HOPSPACK::Matrix& B, double alpha)
{
  for (int i = 0; i < B.getNrows(); i ++)
    addRow(B.getRow(i), alpha);
  matrixChanged();
}

void HOPSPACK::Matrix::addMatrix(const HOPSPACK::Matrix& B, const HOPSPACK::Vector& s)
{
  for (int i = 0; i < B.getNrows(); i ++)
  {
    addRow(B.getRow(i));
    matrix[matrix.size()-1].scale(s);
  }
  matrixChanged();
}

void HOPSPACK::Matrix::addUniqueRows(const HOPSPACK::Matrix& B, double epsilon)
{
  // Record number of rows in A now, once we are appending we don't want
  // to end up looping over the columns of B also.
  int nArows = getNrows();
  int nBrows = B.getNrows();
  Vector diff(getNcols());
  for (int i = 0; i < nBrows; i++)
  {
    const Vector& bi = B.getRow(i);
    bool bi_unique = true;
    for (int j = 0; j < nArows; j++)
    {
      diff = getRow(j);          // diff = aj.
      diff -= bi;                // diff = aj - bi.
      if (diff.norm() < epsilon) // Is ||aj - bi|| < epsilon?
      {
        // row is already in A.
        bi_unique = false;
        break;
      }
    }
    if (bi_unique)
      addRow(bi);
  }
  
  matrixChanged();
}

void HOPSPACK::Matrix::transpose(const Matrix& source)
{
  int nrows = source.getNcols();
  int ncols = source.getNrows();

  resize(nrows,ncols);

  for (int i = 0; i < nrows; i ++)
    for (int j = 0; j < ncols; j ++)
      matrix[i][j] = source.matrix[j][i];
  matrixChanged();
}

void HOPSPACK::Matrix::scale(const HOPSPACK::Matrix& source, const Vector& scaling)
{
  // Copy source
  matrix = source.matrix; 
  
  // Scale columns
  scale(scaling);
  matrixChanged();
}

void HOPSPACK::Matrix::scale(const Vector& scaling)
{
  // Scale columns
  int nrows = getNrows();
  int ncols = getNcols();
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      matrix[i][j] = matrix[i][j] * scaling[j];
  matrixChanged();
}

void HOPSPACK::Matrix::setToIdentity(int n)
{
  resize(n,n);
  for (int i = 0; i < n; i++)
  {
    matrix[i].zero();
    matrix[i][i] = 1.0;
  }
  matrixChanged();
}

void HOPSPACK::Matrix::normalize()
{
  for (int i = 0; i < getNrows(); i++)
  {
    double cnorm = matrix[i].norm();
    if (cnorm == 0)
      deleteRow(i);
    else
      matrix[i].scale(1/cnorm);
  }
  matrixChanged();
}

void  HOPSPACK::Matrix::multVec(const HOPSPACK::Vector& x, HOPSPACK::Vector& y,
                                TransposeType ttype) const
{
  // Error check.
  if (ttype == TRANSPOSE)
  {
    if (x.size() != getNrows())
    {
      cerr << "ERROR: Matrix size mismatch with input vector x"
           << "  <HOPSPACK::Matrix.multVec()>" << endl;
      throw INTERNAL_ERROR;
    }
    
    if (y.size() != getNcols())
    {
      cerr << "ERROR: Matrix size mismatch with input vector y"
           << "  <HOPSPACK::Matrix.multVec()>" << endl;
      throw INTERNAL_ERROR;
    }
  }
  else
  {
    if (x.size() != getNcols())
    {
      cerr << "ERROR: Matrix size mismatch with input vector x"
           << "  <HOPSPACK::Matrix.multVec()>" << endl;
      throw INTERNAL_ERROR;
    }
    
    if (y.size() != getNrows())
    {
      cerr << "ERROR: Matrix size mismatch with input vector y"
           << "  <HOPSPACK::Matrix.multVec()>" << endl;
      throw INTERNAL_ERROR;
    }
  }

  #if defined(HAVE_LAPACK)
    multVecWithBlas(x,y,ttype);
  #else
    multVecWithoutBlas(x,y,ttype);
  #endif

  return;
}

void HOPSPACK::Matrix::svd(HOPSPACK::Matrix& U,
                           HOPSPACK::Vector &s,
                           HOPSPACK::Matrix& VT) const 
{
  int m = getNrows();
  int n = getNcols();

  // Convert matrix to a FORTRAN vector
  Vector Avec = getMatrixVector();

  // Set up vectors to be filled by the SVD routine

  s.resize(m);          // Store diag(Sigma)
  Vector Uvec(m*m);     // Store U
  Vector VTvec(n*n);    // Store V^T

  bool  bSuccess
      = LapackWrappers::getTheInstance().dgesvd ('A', 'A',
                                                 m, n,
                                                 &(Avec[0]),
                                                 &(s[0]),
                                                 &(Uvec[0]), m,
                                                 &(VTvec[0]), n);
  if (bSuccess == false)
  {
    cerr << "ERROR: Call to LAPACK function dgesvd failed" << endl;
    throw INTERNAL_ERROR;
  }

  U.copyFromFortranVector(Uvec,m,m);
  VT.copyFromFortranVector(VTvec,n,n);
  return;
}

void HOPSPACK::Matrix::multMat(const HOPSPACK::Matrix& B,
                               TransposeType ttype)
{
  Matrix C;
  multMat(B,C,ttype);
  matrix = C.matrix;
  matrixChanged();
}

void HOPSPACK::Matrix::multMat(const HOPSPACK::Matrix& B, HOPSPACK::Matrix& C,
                               TransposeType ttype) const
{
  if (ttype == TRANSPOSE)
  {
    // Peforming mult of form (m x k)*(n x k)' = (m x n).
    if (getNcols() != B.getNcols())
    {
      cerr << "ERROR: Matrix has wrong number of columns"
           << "  <HOPSPACK::Matrix.multMat()>" << endl;
      throw INTERNAL_ERROR;
    }
  }
  else
  {
    // Peforming mult of form (m x k)*(k x n) = (m x n).
    if (getNcols() != B.getNrows())
    {
      cerr << "ERROR: Matrix has wrong number of rows"
           << "  <HOPSPACK::Matrix.multMat()>" << endl;
      throw INTERNAL_ERROR;
    }
  }

  #if defined(HAVE_LAPACK)
    multMatWithBlas(B,C,ttype);
  #else
    multMatWithoutBlas(B,C,ttype);
  #endif

  return;
}

void HOPSPACK::Matrix::nullSpace(HOPSPACK::Matrix& ZT, double tol) const
{
  // Determine Z such that the columns of Z span the nullspace of this matrix. 

  int m = getNrows();
  int n = getNcols();

  if ((m == 0) || (n == 0))
  {
    cerr << "ERROR: Input matrix is empty"
         << "  <HOPSPACK::Matrix.nullSpace()>" << endl;
    throw INTERNAL_ERROR;
  }

  Vector s;
  Matrix U;
  Matrix VT;

  // Get SVD decomposition
  svd(U,s,VT);

  // Determine rank
  int rank = s.size();
  for (int i = 0; i < s.size(); i++)
    if (s[i] < tol)
    {
      rank = i;
      break;
    }

  // Form Z' by adding rows of V' that correspond to the zero singular values.
  ZT.copySubMatrix(rank,n-rank, VT);
}

bool HOPSPACK::Matrix::getRightInvAndNullBasis(HOPSPACK::Matrix& RT, HOPSPACK::Matrix& NT,
                                               double tol) const
{
  RT.clear();
  NT.clear();

  int m = getNrows();
  int n = getNcols();

  if (m > n)
    return false; // Pseudo-inverse does not exist.

  // First compute the singular value decomposition of A.
  Matrix U;
  Vector s;
  Matrix VT;
  svd(U,s,VT);
  
  for (int i = 0; i < s.size(); i++)
    if (s[i] < tol)
      return false; // Right inverse does not exists.
  
  // First, let's partition V to be V = [Vr N], i.e., the columns that
  // span the range and nullspaces.
  Matrix VrT;
  VrT.copySubMatrix(0, m, VT);
  NT.copySubMatrix(m, n-m, VT);

  // We know that R = Vr*inv(S)*UT => RT = U*inv(S)*VrT.
  for (int i = 0; i < s.size(); i++)
    s[i] = 1 / s[i];
  
  RT = U;
  RT.scale(s);
  RT.multMat(VrT);

  return true;
}

void HOPSPACK::Matrix::pruneDependentRows(HOPSPACK::Vector& b, double epsilon)
{
  int m = getNrows();
  int n = getNcols();

  // Convert matrix to a FORTRAN vector
  Vector Avec = getMatrixVector();

  // Set up tau to be filled by the LQ routine
  Vector tau (m, 0.0);

  bool  bSuccess
      = LapackWrappers::getTheInstance().dgelqf (m, n,
                                                 &(Avec[0]),
                                                 &(tau[0]));
  if (bSuccess == false)
  {
    cerr << "ERROR: Call to LAPACK function dgelqf failed" << endl;
    throw INTERNAL_ERROR;
  }

  // Delete rows of A and entries of b if necesary.
  // Iterate backwards as matrix size may change, i.e. if we were to iterate
  // forward, i for Avec will be different than the i for A and b.
  for (int i = m-1; i >= 0; i--)
  {
    if (fabs(Avec[i+i*m]) < epsilon)
    {
      deleteRow(i);
      b.erase(i);
      matrixChanged();
    }
  }

  return;
}


bool  HOPSPACK::Matrix::specialConstrainedLSQR
          (      HOPSPACK::Vector &  cX,
           const HOPSPACK::Vector &  cB) const
{  
    if (empty())
        return( true );

    //---- IN dgglse TERMINOLOGY, MATRIX "A" IS THE IDENTITY, MATRIX "B" IS this,
    //---- VECTOR "c" IS THE POINT cX, AND VECTOR "d" IS cB.

    int  nRows = getNrows();
    int  nCols = getNcols();
    if (nRows > nCols)
    {
        cerr << "ERROR: Cannot solve least squares overdetermined system"
             << endl;
        cerr << "       num rows = " << nRows << " is > num cols = " << nCols
             << endl;
        return( false );
    }

    Vector Avec = getMatrixVector();
    Vector bCopy = cB;
    Vector xSolution (nCols);

    Vector Identity (nCols*nCols, 0.0);
    for (int  i = 0; i < nCols; i++)
        Identity[i + i*nCols] = 1.0;

    bool  bSuccess
        = LapackWrappers::getTheInstance().dgglse (nCols, nCols, nRows,
                                                   &(Identity[0]),
                                                   &(Avec[0]),
                                                   &(cX[0]),
                                                   &(bCopy[0]),
                                                   &(xSolution[0]));
    if (bSuccess == false)
    {
        cerr << "ERROR: Call to LAPACK function dgglse failed" << endl;
        return( false );
    }

    cX = xSolution;
    return( true );
}


bool  HOPSPACK::Matrix::generalConstrainedLSQR
          (const HOPSPACK::Vector &  cC,
           const HOPSPACK::Vector &  cD,
           const HOPSPACK::Vector &  cB,
                 HOPSPACK::Vector &  cX) const
{  
    if (empty())
        //---- NOTE THIS RETURNS cX AS THE SOLUTION; OTHERWISE, cX IS IGNORED.
        return( true );

    //---- IN dgglse TERMINOLOGY, MATRIX "A" HAS DIAGONAL cD, MATRIX "B" IS this,
    //---- VECTOR "c" (OBJ) IS THE POINT cC, AND VECTOR "d" (RHS) IS cB.

    int  nRows = getNrows();
    int  nCols = getNcols();
    if (nRows > nCols)
    {
        cerr << "ERROR: Cannot solve least squares overdetermined system"
             << endl;
        cerr << "       num rows = " << nRows << " is > num cols = " << nCols
             << endl;
        return( false );
    }

    Vector Avec = getMatrixVector();
    Vector cCopy = cC;
    Vector bCopy = cB;

    Vector Dvec (nCols*nCols, 0.0);
    for (int  i = 0; i < nCols; i++)
        Dvec[i + i*nCols] = cD[i];

    bool  bSuccess
        = LapackWrappers::getTheInstance().dgglse (nCols, nCols, nRows,
                                                   &(Dvec[0]),
                                                   &(Avec[0]),
                                                   &(cCopy[0]),
                                                   &(bCopy[0]),
                                                   &(cX[0]));
    if (bSuccess == false)
    {
        cerr << "ERROR: Call to LAPACK function dgglse failed" << endl;
        return( false );
    }

    return( true );
}


bool  HOPSPACK::Matrix::generalLS (const HOPSPACK::Vector &  cC,
                                         HOPSPACK::Vector &  cX) const
{
    if (empty())
        return( false );

    int  nRows = getNrows();
    int  nCols = getNcols();
    Vector Avec = getMatrixVector();
    Vector cCopy = cC;

    bool  bSuccess
        = LapackWrappers::getTheInstance().dgelss (nRows, nCols,
                                                   &(Avec[0]),
                                                   &(cCopy[0]),
                                                   &(cX[0]));
    if (bSuccess == false)
    {
        cerr << "ERROR: Call to LAPACK function dgelss failed" << endl;
        return( false );
    }

    return( true );
}


void HOPSPACK::Matrix::formattedPrint (const string  &  sPrefix,
                                             ostream &  stream) const
{
  int  nPrec = Print::getPrecision();
  stream.setf(ios::scientific);
  stream.precision(nPrec);

  stream << endl;
  stream << sPrefix << "[" << endl;
  for (int i = 0; i < getNrows(); i++) 
  {
    stream << sPrefix;
    for (int j = 0; j < getNcols(); j++) 
    {
      double  d = matrix[i][j];
      if (exists(d))
        #if defined(WIN32)
          //---- WINDOWS USES FORMAT SX.XXXeSXXX.
          stream << setw(nPrec + 8) << d << " ";
        #else
          //---- UNIX USES FORMAT    SX.XXXeSXX.
          stream << setw(nPrec + 7) << d << " ";
        #endif
      else
      {
        stream << " DNE";
        #if defined(WIN32)
          for (int k = 0; k < nPrec + 5; k++)
            stream << " ";
        #else
          for (int k = 0; k < nPrec + 4; k++)
            stream << " ";
        #endif
      }
    }
    stream << endl;
  }
  stream << sPrefix << "]";

  stream.unsetf(ios::scientific);
  stream.flush();
  return;
}

void HOPSPACK::Matrix::matrixChanged()
{
  fmatvecTSet = false;
  fmatvecSet = false;
}

const HOPSPACK::Vector& HOPSPACK::Matrix::getMatrixVector(TransposeType ttype) const
{
  if (ttype == TRANSPOSE)
  {  
    if (!fmatvecTSet) // Vectorized matrix not current.  Update.
    {
      copyToFortranVector(fmatvecT, ttype);
      fmatvecTSet = true;
    }
    return fmatvecT;
  }
  else
  {
    if (!fmatvecSet) // Vectorized matrix not current.  Update.
    {
      copyToFortranVector(fmatvec, ttype);
      fmatvecSet = true;
    }
    return fmatvec;
  }
}

void HOPSPACK::Matrix::copyToFortranVector(HOPSPACK::Vector& Avec, TransposeType ttype) const
{
  // Fortran stores its matrices stored columnwise.
  int nrows = getNrows();
  int ncols = getNcols();

  Avec.resize(0);
  Avec.reserve(nrows * ncols);
  
  if (ttype == TRANSPOSE)
  {  
    // Store transpose of A. Elements can be added an entire row at a time.
    for (int i = 0; i < nrows; i++)
      Avec.append(matrix[i]);
  }
  else
  {
    // Store A. Elements must be added individually.
    for (int j = 0; j < ncols; j++)
      for (int i = 0; i < nrows; i ++)
        Avec.push_back(matrix[i][j]);
  }
  
}

void HOPSPACK::Matrix::copyFromFortranVector(const HOPSPACK::Vector& Avec,
                                             int nrows, int ncols,
                                             TransposeType ttype)
{
  // Fortran stores its matrices stored columnwise.

  resize(nrows,ncols);
  
  if (ttype == TRANSPOSE)
  {
    // Avec contains a transposed version of A. Unpack row-by-row.
    int k = 0;
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < ncols; j++)
        matrix[i][j] = Avec[k++];
  }
  else
  {
    // Avec contains A. Unpack column-by-column.
    int k = 0;
    for (int j = 0; j < ncols; j++)
      for (int i = 0; i < nrows; i++)
        matrix[i][j] = Avec[k++];
  }
  
  matrixChanged();
}

// Private
void HOPSPACK::Matrix::resize(int nrows, int ncols)
{
  matrix.resize(nrows);
  for (int i = 0; i < nrows; i++)
    matrix[i].resize(ncols);
  matrixChanged();
}

// Private
void HOPSPACK::Matrix::copySubMatrix(int istart, int iextent,
                                     const HOPSPACK::Matrix& B)
{
  if ((istart + iextent) > B.getNrows()) 
  {
    cerr << "ERROR: Bad submatrix size "
         << (istart + iextent) << " vs " << B.getNrows()
         << "  <HOPSPACK::Matrix.copySubMatrix()>" << endl;
    throw INTERNAL_ERROR;
  }

  clear();
  for (int i = 0; i < iextent; i++)
    addRow(B.getRow(i+istart));
  matrixChanged();
}


//......................................................................
#if defined(HAVE_LAPACK)

// Private
void  HOPSPACK::Matrix::multVecWithBlas(const HOPSPACK::Vector& x,
                                        HOPSPACK::Vector& y,
                                        TransposeType ttype) const
{
  char trans;
  if (ttype == TRANSPOSE)
    trans='T'; // Performing tranpose mult.
  else
    trans='N'; // Performing normal mult.

  Vector &  Acopy = const_cast< Vector & >(getMatrixVector());
  Vector &  xcopy = const_cast< Vector & >(x);
  LapackWrappers::getTheInstance().dgemv (trans,
                                          getNrows(),
                                          getNcols(),
                                          1.0,
                                          &(Acopy[0]),
                                          &(xcopy[0]),
                                          0.0,
                                          &(y[0]));
  return;
}

#else     //-- NOT HAVE_LAPACK

// Private
void  HOPSPACK::Matrix::multVecWithoutBlas(const HOPSPACK::Vector& x,
                                           HOPSPACK::Vector& y,
                                           TransposeType ttype) const
{
  if (ttype == TRANSPOSE)
  {
    // Compute y = A' * x, i.e. y_j = sum_i A_{ij} x_i
    y.zero();
    for (int i = 0; i < getNrows(); i++)
      for (int j = 0; j < getNcols(); j++)
        y[j] += matrix[i][j] * x[i];
  }
  else
  {
    // Compute y = A * x, i.e., y_i = sum_j A_{ij} x_j
    y.zero();
    for (int i = 0; i < getNrows(); i++)
      for (int j = 0; j < getNcols(); j++)
        y[i] += matrix[i][j] * x[j];
  } 
}

#endif     //-- HAVE_LAPACK
//......................................................................


//......................................................................
#if defined(HAVE_LAPACK)

// Private
void HOPSPACK::Matrix::multMatWithBlas(const HOPSPACK::Matrix& B,
                                             HOPSPACK::Matrix& C,
                                       const TransposeType ttype) const
{  
  int  m = getNrows();
  int  k = getNcols();

  int  n;      // n denotes number of columns of C.
  char trans;
  if (ttype == TRANSPOSE)
  {
    // AS EXPLAINED BELOW, PERFORM MULTIPLY OF FORM (m x k)*(n x k)' = (m x n).
    trans = 'N';
    n     = B.getNrows();
  }
  else
  {
    // AS EXPLAINED BELOW, PERFORM MULTIPLY OF FORM (m x k)*(k x n) = (m x n).
    trans = 'T';
    n     = B.getNcols();
  }

  //---- ACTUALLY COMPUTE C = (A')' * (B')' BECAUSE IT'S MORE EFFICIENT TO
  //---- CONVERT A' AND B' TO VECTORS.
  Vector &  ATcopy = const_cast< Vector & >(getMatrixVector (TRANSPOSE));
  Vector &  BTcopy = const_cast< Vector & >(B.getMatrixVector (TRANSPOSE));

  Vector  Cvec (m * n);
  LapackWrappers::getTheInstance().dgemm ('T', trans,
                                          m, n, k,
                                          1.0,
                                          &(ATcopy[0]),
                                          &(BTcopy[0]),
                                          0.0,
                                          &(Cvec[0]));
  C.copyFromFortranVector (Cvec, m, n);
  return;
}

#else     //-- NOT HAVE_LAPACK

// Private
void HOPSPACK::Matrix::multMatWithoutBlas(const HOPSPACK::Matrix& B,
                                          HOPSPACK::Matrix& C,
                                          TransposeType ttype) const
{
  int m = getNrows();
  int k = getNcols();
  if (ttype == TRANSPOSE)
  {
    // Peforming mult of form (m x k)*(n x k)' = (m x n).
    int n = B.getNrows();
    C.resize(m,n);
    // C(i,j) = sum_{l=1}^k A(i,l)*B(j,l)
    for( int i = 0; i < m; i++)
      for( int j = 0; j < n; j++)
      {
        C.matrix[i][j] = 0;
        for( int l = 0; l < k; l++)
          C.matrix[i][j] += matrix[i][l] * B.matrix[j][l];
      }
  }
  else
  {
    // Peforming mult of form (m x k)*(k x n) = (m x n).
    int n = B.getNcols();
    C.resize(m,n);
    // C(i,j) = sum_{l=1}^k A(i,l)*B(l,j)
    for( int i = 0; i < m; i++)
      for( int j = 0; j < n; j++)
      {
        C.matrix[i][j] = 0;
        for( int l = 0; l < k; l++)
          C.matrix[i][j] += matrix[i][l] * B.matrix[l][j];
      }
  }
}

#endif     //-- HAVE_LAPACK
//......................................................................
