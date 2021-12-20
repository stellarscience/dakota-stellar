//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2017 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H
#define UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H

#include <queso/VectorRV.h>

namespace QUESO {

class GslVector;
class GslMatrix;
template <class V, class M> class VectorSet;

/*!
 * \class InvLogitGaussianVectorRV
 * \brief A class representing a (transformed) Gaussian vector RV with bounds
 *
 * This class allows the user to compute the value of a (transoformed) Gaussian
 * PDF and to generate realizations (samples) from it.
 */

template <class V = GslVector, class M = GslMatrix>
class InvLogitGaussianVectorRV : public BaseVectorRV<V, M> {
public:
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*!
   * Construct a (transformed) Gaussian vector RV with mean \c lawExpVector
   * (of the Gaussian, not the transformed Gaussian) and diagonal covariance
   * matrix \c lawVarVector whose variates live in \c imageSet.
   */
  InvLogitGaussianVectorRV(const char * prefix,
      const VectorSet<V, M> & imageSet, const V & lawExpVector,
      const V & lawVarVector);

  //! Constructor
  /*!
   * Construct a (transformed) Gaussian vector RV with mean \c lawExpVector
   * (of the Gaussian, not the transformed Gaussian) and covariance matrix
   * \c lawCovMatrix whose variates live in \c imageSet.
   */
  InvLogitGaussianVectorRV(const char * prefix,
      const VectorSet<V, M> & imageSet, const V & lawExpVector,
      const M & lawCovMatrix);

  //! Virtual destructor
  virtual ~InvLogitGaussianVectorRV();
  //@}

  //! @name Statistical methods
  //@{
  //! Updates the vector that contains the mean values for the underlying Gaussian.
  void updateLawExpVector(const V & newLawExpVector);

  //! Updates the covariance matrix.
  /*!
   * This method tries to use Cholesky decomposition; and if it fails, the
   * method then calls a SVD decomposition.
   */
  void updateLawCovMatrix(const M & newLawCovMatrix);
  //@}

  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream & os) const;
 //@}

private:
  using BaseVectorRV<V,M>::m_env;
  using BaseVectorRV<V,M>::m_prefix;
  using BaseVectorRV<V,M>::m_imageSet;
  using BaseVectorRV<V,M>::m_pdf;
  using BaseVectorRV<V,M>::m_realizer;
  using BaseVectorRV<V,M>::m_subCdf;
  using BaseVectorRV<V,M>::m_unifiedCdf;
  using BaseVectorRV<V,M>::m_mdf;
};

}  // End namespace QUESO

#endif // UQ_INVLOGIT_GAUSSIAN_VECTOR_RV_H
