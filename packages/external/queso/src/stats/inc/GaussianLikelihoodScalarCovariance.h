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

#ifndef UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H
#define UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H

#include <queso/LikelihoodBase.h>

namespace QUESO {

class GslVector;
class GslMatrix;

/*!
 * \file GaussianLikelihoodScalarCovariance.h
 *
 * \class GaussianLikelihoodScalarCovariance
 * \brief A class that represents a Gaussian likelihood with scalar covariance
 */

template <class V = GslVector, class M = GslMatrix>
class GaussianLikelihoodScalarCovariance : public LikelihoodBase<V, M> {
public:
  //! @name Constructor/Destructor methods.
  //@{
  //! Default constructor.
  /*!
   * Instantiates a Gaussian likelihood function, given a prefix, its domain,
   * a set of observations and a scalar covariance matrix.  The scalar
   * 'covariance matrix' is just passed as a \c double.
   */
  GaussianLikelihoodScalarCovariance(const char * prefix,
      const VectorSet<V, M> & domainSet, const V & observations,
      double covariance);

  //! Destructor
  virtual ~GaussianLikelihoodScalarCovariance();
  //@}

  //! Logarithm of the value of the scalar function.
  virtual double lnValue(const V & domainVector) const;

  using LikelihoodBase<V, M>::lnValue;

private:
  double m_covariance;
};

}  // End namespace QUESO

#endif  // UQ_GAUSSIAN_LIKELIHOOD_SCALAR_COV_H
