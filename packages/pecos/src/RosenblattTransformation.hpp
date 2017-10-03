/*  _______________________________________________________________________

 PECOS: Parallel Environment for Creation Of Stochastics
 Copyright (c) 2011, Sandia National Laboratories.
 This software is distributed under the GNU Lesser General Public License.
 For more information, see the README file in the top Pecos directory.
 _______________________________________________________________________ */

#ifndef ROSENBLATT_TRANSFORMATION_HPP
#define ROSENBLATT_TRANSFORMATION_HPP

#include "ProbabilityTransformation.hpp"
#include "pecos_data_types.hpp"
#include "DensityEstimator.hpp"

namespace Pecos {

/// Class for Rosenblatt nonlinear distribution transformation.
class RosenblattTransformation: public ProbabilityTransformation {
public:

    //
    //- Heading: Constructors and destructor
    //

    /// default constructor
    RosenblattTransformation();
    /// destructor
    virtual ~RosenblattTransformation();

    void initialize(DensityEstimator& density);

    //
    //- Heading: Virtual function redefinitions
    //

    /// Transformation routine from x-space of correlated random variables
    /// to u-space of uncorrelated standard normal variables
    virtual void trans_X_to_U(const RealVector& x_vars, RealVector& u_vars);

    /// Transformation routine from u-space of uncorrelated standard normal
    /// variables to x-space of correlated random variables
    virtual void trans_U_to_X(const RealVector& u_vars, RealVector& x_vars);

    /// get the maximum error made during inversion
    Real getMaxInversionError();

    /// get the density estimator on which the Rosenblatt transformation is
    /// defined on
    DensityEstimator* getDensityEstimator();

private:
    Real trans_X_to_U_1d(const Real x, DensityEstimator& cond,
            size_t numQuadIntervals = 100);
    Real trans_U_to_X_1d(const Real u, DensityEstimator& cond,
            size_t maxIterations = 20);

    DensityEstimator densityEstimator;
    std::vector<DensityEstimator*> marginals;

protected:
    /// maximum allowed inversion error
    Real inversionEpsilon;
};

} // namespace Pecos

#endif
