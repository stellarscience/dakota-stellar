/*  _______________________________________________________________________

 PECOS: Parallel Environment for Creation Of Stochastics
 Copyright (c) 2011, Sandia National Laboratories.
 This software is distributed under the GNU Lesser General Public License.
 For more information, see the README file in the top Pecos directory.
 _______________________________________________________________________ */

#include "Teuchos_SerialDenseHelpers.hpp"
#include "pecos_stat_util.hpp"
#include <iostream>

#include "RosenblattTransformation.hpp"

static const char rcsId[] =
        "@(#) $Id: RosenblattTransformation.cpp 4768 2015-04-07 14:55:32Z franzefn $";

//#define DEBUG

using namespace std;

namespace Pecos {

// -------------------- constructors and desctructors --------------------
RosenblattTransformation::RosenblattTransformation() :
        inversionEpsilon(1e-10) {
#ifdef REFCOUNT_DEBUG
    PCout << "RosenblattTransformation::RosenblattTransformation() called to "
    << "build." << std::endl;
#endif
}

RosenblattTransformation::~RosenblattTransformation() {
    int ndim = marginals.size();
    // do not delete the last marginal since it is the original
    // density estimator
    for (int idim = 0; idim < ndim - 1; idim++) {
        delete marginals[idim];
    }

#ifdef REFCOUNT_DEBUG
    PCout << "RosenblattTransformation deleted"
    << std::endl;
#endif
}

void RosenblattTransformation::initialize(DensityEstimator& density) {
    // prepare the marginalized densities
    densityEstimator = density;
    size_t ndim = densityEstimator.getDim();
    marginals.resize(ndim);
    marginals[ndim - 1] = &densityEstimator;
    for (int idim = ndim - 2; idim >= 0; idim--) {
        marginals[idim] = new DensityEstimator(densityEstimator.getType());
        marginals[idim + 1]->marginalize(idim + 1, *marginals[idim]);
    }
}

// ----------------------------------------------------------------------

Real RosenblattTransformation::getMaxInversionError() {
    return inversionEpsilon;
}

DensityEstimator* RosenblattTransformation::getDensityEstimator() {
    return &densityEstimator;
}

/// Transformation routine from x-space of correlated random variables
/// to u-space of uncorrelated uniform variables
void RosenblattTransformation::trans_X_to_U(const RealVector& x_vars,
        RealVector& u_vars) {
    // helper variables
    size_t ndim = densityEstimator.getDim();

    // 1) do the conditionalization
    vector<DensityEstimator*> cond(ndim);
    cond[0] = marginals[0];
    for (size_t idim = 1; idim < ndim; idim++) {
        cond[idim] = new DensityEstimator(marginals[idim]->getType());
        marginals[idim]->condToDimX(x_vars, idim, *cond[idim]);
    }

    // 2) do the integration
    for (size_t idim = 0; idim < ndim; idim++) {
        u_vars[idim] = trans_X_to_U_1d(x_vars[idim], *cond[idim]);
    }

    // 3) clean up all but cond at position 0
    for (size_t idim = 1; idim < ndim; idim++) {
        delete cond[idim];
    }
}

Real RosenblattTransformation::trans_X_to_U_1d(const Real x,
        DensityEstimator& cond, size_t numQuadIntervals) {
    // -------------------------------------------------------------
    // TODO: better quadrature rule with error control
    // trapezoidal rule
    // -------------------------------------------------------------
    Real y = 0.0;
    Real xlower = -1.0;
    Real dx = (x - xlower) / numQuadIntervals;
    RealVector xi(1);

    // boundary points
    xi[0] = xlower;
    y += cond.pdf(xi);
    xi[0] = x;
    y += cond.pdf(xi);

    // inner points
    xi[0] = xlower + dx;
    while (xi[0] < x) {
        y += 2 * cond.pdf(xi);
        xi[0] += dx;
    }

    y *= dx / 2.0;
    // -------------------------------------------------------------
    return y;
}

/// Transformation routine from u-space of uncorrelated uniform
/// variables to x-space of correlated random variables
void RosenblattTransformation::trans_U_to_X(const RealVector& u_vars,
        RealVector& x_vars) {
    x_vars.putScalar(0.0);

    // helper variables
    size_t ndim = densityEstimator.getDim();
    DensityEstimator* cond = marginals[0];

    for (size_t idim = 0; idim < ndim; idim++) {
        x_vars[idim] = trans_U_to_X_1d(u_vars[idim], *cond);

        // update the conditionalized pdf
        if (idim + 1 < ndim) {
            marginals[idim + 1]->condToDimX(x_vars, idim + 1, *cond);
        }
    }
}

Real RosenblattTransformation::trans_U_to_X_1d(const Real u,
        DensityEstimator& cond, size_t maxIterations) {
    // -------------------------------------------------------------
    // TODO: better root search algorithm
    // bisection
    // -------------------------------------------------------------
    Real xlower = -1e2;
    Real xupper = 1e2;
    Real x = 0.0;
    size_t ii = 0;
    Real ui = 0.0;
    Real xerr = xupper - xlower;

    do {
        // evaluate u_i = F(x_i)
        ui = trans_X_to_U_1d(x, cond);
        // truncate the search space
        if (u < ui) {
            xupper = x;
        } else {
            xlower = x;
        }

        // select new center posize_t
        x = (xlower + xupper) / 2.0;
        xerr = xupper - xlower;
        ii++;
    } while ((xerr > inversionEpsilon) && (ii < maxIterations));

    return x;
}

}
// namespace Pecos
