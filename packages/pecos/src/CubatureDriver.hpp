/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 CubatureDriver
//- Description: Wrapper class for cubature components within VPISparseGrid
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef CUBATURE_DRIVER_HPP
#define CUBATURE_DRIVER_HPP

#include "IntegrationDriver.hpp"

namespace Pecos {

class MultivariateDistribution;


/// Generates N-dimensional cubature grids for numerical evaluation of
/// expectation integrals over independent standard random variables.

/** Includes Stroud rules and extensions.  This class is used by
    Dakota::NonDCubature, but could also be used for general numerical
    integration of moments. */

class CubatureDriver: public IntegrationDriver
{
public:

  //
  //- Heading: Constructors and destructor
  //

  CubatureDriver();                               ///< default constructor
  CubatureDriver(unsigned short integrand_order); ///< constructor
  ~CubatureDriver();                              ///< destructor

  //
  //- Heading: Virtual function redefinitions
  //

  void initialize_grid(const MultivariateDistribution& mv_dist,
		       unsigned short order, unsigned short rule);
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);
  void initialize_grid_parameters(const MultivariateDistribution& mv_dist);

  const RealMatrix& variable_sets() const;
  const RealVector& type1_weight_sets() const;
  //const RealMatrix& type2_weight_sets() const;
  //const ActiveKey& maximal_grid() const;
  int grid_size();
  void compute_grid();
  void compute_grid(RealMatrix& var_sets);

  //
  //- Heading: Member functions
  //

  /// set integrandOrder
  void integrand_order(unsigned short order);
  /// get integrandOrder
  unsigned short integrand_order() const;

private:

  //
  //- Heading: Convenience functions
  //

  /// verify that all values within params are identical
  template <typename T>
  bool verify_homogeneity(const std::vector<T>& params,
			  const BitArray& active_subset = BitArray()) const;

  /// size collocRules and set first entry
  void collocation_rule(unsigned short rule);

  //
  //- Heading: Data
  //

  /// integrand order
  unsigned short integrandOrder;
  /// the current number of unique points in the grid
  int numPts;

  /// the set of collocation points in the cubature grid
  RealMatrix variableSets;
  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the Cubature grid
  RealVector type1WeightSets;
  // the set of type2 weights (for integration of gradient interpolants)
  // for each derivative component and for each point in the grid
  //RealMatrix type2WeightSets;
};


inline CubatureDriver::CubatureDriver(): IntegrationDriver(BaseConstructor()),
  integrandOrder(0), numPts(0)
{ }


inline CubatureDriver::CubatureDriver(unsigned short integrand_order):
  IntegrationDriver(BaseConstructor()), integrandOrder(integrand_order),
  numPts(0)
{ }


inline CubatureDriver::~CubatureDriver()
{ }


inline void CubatureDriver::compute_grid(RealMatrix& var_sets)
{
  compute_grid();
  var_sets = variableSets; // copy
}


inline void CubatureDriver::integrand_order(unsigned short order)
{
  if (integrandOrder != order) {
    integrandOrder = order;
    numPts = 0; // special value indicates update required
  }
}


inline unsigned short CubatureDriver::integrand_order() const
{ return integrandOrder; }


inline const RealMatrix& CubatureDriver::variable_sets() const
{ return variableSets; }


inline const RealVector& CubatureDriver::type1_weight_sets() const
{ return type1WeightSets; }


//inline const RealMatrix& CubatureDriver::type2_weight_sets() const
//{ return type2WeightSets; }


inline void CubatureDriver::collocation_rule(unsigned short rule)
{ collocRules.resize(1); collocRules[0] = (short)rule; }


template <typename T>
bool CubatureDriver::
verify_homogeneity(const std::vector<T>& params,
		   const BitArray& active_subset) const
{
  bool err_flag = false;
  size_t len = params.size();
  if (len <= 1) return false;

  if (active_subset.empty()) { // no subset
    const T& param0 = params[0];
    for (size_t i=1; i<len; ++i)
      if (params[i] != param0)
	{ err_flag = true; break; }
  }
  else {
    size_t i, first_i = len;
    for (i=0; i<len; ++i)
      if (active_subset[i])
	{ first_i = i; break; }
    if (first_i < len) {
      const T& param = params[first_i];
      for (i=first_i+1; i<len; ++i)
	if (active_subset[i] && params[i] != param)
	  { err_flag = true; break; }
    }
  }
  return err_flag;
}

} // namespace Pecos

#endif
