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

class AleatoryDistParams;


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

  CubatureDriver();                                ///< default constructor
  CubatureDriver(unsigned short integrand_order);  ///< constructor
  ~CubatureDriver();                               ///< destructor

  //
  //- Heading: Member functions
  //

  /// initialize cubature settings except for distribution params
  void initialize_grid(const ShortArray& u_types, unsigned short order,
		       unsigned short rule);
  /// initialize all cubature settings (distribution params already
  /// set within poly_basis)
  void initialize_grid(const std::vector<BasisPolynomial>& poly_basis);

  /// initialize settings for parameterized cubature rules
  void initialize_grid_parameters(const ShortArray& u_types,
				  const AleatoryDistParams& dp);

  /// set integrandOrder
  void integrand_order(unsigned short order);
  /// get integrandOrder
  unsigned short integrand_order() const;

  /// return type1WeightSets
  const RealVector& type1_weight_sets() const;
  // return type2WeightSets
  //const RealMatrix& type2_weight_sets() const;

  // return index of the maximal stored grid state (_NPOS if current grid)
  //size_t maximal_grid() const;

  /// number of collocation points with duplicates removed
  int grid_size();
  /// compute scaled variable and weight sets for the cubature grid
  void compute_grid(RealMatrix& variable_sets);

private:

  //
  //- Heading: Convenience functions
  //

  /// verify that all values within params are identical
  bool verify_homogeneity(const RealVector& params) const;
  /// verify that all vectors within params are identical
  bool verify_homogeneity(const RealRealMapArray& params) const;

  /// size collocRules and set first entry
  void collocation_rule(unsigned short rule);

  //
  //- Heading: Data
  //

  /// integrand order
  unsigned short integrandOrder;
  /// the current number of unique points in the grid
  int numPts;
  /// flag indicating when numPts needs to be recomputed due to an
  /// update to the cubature settings
  bool updateGridSize;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the {TPQ,SSG,Cub} grid
  RealVector type1WeightSets;
  // the set of type2 weights (for integration of gradient interpolants)
  // for each derivative component and for each point in the {TPQ,SSG} grid
  //RealMatrix type2WeightSets;
};


inline CubatureDriver::CubatureDriver(): IntegrationDriver(BaseConstructor()),
  integrandOrder(0), numPts(0), updateGridSize(true)
{ }


inline CubatureDriver::CubatureDriver(unsigned short integrand_order):
  IntegrationDriver(BaseConstructor()), integrandOrder(integrand_order),
  numPts(0), updateGridSize(true)
{ }


inline CubatureDriver::~CubatureDriver()
{ }


inline void CubatureDriver::integrand_order(unsigned short order)
{
  if (integrandOrder != order)
    { integrandOrder = order; updateGridSize = true; }
}


inline unsigned short CubatureDriver::integrand_order() const
{ return integrandOrder; }


inline const RealVector& CubatureDriver::type1_weight_sets() const
{ return type1WeightSets; }


//inline const RealMatrix& CubatureDriver::type2_weight_sets() const
//{ return type2WeightSets; }


inline void CubatureDriver::collocation_rule(unsigned short rule)
{ collocRules.resize(1); collocRules[0] = (short)rule; }


inline bool CubatureDriver::verify_homogeneity(const RealVector& params) const
{
  bool err_flag = false;
  if (!params.empty()) {
    const Real& param0 = params[0];
    for (size_t i=1; i<numVars; ++i)
      if (params[i] != param0)
	{ err_flag = true; break; }
  }
  return err_flag;
}


inline bool CubatureDriver::
verify_homogeneity(const RealRealMapArray& params) const
{
  bool err_flag = false;
  if (!params.empty()) {
    const RealRealMap& param0 = params[0];
    for (size_t i=1; i<numVars; ++i)
      if (params[i] != param0)
	{ err_flag = true; break; }
  }
  return err_flag;
}

} // namespace Pecos

#endif
