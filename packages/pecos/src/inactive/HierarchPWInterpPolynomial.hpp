/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        HierarchicalBasis
//- Description:  Class for 1-D HierarchicalLinearBasis
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#ifndef HIERARCHICABASIS_HPP
#define HIERARCHICABASIS_HPP

#include "PiecewiseInterpPolynomial.hpp"
#include "pecos_data_types.hpp"
#include "RefinablePointSet.hpp"

namespace Pecos {

/// Derived basis for 1-D piecewise linear interpolants expressed in a
/// hierarchical basis.

/** The HierarchicalBasis class evaluates a univariate piecewise
    polynomial interpolation function defined on a hierarchical grid.
    The grid is provided by the RefinablePointSet passed to the
    constructor.

    For additional information see:
    X. Ma and N. Zabras, "An adaptive hierarchical sparse grid 
    collocation algorithm for the solution of stochastic 
    differential equations", J. Comput. Phys. 228:3084--3113, 2009.
*/

class HierarchPWInterpPolynomial: public PiecewiseInterpPolynomial
{
public:

  //
  //- Heading: Constructor and Destructor
  // 

  // These forms preferred for use from BasisPolynomial::get_polynomial():
  //
  // Default constructor (can't do due to reference init)
  //HierarchPWInterpPolynomial();
  // Standard constructor (can't do due to reference init)
  //HierarchPWInterpPolynomial(short basisPolyType_,
  //			       short collocRule_ = CLENSHAW_CURTIS);

  /// constructor with point set
  HierarchPWInterpPolynomial(RefinablePointSet& pointSet_,
			     short basisPolyType_ = PIECEWISE_LINEAR_INTERP,
			     short    collocRule_ = CLENSHAW_CURTIS);

  /// Destructor
  virtual ~HierarchPWInterpPolynomial();

  //
  //- Heading: Virtual function redefinitions
  //

  Real type1_value(const Real& x, const unsigned int i);
  Real type2_value(const Real& x, const unsigned int i);

  Real type1_gradient(const Real& x, const unsigned int i);
  Real type2_gradient(const Real& x, const unsigned int i);
       
protected:
    
  //
  //- Heading: Data
  //

  /// Grid of definition
  RefinablePointSet& pointSet;

private:

}; //End class definition


//inline HierarchPWInterpPolynomial::HierarchPWInterpPolynomial():
//  PiecewiseInterpPolynomial()
//{ }

//inline HierarchPWInterpPolynomial::
//HierarchPWInterpPolynomial(short basisPolyType_, short collocRule_):
//  PiecewiseInterpPolynomial(basisPolyType_, collocRule_))
//{ }

/// throws an exception if left_end_ >= right_end_.
inline HierarchPWInterpPolynomial::
HierarchPWInterpPolynomial(RefinablePointSet& pointSet_, short basisPolyType_,
			   short collocRule_):
  PiecewiseInterpPolynomial(pointSet_.get_interp_points(), basisPolyType_,
			    collocRule_),
  pointSet(pointSet_)
{ }

inline HierarchPWInterpPolynomial::~HierarchPWInterpPolynomial()
{ }

} //End namespace Pecos

#endif
