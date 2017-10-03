/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        LRDHierarchInterpPolyApproximation
//- Description:  Implementation code for InterpPolyApproximation class
//-               
//- Owner:        Chris Miller

#include "LRDHierarchInterpPolyApproximation.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

namespace Pecos {


LRDHierarchInterpPolyApproximation::
LRDHierarchInterpPolyApproximation(short basis_type, size_t num_vars,
				bool use_derivs):
  InterpPolyApproximation(basis_type,num_vars,use_derivs)
{}


LRDHierarchInterpPolyApproximation::~LRDHierarchInterpPolyApproximation()
{}


Real LRDHierarchInterpPolyApproximation::value(const RealVector& x)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "LRDHierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }
  Real approx_val = 0.;
  int num_colloc_points = driverRep->grid_size();

  //Grab the supports of the hierarchical basis functions.
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();

  //Determine which basis functions share their support with x.
  const IntArray& supportIndicator = in_support_of(x);
  const std::vector<CollocationPoint>& colloc_pts = lr_driver->get_collocation_points();
  switch ( basisConfigOptions.useDerivs ) {
  case false:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
	
      Real local_val = 1.0;
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
      for ( unsigned int dim_idx = 0; dim_idx < numVars ; ++dim_idx ) {
	const unsigned int this_dim_level = level_index[dim_idx][0];
	//const unsigned int this_dim_index = level_index[dim_idx][1];
	if ( this_dim_level == 1 )
	  local_val *= 1;  //constant case
	else { //Linear hat functions
	  if ( x[dim_idx] == point[dim_idx] )
	    local_val *= 1;
	  else if (x[dim_idx] < point[dim_idx])
	    local_val *= ( x[dim_idx] - this_point_support[0][dim_idx] ) / 
	      ( point[dim_idx] - this_point_support[0][dim_idx] );
	  else
	    local_val *= ( this_point_support[1][dim_idx] - x[dim_idx] ) /
	      ( this_point_support[1][dim_idx] - point[dim_idx] );
	}
      }
      approx_val += expansionType1Coeffs[supportIndicator[i]] * local_val;
    }
    break;
  case true:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i ) {
	
      Real local_val = 1.0;
      const Real*     coeff2_i = expansionType2Coeffs[supportIndicator[i]]; // column vector
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
      RealVector terms(numVars+1);
      terms[0] = expansionType1Coeffs[supportIndicator[i]];
      for ( unsigned int j = 0; j < numVars; ++j) {
	terms[j+1] = coeff2_i[j];
      }
      for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	const unsigned int this_dim_level = level_index[dim_idx][0];
	const unsigned int this_dim_index = level_index[dim_idx][1];
	if ( this_dim_level == 1 ) {
	  terms[0] *= 1;
	  if (x[dim_idx] == point[dim_idx]) {
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 0 : 1;
	    }
	  } else if (x[dim_idx] < point[dim_idx]) {
	    const Real t = (x[dim_idx] - this_point_support[0][dim_idx])/
	      (point[dim_idx] - this_point_support[0][dim_idx]);
	    const Real dx_by_dt = 
	      point[dim_idx]-this_point_support[0][dim_idx];
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		t*t*(t-1)*dx_by_dt : 1;
	    }
	  } else {
	    const Real t = (x[dim_idx] - point[dim_idx])/
	      (this_point_support[1][dim_idx] - point[dim_idx]);
	    const Real dx_by_dt = 
	      (this_point_support[1][dim_idx]-point[dim_idx]);
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		t*(1-t)*(1-t)*dx_by_dt : 1;
	    }
	  }
	} else {
	  if (x[dim_idx] == point[dim_idx]) {
	    terms[0] *= 1;
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 0 : 1;
	    }
	  } else if (x[dim_idx] < point[dim_idx]) {
	    const Real t = (x[dim_idx] - this_point_support[0][dim_idx])/
	      (point[dim_idx] - this_point_support[0][dim_idx]);
	    const Real dx_by_dt = (point[dim_idx] - this_point_support[0][dim_idx]);
	    terms[0] *= t*t*(3-2*t);
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? t*t*(t-1)*dx_by_dt : t*t*(3-2*t);
	    }
	  } else {
	    const Real t = (x[dim_idx] - point[dim_idx])/
	      (this_point_support[1][dim_idx] - point[dim_idx]);
	    const Real dx_by_dt = (this_point_support[1][dim_idx] - point[dim_idx]);
	    terms[0] *= (1+2*t)*(1-t)*(1-t);
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? t*(1-t)*(1-t)*dx_by_dt : (1+2*t)*(1-t)*(1-t);
	    }
	  }
	}
      }
      for ( unsigned int accumulator = 0; accumulator < terms.length(); ++accumulator ) {
	approx_val += terms[accumulator];
      }
      //Real const * fakeout = terms.RealMatrix::operator[](0);
      //approx_val += std::accumulate(fakeout, fakeout+terms.length(), 0);
    }
    break;
  }
  return approx_val;
}
  

const RealVector& LRDHierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "LRDHierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }
  approxGradient.size(numVars);
  for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx )
    approxGradient[dim_idx] = 0; 
    
  int num_colloc_points = driverRep->grid_size();

  //Grab the supports of the hierarchical basis functions.
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();

  //Determine which basis functions share their support with x.
  const IntArray& supportIndicator = in_support_of(x);
  const std::vector<CollocationPoint>& colloc_pts = lr_driver->get_collocation_points();
  switch ( basisConfigOptions.useDerivs ) {
  case false:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
      for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	Real localGrad_i = 1;
	const unsigned int this_dim_level = level_index[dim_idx][0];
	const unsigned int this_dim_index = level_index[dim_idx][1];
	if ( this_dim_level == 1 ) {
	  localGrad_i *= 0;  //constant case
	} else { //Linear hat functions
	  if ( (x[dim_idx] == point[dim_idx]) || 
	       (x[dim_idx] == this_point_support[0][dim_idx]) || 
	       (x[dim_idx] == this_point_support[1][dim_idx]) ) {
	    localGrad_i *= 0;
	  } else if (x[dim_idx] < point[dim_idx]) {
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      localGrad_i *= (dim_idx == dim_idx2) ? 
		1 / ( point[dim_idx] - this_point_support[0][dim_idx] ) :
		(x[dim_idx2] - this_point_support[0][dim_idx2]) / 
		( point[dim_idx2] - this_point_support[0][dim_idx2] );
	    }
	  } else {
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      localGrad_i *= (dim_idx == dim_idx2) ? 
		-1 / ( this_point_support[1][dim_idx2] - point[dim_idx2] ) :
		(point[dim_idx2] - x[dim_idx2])/( this_point_support[1][dim_idx2] - point[dim_idx2] );
	    }
	  }   
	}
	approxGradient[dim_idx] = approxGradient[dim_idx] + 
	  expansionType1Coeffs[supportIndicator[i]] * localGrad_i;
      }
    }
    break;
  case true:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      const Real*     coeff2_i = expansionType2Coeffs[supportIndicator[i]];
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
      // These nested loops are filthy.
      // dim_idx keeps track of which gradient entry we're filling.
      // dim_idx2 keeps track of which variable we're thinking about.
      // dim_idx3 tess us which term we're working on. 
      for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	RealVector terms(numVars+1);
	terms[0] = expansionType1Coeffs[supportIndicator[i]];
	for ( unsigned int j = 0; j < numVars; ++j ) {
	  terms[j+1] = coeff2_i[j];
	}
	for ( unsigned int dim_idx2 = 0; dim_idx2< numVars; ++dim_idx2 ) { 
	  const unsigned int this_dim_level = level_index[dim_idx2][0];
	  const unsigned int this_dim_index = level_index[dim_idx2][1];
	  if ( this_dim_level == 1 ) {
	    terms[0] *= (dim_idx == dim_idx2) ? 0 : 1;  
	    if ( x[dim_idx2] == point[dim_idx2] ) {
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  1 : 0;  // Take derivative of basis else don't
		} else {  //Use type 1.
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else if ( x[dim_idx2] < point[dim_idx2] ) {
	      const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		(point[dim_idx2] - this_point_support[0][dim_idx2]);
	      const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		    (3*t*t-2*t) : 
		    dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		    0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else {
	      const Real t = (x[dim_idx2] - point[dim_idx2])/
		(this_point_support[1][dim_idx2] - point[dim_idx2]);
	      const Real dx_by_dt = this_point_support[1][dim_idx2] - point[dim_idx2];
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		    (3*t*t - 4*t + 1) : 
		    dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		    0 : 1; // Take derivative of basis else don't
		}
	      }
	    }
	  } else {  
	    if ( x[dim_idx2] == point[dim_idx2] ) {
	      terms[0] *= (dim_idx == dim_idx2) ? 0 : 1;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  1 : 0;  // Take derivative of basis else don't
		} else {  //Use type 1.
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else if ( x[dim_idx2] < point[dim_idx2] ) {
	      const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		(point[dim_idx2] - this_point_support[0][dim_idx2]);
	      const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      terms[0] *= (dim_idx == dim_idx2) ?  dt_by_dx*(-6*t*t + 6*t) : -2*t*t*t + 3*t*t;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		    (3*t*t-2*t) : 
		    dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		    dt_by_dx*(-6*t*t+6*t) : t*t*(3-2*t); // Take derivative of basis else don't
		}
	      }
	    } else {
	      const Real t = (x[dim_idx2] - point[dim_idx2])/
		(this_point_support[1][dim_idx2] - point[dim_idx2]);
	      const Real dx_by_dt = this_point_support[0][dim_idx2] - point[dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      terms[0] *= (dim_idx == dim_idx2) ?  dt_by_dx*(6*t*t - 6*t) : 2*t*t*t - 3*t*t + 1;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		    (3*t*t - 4*t + 1) : 
		    dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		    dt_by_dx*(6*t*t-6*t) : (1+2*t)*(1-t)*(1-t); // Take derivative of basis else don't
		}
	      }
	    }
	  }
	}
	for ( unsigned int accumulator = 0; accumulator < terms.length(); ++accumulator ) {
	  approxGradient[dim_idx] += terms[accumulator];
	}
	//approxGradient[dim_idx] += std::accumulate(terms.RealMatrix::operator[](0),
	//					     terms.RealMatrix::operator[](0) + terms.length(),
	//					     0);
      }
    }
    break;
  }
  return approxGradient;
}


const RealVector& LRDHierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, const SizetArray& dvv)
{
  size_t num_deriv_vars = dvv.size();
  if ( approxGradient.length() != num_deriv_vars )
    approxGradient.sizeUninitialized(num_deriv_vars);
    
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "LRDHierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }
    
  for ( unsigned int dim_idx = 0; dim_idx < num_deriv_vars; ++dim_idx )
    approxGradient[dim_idx] = 0; 
    
  int num_colloc_points = driverRep->grid_size();

  //Grab the supports of the hierarchical basis functions.
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();

  //Determine which basis functions share their support with x.
  const IntArray& supportIndicator = in_support_of(x);
  const std::vector<CollocationPoint>& colloc_pts = lr_driver->get_collocation_points();
  switch ( basisConfigOptions.useDerivs ) {
  case false:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
      for ( unsigned int dim_idx = 0; dim_idx < num_deriv_vars; ++dim_idx ) {
	Real localGrad_i = 1;
	const size_t deriv_index = dvv[dim_idx] - 1;
	const unsigned int this_dim_level = level_index[deriv_index][0];
	const unsigned int this_dim_index = level_index[deriv_index][1];
	if ( this_dim_level == 1 ) {
	  localGrad_i *= 0;  //constant case
	} else { //Linear hat functions
	  if ( (x[deriv_index] == point[deriv_index]) || 
	       (x[deriv_index] == this_point_support[0][deriv_index]) || 
	       (x[deriv_index] == this_point_support[1][deriv_index]) ) {
	    localGrad_i *= 0;
	  } else if (x[deriv_index] < point[deriv_index]) {
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      localGrad_i *= (deriv_index == dim_idx2) ? 
		1 / ( point[deriv_index] - this_point_support[0][deriv_index] ) :
		(x[dim_idx2] - this_point_support[0][dim_idx2]) / 
		( point[dim_idx2] - this_point_support[0][dim_idx2] );
	    }
	  } else {
	    for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
	      localGrad_i *= (deriv_index == dim_idx2) ? 
		-1 / ( this_point_support[1][deriv_index] - point[deriv_index] ) :
		(point[dim_idx2] - x[dim_idx2])/( this_point_support[1][dim_idx2] - point[dim_idx2] );
	    }
	  }   
	}
	approxGradient[dim_idx] = approxGradient[dim_idx] + 
	  expansionType1Coeffs[supportIndicator[i]] * localGrad_i;
      }
    }
    break;
  case true:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      const Real*     coeff2_i = expansionType2Coeffs[supportIndicator[i]];
      const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
      const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
      const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
      // These nested loops are filthy.
      // dim_idx keeps track of which gradient entry we're filling.
      // dim_idx2 keeps track of which variable we're thinking about.
      // dim_idx3 tess us which term we're working on. 
      for ( unsigned int dim_idx = 0; dim_idx < num_deriv_vars; ++dim_idx ) {
	RealVector terms(numVars+1);
	const size_t deriv_index = dvv[dim_idx] - 1;
	terms[0] = expansionType1Coeffs[supportIndicator[i]];
	for ( unsigned int j = 0; j < numVars; ++j ) {
	  terms[j+1] = coeff2_i[j];
	}
	for ( unsigned int dim_idx2 = 0; dim_idx2< numVars; ++dim_idx2 ) { 
	  const unsigned int this_dim_level = level_index[dim_idx2][0];
	  const unsigned int this_dim_index = level_index[dim_idx2][1];
	  if ( this_dim_level == 1 ) {
	    terms[0] *= (deriv_index == dim_idx2) ? 0 : 1;  
	    if ( x[dim_idx2] == point[dim_idx2] ) {
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  1 : 0;  // Take derivative of basis else don't
		} else {  //Use type 1.
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else if ( x[dim_idx2] < point[dim_idx2] ) {
	      const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		(point[dim_idx2] - this_point_support[0][dim_idx2]);
	      const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  
		    (3*t*t-2*t) : 
		    dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 
		    0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else {
	      const Real t = (x[dim_idx2] - point[dim_idx2])/
		(this_point_support[1][dim_idx2] - point[dim_idx2]);
	      const Real dx_by_dt = this_point_support[1][dim_idx2] - point[dim_idx2];
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  
		    (3*t*t - 4*t + 1) : 
		    dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 
		    0 : 1; // Take derivative of basis else don't
		}
	      }
	    }
	  } else {  
	    if ( x[dim_idx2] == point[dim_idx2] ) {
	      terms[0] *= (deriv_index == dim_idx2) ? 0 : 1;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  1 : 0;  // Take derivative of basis else don't
		} else {  //Use type 1.
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 0 : 1; // Take derivative of basis else don't
		}
	      }
	    } else if ( x[dim_idx2] < point[dim_idx2] ) {
	      const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		(point[dim_idx2] - this_point_support[0][dim_idx2]);
	      const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      terms[0] *= (deriv_index == dim_idx2) ?  dt_by_dx*(-6*t*t + 6*t) : -2*t*t*t + 3*t*t;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  
		    (3*t*t-2*t) : 
		    dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 
		    dt_by_dx*(-6*t*t+6*t) : t*t*(3-2*t); // Take derivative of basis else don't
		}
	      }
	    } else {
	      const Real t = (x[dim_idx2] - point[dim_idx2])/
		(this_point_support[1][dim_idx2] - point[dim_idx2]);
	      const Real dx_by_dt = this_point_support[0][dim_idx2] - point[dim_idx2];
	      const Real dt_by_dx = 1/dx_by_dt;
	      terms[0] *= (deriv_index == dim_idx2) ?  dt_by_dx*(6*t*t - 6*t) : 2*t*t*t - 3*t*t + 1;
	      for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ?  
		    (3*t*t - 4*t + 1) : 
		    dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		} else { //use type 1
		  terms[dim_idx3+1] *= (dim_idx2 == deriv_index) ? 
		    dt_by_dx*(6*t*t-6*t) : (1+2*t)*(1-t)*(1-t); // Take derivative of basis else don't
		}
	      }
	    }
	  }
	}
	for ( unsigned int accumulator = 0; accumulator < terms.length(); ++accumulator ) {
	  approxGradient[dim_idx] += terms[accumulator];
	}
      }
    }
    break;
  }
  return approxGradient;
}

const RealVector& LRDHierarchInterpPolyApproximation::
gradient_nonbasis_variables(const RealVector& x)
{ return approxGradient; /* TO DO */ }

Real LRDHierarchInterpPolyApproximation::stored_value(const RealVector& x)
{ return 0.; /* TO DO */ }

const RealVector& LRDHierarchInterpPolyApproximation::
stored_gradient_basis_variables(const RealVector& x)
{ return approxGradient; /* TO DO */ }

const RealVector& LRDHierarchInterpPolyApproximation::
stored_gradient_nonbasis_variables(const RealVector& x)
{ return approxGradient; /* TO DO */ }

Real LRDHierarchInterpPolyApproximation::mean()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::mean()" << std::endl;
    abort_handler(-1);
  }
  Real mean = 0.;
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  switch (basisConfigOptions.useDerivs) {
  case false:
    for ( unsigned int i = 0; i<numCollocPts; ++i)
      mean += expansionType1Coeffs[i] * t1_wts[i];
    break;
  case true:
    const RealMatrix& t2_wts = driverRep->type2_weight_sets();
    for ( unsigned int i = 0; i< numCollocPts; ++i) {
      mean += expansionType1Coeffs[i] * t1_wts[i];
      const Real* coeff2_i = expansionType2Coeffs[i];
      const Real* t2_wt_i = t2_wts[i];
      for (unsigned int j = 0; j< numVars; ++j)
	mean += coeff2_i[j] * t2_wt_i[j];
    }
    break;
  }
  return mean;
}

Real LRDHierarchInterpPolyApproximation::mean(const RealVector& x)
{
  std::cout << "TODO: mean in all variables mode";
  return numericalMoments[0];
}

const RealVector& LRDHierarchInterpPolyApproximation::
mean_gradient()
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffGradFlag) {
    PCerr << "Error: expansion coefficient gradients not defined in Nodal"
	  << "InterpPolyApproximation::mean_gradient()." << std::endl;
    abort_handler(-1);
  }
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  size_t i, j, num_deriv_vars = expansionType1CoeffGrads.numRows();
  if (meanGradient.length() != num_deriv_vars)
    meanGradient.sizeUninitialized(num_deriv_vars);
  meanGradient = 0.;
  for (i=0; i<numCollocPts; ++i) {
    const Real& t1_wt_i = t1_wts[i];
    for (j=0; j<num_deriv_vars; ++j)
      meanGradient[j] += expansionType1CoeffGrads(j,i) * t1_wt_i;
  }
  return meanGradient;
}

const RealVector& LRDHierarchInterpPolyApproximation::
mean_gradient(const RealVector& x, const SizetArray& dvv)
{
  std::cout << "TODO: mean_gradient in all variables mode" << std::endl;
  return meanGradient;
}

Real LRDHierarchInterpPolyApproximation::variance()
{ return covariance(this); }

Real LRDHierarchInterpPolyApproximation::variance(const RealVector& x)
{
  //TODO
  PCerr << "TODO: variance in all variables mode" << std::endl;
  return numericalMoments[1];
}

const RealVector& LRDHierarchInterpPolyApproximation::variance_gradient()
{
  //TODO
  PCerr << "TODO: variance_gradient()" << std::endl;
  return varianceGradient;
}

const RealVector& LRDHierarchInterpPolyApproximation::
variance_gradient(const RealVector& x, const SizetArray& dvv)
{
  //TODO
  PCerr << "TODO: variance_gradient in all variables mode" << std::endl;
  return varianceGradient;
}

Real LRDHierarchInterpPolyApproximation::
covariance(PolynomialApproximation* poly_approx_2)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "NodalInterpPolyApproximation::covariance()" << std::endl;
    abort_handler(-1);
  }
  LRDHierarchInterpPolyApproximation* hip_approx_2 = 
    static_cast<LRDHierarchInterpPolyApproximation*>(poly_approx_2);
  Real mean_1 = mean(), mean_2 = hip_approx_2->mean();
  const RealVector& t1_coeffs_2 = hip_approx_2->expansionType1Coeffs;
  const RealVector& t1_wts = driverRep->type1_weight_sets();
  Real covar = 0.0;
  switch (basisConfigOptions.useDerivs) {
  case false:
    for ( unsigned int i=0; i<numCollocPts; ++i)
      covar += (expansionType1Coeffs[i] - mean_1) * (t1_coeffs_2[i] - mean_2)
	*  t1_wts[i];
    break;
  case true:
    const RealMatrix& t2_coeffs_2 = hip_approx_2->expansionType2Coeffs;
    const RealMatrix& t2_wts = driverRep->type2_weight_sets();
    for ( unsigned int i = 0; i < numCollocPts; ++i) {
      // type1 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      Real coeff1_i_mm1 = expansionType1Coeffs[i] - mean_1;
      Real coeff1_2i_mm2 = t1_coeffs_2[i]          - mean_2;
      covar += coeff1_i_mm1 * coeff1_2i_mm2 * t1_wts[i];
      // type2 interpolation of (R_1 - \mu_1) (R_2 - \mu_2)
      // --> interpolated gradients are (R_1-\mu_1) * R_2' + (R_2-\mu_2) * R_1'
      const Real *coeff2_i  = expansionType2Coeffs[i];
      const Real *coeff2_2i = t2_coeffs_2[i], *t2_wt_i = t2_wts[i];
      for (unsigned int j=0; j<numVars; ++j)
	covar  += (coeff1_i_mm1 * coeff2_2i[j] + coeff1_2i_mm2 * coeff2_i[j])
	  *  t2_wt_i[j];
    }
    break;
  }
  return covar;
}

Real LRDHierarchInterpPolyApproximation::
covariance(const RealVector& x, PolynomialApproximation* poly_approx_2)
{
  //TODO
  PCerr << "TODO: covariance in all variables" << std::endl;
  Real var;
  return var;
}

const IntArray& LRDHierarchInterpPolyApproximation::
in_support_of(const RealVector& x)
{
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  supportIndicator.clear();
  supportIndicator.reserve( driverRep->grid_size() );
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();
  bool is_in_support;

  for ( unsigned int idx = 0; idx < numCollocPts; ++idx ) {
    is_in_support = true;
    for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
      if ( ( x[dim_idx] > supports[idx][1][dim_idx] ) || ( x[dim_idx] < supports[idx][0][dim_idx] ) ) {
	is_in_support = false;
	break;
      }
    }
    if ( is_in_support ) supportIndicator.push_back(idx);
  }
  return supportIndicator;
}


void LRDHierarchInterpPolyApproximation::allocate_expansion_coefficients()
{
  //std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  expansionType1Coeffs.resize(surrData.size());
  if (basisConfigOptions.useDerivs)
    expansionType2Coeffs.shapeUninitialized(numVars,surrData.size());
}


void LRDHierarchInterpPolyApproximation::compute_expansion_coefficients()
{
  //std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  assert( driverRep != NULL );
  assert(numCollocPts == driverRep->grid_size());
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<CollocationPoint>& col_pts = 
    lr_driver->get_collocation_points();
  numVars = surrData.continuous_variables(0).length();
  expansionType1Coeffs[0] = surrData.response_function(0);
  Teuchos::setCol(surrData.response_gradient(0),0,expansionType2Coeffs);
  unsigned int level = 1;
  for (int i=1; i<numCollocPts; ++i) {
    if (expConfigOptions.expansionCoeffFlag) {
      const Int2DArray& level_index = col_pts[i].get_level_index();
      const RealVector& point = col_pts[i].get_point();
      unsigned int this_point_level = col_pts[i].get_level();
      if ( this_point_level > level ) {
	level = this_point_level;
      }
      // for the Hierarchical case the coefficient is the function value 
      //minus the intepolant value from the
      // previous approximation level.  This is the so called 'hierarchical 
      //surplus.'
      expansionType1Coeffs[i] = surrData.response_function(i) - 
	this->value(surrData.continuous_variables(i),level - 1);
      if (basisConfigOptions.useDerivs) {
	RealVector trueGrad = surrData.response_gradient(i);
	RealVector approxGrad = this->
	  gradient_basis_variables(surrData.continuous_variables(i), level-1);
	for ( unsigned int idx = 0; idx < numVars; ++idx ){
	  trueGrad[idx] -= approxGrad[idx];
	}
	Teuchos::setCol(trueGrad, i, expansionType2Coeffs);
      }
    }
    //if (expConfigOptions.expansionCoeffGradFlag)
    //Teuchos::setCol(it->response_gradient(), i, expansionType1CoeffGrads);
  }
  maxComputedCoeff = numCollocPts-1;
}

void LRDHierarchInterpPolyApproximation::restore_expansion_coefficients()
{
  size_t new_colloc_pts = surrData.size();
  /*
  if (expConfigOptions.expansionCoeffFlag) {
    expansionType1Coeffs.resize(new_colloc_pts);
    if (basisConfigOptions.useDerivs) {
      size_t num_deriv_vars = expansionType2Coeffs.numRows();
      expansionType2Coeffs.reshape(num_deriv_vars, new_colloc_pts);
    }
  }
  if (expConfigOptions.expansionCoeffGradFlag) {
    size_t num_deriv_vars = expansionType1CoeffGrads.numRows();
    expansionType1CoeffGrads.reshape(num_deriv_vars, new_colloc_pts);
  }

  for (int i=numCollocPts; i<new_colloc_pts; ++i) {
    if (expConfigOptions.expansionCoeffFlag) {
      expansionType1Coeffs[i] = surrData.response_function(i);
      if (basisConfigOptions.useDerivs)
	Teuchos::setCol(surrData.response_gradient(i), i, expansionType2Coeffs);
    }
    if (expConfigOptions.expansionCoeffGradFlag)
      Teuchos::setCol(surrData.response_gradient(i),i,expansionType1CoeffGrads);
  }
  */
  numCollocPts = new_colloc_pts;
}

void LRDHierarchInterpPolyApproximation::increment_coefficients()
{
  //std::vector<SurrogateDataPoint>::iterator it = dataPoints.begin();
  expansionType1Coeffs.resize(surrData.size());
  numCollocPts = driverRep->grid_size();
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<CollocationPoint>& col_pts = 
    lr_driver->get_collocation_points();
  if ( maxComputedCoeff == numCollocPts ) return;
  assert(numCollocPts == surrData.size());    
  unsigned int level = 1;
  for (int i=maxComputedCoeff+1; i<numCollocPts; ++i) {
    if (expConfigOptions.expansionCoeffFlag) {
      const Int2DArray& level_index = col_pts[i].get_level_index();
      const RealVector& point = col_pts[i].get_point();
      unsigned int this_point_level = col_pts[i].get_level();
      if ( this_point_level > level ) {
	level = this_point_level;
      }
      // for the Hierarchical case the coefficient is the function value 
      //minus the intepolant value from the previous approximation level.  
      //This is the so called 'hierarchical surplus.'
      expansionType1Coeffs[i] = surrData.response_function(i) - 
	this->value(surrData.continuous_variables(i),level - 1);
      if (basisConfigOptions.useDerivs) {
	RealVector trueGrad = surrData.response_gradient(i);
	RealVector approxGrad = this->
	  gradient_basis_variables(surrData.continuous_variables(i), level-1);
	for ( unsigned int idx = 0; idx < numVars; ++idx ){
	  trueGrad[idx] -= approxGrad[idx];
	}
	Teuchos::setCol(trueGrad, i, expansionType2Coeffs);
      }
    }
    //if (expConfigOptions.expansionCoeffGradFlag)
    //Teuchos::setCol(it->response_gradient(), i, expansionType1CoeffGrads);
  }
  maxComputedCoeff = numCollocPts-1;
}

Real LRDHierarchInterpPolyApproximation::
value(const RealVector& x, unsigned int max_level)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "LRDHierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }
  Real approx_val = 0.;
  int num_colloc_points = driverRep->grid_size();

  //Grab the supports of the hierarchical basis functions.
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();

  //Determine which basis functions share their support with x.
  const IntArray& supportIndicator = in_support_of(x);
  const std::vector<CollocationPoint>& colloc_pts = lr_driver->get_collocation_points();
  switch ( basisConfigOptions.useDerivs ) {
  case false:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
	
      unsigned int this_point_level = colloc_pts[supportIndicator[i]].get_level();
      if ( this_point_level <= max_level ) {
	Real local_val = 1.0;
	const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
	const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
	const Real2DArray& this_point_support = supports[supportIndicator[i]];
	
	for ( unsigned int dim_idx = 0; dim_idx < numVars ; ++dim_idx ) {
	  const unsigned int this_dim_level = level_index[dim_idx][0];
	  const unsigned int this_dim_index = level_index[dim_idx][1];
	  if ( this_dim_level == 1 ) {
	    local_val *= 1;  //constant case
	  } else { //Linear hat functions
	    if ( x[dim_idx] == point[dim_idx] ) {
	      local_val *= 1;
	    } else if (x[dim_idx] < point[dim_idx]) {
	      local_val *= ( x[dim_idx] - this_point_support[0][dim_idx] ) / 
		( point[dim_idx] - this_point_support[0][dim_idx] );
	    } else {
	      local_val *= ( this_point_support[1][dim_idx] - x[dim_idx] ) /
		( this_point_support[1][dim_idx] - point[dim_idx] );
	    }
	  }   
	}
	approx_val += expansionType1Coeffs[supportIndicator[i]] * local_val;
      } else break;
    } 
	
    break;
  case true:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i ) {
      unsigned int this_point_level = colloc_pts[supportIndicator[i]].get_level();
      if ( this_point_level <= max_level ) {
	Real local_val = 1.0;
	const Real*     coeff2_i = expansionType2Coeffs[supportIndicator[i]]; // column vector
	const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
	const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
	const Real2DArray& this_point_support = supports[supportIndicator[i]];
	RealVector terms(numVars+1);
	terms[0] = expansionType1Coeffs[supportIndicator[i]];
	for ( unsigned int j = 0; j < numVars; ++j) {
	  terms[j+1] = coeff2_i[j];
	}
	for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	  const unsigned int this_dim_level = level_index[dim_idx][0];
	  const unsigned int this_dim_index = level_index[dim_idx][1];
	  if ( this_dim_level == 1 ) {
	    terms[0] *= 1;
	    if (x[dim_idx] == point[dim_idx]) {
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 0 : 1;
	      }
	    } else if (x[dim_idx] < point[dim_idx]) {
	      Real t = (x[dim_idx] - this_point_support[0][dim_idx])/
		(point[dim_idx] - this_point_support[0][dim_idx]);
	      Real dx_by_dt = 
		point[dim_idx] - this_point_support[0][dim_idx];
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		  t*t*(t-1)*dx_by_dt : 1;
	      }
	    } else {
	      const Real t = (x[dim_idx] - point[dim_idx])/
		(this_point_support[1][dim_idx] - point[dim_idx]);
	      const Real dx_by_dt = 
		this_point_support[1][dim_idx] - point[dim_idx];
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		  dx_by_dt*t*(1-t)*(1-t) : 1;
	      }
	    }
	  } else {
	    if (x[dim_idx] == point[dim_idx]) {
	      terms[0] *= 1;
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 0 : 1;
	      }
	    } else if (x[dim_idx] < point[dim_idx]) {
	      const Real t = (x[dim_idx] - this_point_support[0][dim_idx])/
		(point[dim_idx] - this_point_support[0][dim_idx]);
	      const Real dx_by_dt = 
		(point[dim_idx] - this_point_support[0][dim_idx]);
	      terms[0] *= t*t*(3-2*t);
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		  t*t*(t-1)*dx_by_dt : t*t*(3-2*t);
	      }
	    } else {
	      const Real t = (x[dim_idx] - point[dim_idx])/
		(this_point_support[1][dim_idx] - point[dim_idx]);
	      const Real dx_by_dt = 
		(this_point_support[1][dim_idx] - point[dim_idx]);
	      terms[0] *= (1+2*t)*(1-t)*(1-t);
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		terms[dim_idx2+1] *= (dim_idx == dim_idx2) ? 
		  t*(1-t)*(1-t)*dx_by_dt : (1+2*t)*(1-t)*(1-t);
	      }
	    }
	  }
	}
	for ( unsigned int accumulator = 0; accumulator < terms.length(); ++accumulator ) {
	  approx_val += terms[accumulator];
	}
	//approx_val += std::accumulate(terms.RealMatrix::operator[](0),
	//				 terms.RealMatrix::operator[](0)+terms.length(),
	//				 0);
      } else break;
    }
    break;
  }
  return approx_val;
}

const RealVector& LRDHierarchInterpPolyApproximation::
gradient_basis_variables(const RealVector& x, unsigned int max_level)
{
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in "
	  << "LRDHierarchInterpPolyApproximation::get_value()" << std::endl;
    abort_handler(-1);
  }
  approxGradient.size(numVars);
  for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx )
    approxGradient[dim_idx] = 0; 
    
  int num_colloc_points = driverRep->grid_size();

  //Grab the supports of the hierarchical basis functions.
  LocalRefinableDriver* lr_driver = 
    static_cast<LocalRefinableDriver*>(driverRep);
  const std::vector<Real2DArray>& supports = lr_driver->get_supports();

  //Determine which basis functions share their support with x.
  const IntArray& supportIndicator = in_support_of(x);
  const std::vector<CollocationPoint>& colloc_pts = lr_driver->get_collocation_points();
  switch ( basisConfigOptions.useDerivs ) {
  case false:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      unsigned int this_point_level = colloc_pts[supportIndicator[i]].get_level();
      if ( this_point_level <= max_level ) {
	const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
	const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
	const Real2DArray& this_point_support = supports[supportIndicator[i]];

	for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	  Real localGrad_i = 1;
	  const unsigned int this_dim_level = level_index[dim_idx][0];
	  const unsigned int this_dim_index = level_index[dim_idx][1];
	  if ( this_dim_level == 1 ) {
	    localGrad_i *= 0;  //constant case
	  } else { //Linear hat functions
	    if ( (x[dim_idx] == point[dim_idx]) || 
		 (x[dim_idx] == this_point_support[0][dim_idx]) || 
		 (x[dim_idx] == this_point_support[1][dim_idx]) ) {
	      localGrad_i *= 0;
	    } else if (x[dim_idx] < point[dim_idx]) {
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		localGrad_i *= (dim_idx == dim_idx2) ? 
		  1 / ( point[dim_idx] - this_point_support[0][dim_idx] ) :
		  (x[dim_idx2] - this_point_support[0][dim_idx]) / 
		  ( point[dim_idx2] - this_point_support[0][dim_idx2] );
	      }
	    } else {
	      for ( unsigned int dim_idx2 = 0; dim_idx2 < numVars; ++dim_idx2 ) {
		localGrad_i *= (dim_idx == dim_idx2) ? 
		  -1 / ( this_point_support[1][dim_idx2] - point[dim_idx2] ) :
		  (point[dim_idx2] - x[dim_idx2])/( this_point_support[1][dim_idx2] - point[dim_idx2] );
	      }
	    }   
	  }
	  approxGradient[dim_idx] = approxGradient[dim_idx] + 
	    expansionType1Coeffs[supportIndicator[i]] * localGrad_i;
	}
      } else break;
    }
    break;
  case true:
    // Sum over only those elements whose support contains x
    for ( unsigned int i = 0; i<supportIndicator.size(); ++i) {
      unsigned int this_point_level = colloc_pts[supportIndicator[i]].get_level();
      if ( this_point_level <= max_level ) {
	const Real*     coeff2_i = expansionType2Coeffs[supportIndicator[i]];
	const RealVector& point = colloc_pts[supportIndicator[i]].get_point();
	const Int2DArray& level_index = colloc_pts[supportIndicator[i]].get_level_index();
	const Real2DArray& this_point_support = supports[supportIndicator[i]];

	// These nested loops are filthy.
	// dim_idx keeps track of which gradient entry we're filling.
	// dim_idx2 keeps track of which variable we're thinking about.
	// dim_idx3 tess us which term we're working on. 
	for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	  RealVector terms(numVars+1);
	  terms[0] = expansionType1Coeffs[supportIndicator[i]];
	  for ( unsigned int j = 0; j < numVars; ++j ) {
	    terms[j+1] = coeff2_i[j];
	  }
	  for ( unsigned int dim_idx2 = 0; dim_idx2< numVars; ++dim_idx2 ) { 
	    const unsigned int this_dim_level = level_index[dim_idx2][0];
	    const unsigned int this_dim_index = level_index[dim_idx2][1];
	    if ( this_dim_level == 1 ) {
	      terms[0] *= (dim_idx == dim_idx2) ? 0 : 1;  
	      if ( x[dim_idx2] == point[dim_idx2] ) {
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  1 : 0;  // Take derivative of basis else don't
		  } else {  //Use type 1.
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 0 : 1; // Take derivative of basis else don't
		  }
		}
	      } else if ( x[dim_idx2] < point[dim_idx2] ) {
		const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		  (point[dim_idx2] - this_point_support[0][dim_idx2]);
		const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
		const Real dt_by_dx = 1/dx_by_dt;
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		      (3*t*t-2*t) : 
		      dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		  } else { //use type 1
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		      0 : 1; // Take derivative of basis else don't
		  }
		}
	      } else {
		const Real t = (x[dim_idx2] - point[dim_idx2])/
		  (this_point_support[1][dim_idx2] - point[dim_idx2]);
		const Real dx_by_dt = this_point_support[1][dim_idx2] - point[dim_idx2];
		const Real dt_by_dx = 1/dx_by_dt;
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		      (3*t*t - 4*t + 1) : 
		      dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		  } else { //use type 1
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		      0 : 1; // Take derivative of basis else don't
		  }
		}
	      }
	    } else {  
	      if ( x[dim_idx2] == point[dim_idx2] ) {
		terms[0] *= (dim_idx == dim_idx2) ? 0 : 1;
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  1 : 0;  // Take derivative of basis else don't
		  } else {  //Use type 1.
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 0 : 1; // Take derivative of basis else don't
		  }
		}
	      } else if ( x[dim_idx2] < point[dim_idx2] ) {
		const Real t = (x[dim_idx2] - this_point_support[0][dim_idx2])/
		  (point[dim_idx2] - this_point_support[0][dim_idx2]);
		const Real dx_by_dt = point[dim_idx2] - this_point_support[0][dim_idx2];
		const Real dt_by_dx = 1/dx_by_dt;
		terms[0] *= (dim_idx == dim_idx2) ?  dt_by_dx*(-6*t*t + 6*t) : -2*t*t*t + 3*t*t;
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		      (3*t*t-2*t) : 
		      dx_by_dt*(t*t*t-t*t);  // Take derivative of basis else don't
		  } else { //use type 1
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		      dt_by_dx*(-6*t*t+6*t) : t*t*(3-2*t); // Take derivative of basis else don't
		  }
		}
	      } else {
		const Real t = (x[dim_idx2] - point[dim_idx2])/
		  (this_point_support[1][dim_idx2] - point[dim_idx2]);
		const Real dx_by_dt = this_point_support[1][dim_idx2] - point[dim_idx2];
		const Real dt_by_dx = 1/dx_by_dt;
		terms[0] *= (dim_idx == dim_idx2) ?  dt_by_dx*(6*t*t - 6*t) : 2*t*t*t - 3*t*t + 1;
		for ( unsigned int dim_idx3 = 0; dim_idx3 < numVars; ++dim_idx3 ) {
		  if ( dim_idx2 == dim_idx3 ) { //Use type 2 function
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ?  
		      (3*t*t - 4*t + 1) : 
		      dx_by_dt*t*(1-t)*(1-t);  // Take derivative of basis else don't
		  } else { //use type 1
		    terms[dim_idx3+1] *= (dim_idx2 == dim_idx) ? 
		      dt_by_dx*(6*t*t-6*t) : (1+2*t)*(1-t)*(1-t); // Take derivative of basis else don't
		  }
		}
	      }
	    }
	  }
	  for ( unsigned int accumulator = 0; accumulator < terms.length(); ++accumulator) {
	    approxGradient[dim_idx] += terms[accumulator];
	  }
	  //approxGradient[dim_idx] += std::accumulate(terms.RealMatrix::operator[](0),
	  //					       terms.RealMatrix::operator[](0) + terms.length(),
	  //					       0);
	}
      } else break;
    }
    break;
  }
  return approxGradient;
}


void LRDHierarchInterpPolyApproximation::
compute_numerical_response_moments(size_t num_moments)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_response_moments()" << std::endl;
    abort_handler(-1);
  }

  // stub
}


void LRDHierarchInterpPolyApproximation::
compute_numerical_expansion_moments(size_t num_moments)
{
  // Error check for required data
  if (!expConfigOptions.expansionCoeffFlag) {
    PCerr << "Error: expansion coefficients not defined in InterpPoly"
	  << "Approximation::compute_numerical_expansion_moments()"<< std::endl;
    abort_handler(-1);
  }
  if (expansionMoments.length() != num_moments)
    expansionMoments.sizeUninitialized(num_moments);

  // stub
}


void LRDHierarchInterpPolyApproximation::compute_total_sobol_indices()
{
  // stub
}


void LRDHierarchInterpPolyApproximation::
member_coefficients_weights(int set_value, const UShortArray& quad_order,
			    const UShortArray& lev_index,
			    const UShort2DArray& key,
			    const SizetArray& colloc_index,
			    RealVector& member_coeffs, RealVector& member_wts)
{
  // stub
}

}
