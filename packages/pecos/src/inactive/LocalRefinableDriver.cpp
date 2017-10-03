/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        UniformRefinableGrid
//- Description:  Class implementing a refinable point set using uniformly
//	          spaced points on [a_0,b_0]\times ... \times [a_d-1,b_d-1]
//-               
//- Owner:        Christopher Miller: University of Maryland at College Park
//- Contact:      cmiller@math.umd.edu

#include "LocalRefinableDriver.hpp"


namespace Pecos {

void LocalRefinableDriver::
initialize_grid(const RealArray& lower_bounds,const RealArray& upper_bounds,
		const unsigned int starting_level,  const short poly_type_,
		const bool use_derivs)
{
  if ( (poly_type_ != PIECEWISE_LINEAR_INTERP) && 
       (poly_type_ != PIECEWISE_CUBIC_INTERP) ){
    PCerr << "Polynomial type not supported by this driver"
	  << "defaulting to PIECEWISE_LINEAR_INTERP" << std::endl;
    poly_type = PIECEWISE_LINEAR_INTERP;
  }
  else
    poly_type = poly_type_;

  //Check that lower an upper bounds have same numVars.
  if ( lower_bounds.size() != upper_bounds.size() ) {
    PCerr << "Error: lower_bounds and upper_bounds must be the same size."
	  << " Got lower_bounds.size() == " << lower_bounds.size()
	  << " and upper_bounds.size() == " << upper_bounds.size()
	  << "." << std::endl;
    throw ( std::invalid_argument("") );
  }
  if (use_derivs) computeType2Weights = true;
  numVars = lower_bounds.size();
  bounds.resize(2);
  bounds[0] = lower_bounds;
  bounds[1] = upper_bounds;
  midpoint.size(numVars);

  for ( unsigned int i = 0; i < numVars; ++i ) {
    if ( lower_bounds[i] >= upper_bounds[i] ) {
      PCerr << "Error: Expected lower_bounds[i] < upper_bounds[i]."
	    << " Got lower_bounds[" << i << "] = " << lower_bounds[i]
	    << " and upper_bounds[" << i << "] = " << upper_bounds[i]
	    << "." << std::endl;
      throw ( std::invalid_argument("") );
    }
    else
      midpoint[i] = lower_bounds[i] + (upper_bounds[i] - lower_bounds[i]) / 2.0;
  }

  IntArray oneDLevelIndex(2,0);
  oneDLevelIndex[0] = 1;
  interp_points.resize(1);
    
  interp_points[0] = 
    CollocationPoint( midpoint,Int2DArray(numVars,oneDLevelIndex) );
  highest_level_start = interp_points.begin();

  supports.resize(1);
  supports[0] = bounds;

  current_level = 1;
  current_level_size = 1;
  num_interp_points = 1;

  set_highest_level_weights();
  isInitialized = true;
 
  while ( current_level != starting_level )
    refine_globally();
}


void LocalRefinableDriver::refine_globally()
{ refine_locally( BoolDeque(current_level_size, true) ); }


void LocalRefinableDriver::
refine_locally(const BoolDeque& refinement_selector)
{
  if ( !isInitialized ){
    PCerr << "Need to initialize grid before performing refinements."
	  << std::endl;
    return;
  } 

  if ( refinement_selector.size() != current_level_size ) {
    PCerr << "Error: refinement_selector needs to be the "
          << "same size as current_level_size. Got refinement_selector"
	  << ".size() == " 
          << refinement_selector.size() << " and current_level_size == "
          << current_level_size << "." << std::endl;
    throw ( std::invalid_argument("") );
  }

  //The integer index of the current interp point
  unsigned int masterIndex;
  unsigned int this_dim_level, this_dim_index;
  Int2DArray childPointIndex1, childPointIndex2;

  RealVector newPoint(numVars);

  size_t integer_index = 
    std::distance(interp_points.begin(),interp_points.end());
  for ( unsigned int i = 0; i < current_level_size; ++i ) {
    if ( refinement_selector[i] ) { //refine this point?
      masterIndex = num_interp_points - current_level_size + i;
      childPointIndex1 = interp_points[masterIndex].get_level_index();
      childPointIndex2 = childPointIndex1;
      // loop over numVars
      for ( unsigned int dim_idx = 0; dim_idx < numVars; ++dim_idx ) {
	this_dim_level = childPointIndex1[dim_idx].at(0);
	this_dim_index = childPointIndex1[dim_idx].at(1);
	  
	if ( this_dim_level == 2 ) {
	  childPointIndex1[dim_idx][0] = this_dim_level + 1;
	  childPointIndex1[dim_idx][1] = this_dim_index;
	  this->get_interp_point(childPointIndex1,newPoint);
	  interp_points.
	    push_back(CollocationPoint(newPoint,childPointIndex1));
	  --childPointIndex1[dim_idx][0];
	} else { //General case
	  childPointIndex1[dim_idx][0] = this_dim_level + 1;
	  childPointIndex1[dim_idx][1] = 2*this_dim_index;
	  childPointIndex2[dim_idx][0] = this_dim_level + 1;
	  childPointIndex2[dim_idx][1] = 2*this_dim_index + 1;
	  this->get_interp_point(childPointIndex1,newPoint);
	  interp_points.
	    push_back(CollocationPoint(newPoint,childPointIndex1));
	  this->get_interp_point(childPointIndex2,newPoint);
	  interp_points.
	    push_back(CollocationPoint(newPoint,childPointIndex2));
	  --childPointIndex1[dim_idx][0];
	  --childPointIndex2[dim_idx][0];
	  childPointIndex1[dim_idx][1] /= 2;
	  childPointIndex2[dim_idx][1] = 
	    ( childPointIndex2[dim_idx][1] - 1 ) / 2;
	}
      }
    }
  }
  // do some cleanup.  Remove duplicates and pack interp_points.
  highest_level_start = interp_points.begin();
  std::advance(highest_level_start,integer_index);
  std::sort(highest_level_start,interp_points.end());
  ColPtIterator it = 
    std::unique(highest_level_start,interp_points.end());
  interp_points.resize(std::distance(interp_points.begin(), it));
  current_level_size = 
    std::distance(highest_level_start,interp_points.end());
  num_interp_points += current_level_size;
  ++current_level;
    
  set_highest_level_supports();
  set_highest_level_weights();
}


void LocalRefinableDriver::
get_interp_point(const Int2DArray& index, RealVector& point) const
{
  if (point.length() != numVars) point.resize(numVars);
  assert( index.size() == numVars );
  unsigned int this_dim_level;
  unsigned int this_dim_index;

  for ( size_t idx = 0; idx < numVars; ++idx ) {
    this_dim_level = index[idx][0];
    this_dim_index = index[idx][1];
    if ( this_dim_level == 1 ) point[idx] = midpoint[idx];
    else if( this_dim_level == 2 ) {
      if( this_dim_index == 0 ) point[idx] = bounds[0][idx];
      else point[idx] = bounds[1][idx];
    }
    else {
      const Real offset = ( (bounds[0][idx] + bounds[1][idx]) / 2.0
	- bounds[0][idx]) / std::pow(2.0,(double)this_dim_level-2);
      point[idx] = bounds[0][idx] + (1 + 2*this_dim_index)*offset;
    }
  }
}


void LocalRefinableDriver::undo_refinement()
{
  std::cout << "Not yet implemented" << std::endl;
}


void LocalRefinableDriver::compute_grid(RealMatrix& var_sets)
{
  var_sets.shapeUninitialized(numVars, num_interp_points);
  for ( unsigned int i = 0; i < num_interp_points; ++i ) {
    //interp_points[i].get_point() is a SerialDenseVector whose [] operator is element access
    //To use std copy to pack the vector into the matrix get a pointer to the data by using
    //the parent's [] operator.
    Real const * fakeout =  interp_points[i].get_point().RealMatrix::operator[](0);
    //var_sets[i] = interp_points[i].get_point();
    std::copy(fakeout,fakeout+numVars,var_sets[i]);
    //std::copy(interp_points[i].get_point().begin(),
    //interp_points[i].get_point().end(),
    //	var_sets[i] );
  }
}


void LocalRefinableDriver::
set_highest_level_weights() 
{
  type1WeightSets.resize(type1WeightSets.length() + current_level_size);
  std::vector<Real2DArray>::iterator support_it = supports.begin();
  std::advance(support_it,std::distance(interp_points.begin(),highest_level_start));
  size_t int_idx = std::distance(interp_points.begin(),highest_level_start);
  for ( const_ColPtIterator it = highest_level_start; it != interp_points.end(); ++it, ++support_it, ++int_idx ){
    type1WeightSets[int_idx] = 1;
    const RealVector& point = it->get_point();
    const Int2DArray& level_index = it->get_level_index();
    for ( unsigned int dimIdx = 0; dimIdx < numVars; ++dimIdx ) {
      const Real right_support_width = (*support_it)[1][dimIdx] - point[dimIdx];
      const Real left_support_width = point[dimIdx] - (*support_it)[0][dimIdx];
      if (level_index[dimIdx][0] == 1) {
	type1WeightSets[int_idx] *= right_support_width + left_support_width;
      }
      else {
	type1WeightSets[int_idx] *= 0.5*(left_support_width + right_support_width);
      }
      if (computeType2Weights) {
	type2WeightSets.reshape(numVars, type2WeightSets.numCols() + current_level_size);
	for (unsigned int j = 0; j < current_level_size; ++j) {
	  //unsigned int integer_index = num_interp_points - current_level + j;
	  const RealVector& point = interp_points[int_idx].get_point();
	  const Real2DArray& support = supports[int_idx];
	  Real* t2_wt = type2WeightSets[int_idx];
	  for ( unsigned int k = 0; k<numVars; ++k ) {
	    Real& t2_wt_k = t2_wt[k]; t2_wt_k = 1.0;
	    for (unsigned int l = 0; l < numVars; ++l) {
	      t2_wt_k *= (l==k) ? (support[0][k] + support[1][k] -2*point[k]) / 12 :
		(point[k] - support[0][l])/2 + (support[1][k] - point[k])/2;
	    }
	  }
	}
      } 
    }
  }
}


std::ostream&
operator<<(std::ostream& ostr, const LocalRefinableDriver& pointSet)
{

  ostr << "Index" << "\t" << "[level/index]" << "\t" << "Point" << std::endl;
  ostr << "-----" << "\t" << "-------------" << "\t" << "-----" << std::endl;
  Int2DArray const* levelIndex;
  RealVector const* interpPoint;
  for ( unsigned int i = 0; i < pointSet.num_interp_points; ++i ) {
    levelIndex = &(pointSet.interp_points[i].get_level_index());
    interpPoint = &(pointSet.interp_points[i].get_point());
    ostr << i << "\t[" << (*levelIndex)[0][0] << "," << (*levelIndex)[0][1]
	 << "\t\t[" << (*interpPoint)[0] << std::endl;
    for ( unsigned int j = 1; j < pointSet.numVars-1; ++j ) {
      ostr << " " << "\t" << (*levelIndex)[j][0] << "," << (*levelIndex)[j][1]
	   << "\t\t" << (*interpPoint)[j] << std::endl;
    }
    ostr << " " << "\t" << (*levelIndex)[pointSet.numVars-1][0] 
	 << "," << (*levelIndex)[pointSet.numVars-1][1]
	 << "]\t\t" << (*interpPoint)[pointSet.numVars-1] 
	 << "]" << std::endl;
  }
  ostr << "Current level = " << pointSet.current_level << std::endl;
  ostr << "Num collocation points = " << pointSet.num_interp_points << std::endl;
  ostr << "Num points at highest level = " << pointSet.current_level_size << std::endl;
  return ostr;
}


void LocalRefinableDriver::set_highest_level_supports()
{
  supports.resize(supports.size() + current_level_size, Real2DArray(2,RealArray(numVars,0)));
  std::vector<Real2DArray>::iterator support_it = supports.begin();
  std::advance(support_it,std::distance(interp_points.begin(),highest_level_start));
  for ( const_ColPtIterator it = highest_level_start; it != interp_points.end(); ++it ){
    const Int2DArray& level_index = it->get_level_index();
    const RealVector& point = it->get_point();
    for ( unsigned int dimIdx = 0; dimIdx < numVars; ++dimIdx ){
      if ( level_index[dimIdx][0] == 1 ) {
	(*support_it)[0][dimIdx] = bounds[0][dimIdx];
	(*support_it)[1][dimIdx] = bounds[1][dimIdx];
      } else if ( level_index[dimIdx][0] == 2 ) {
	if ( level_index[dimIdx][1] == 0 ) {
	  (*support_it)[0][dimIdx] = bounds[0][dimIdx];
	  (*support_it)[1][dimIdx] = midpoint[dimIdx];
	} else {
	  (*support_it)[0][dimIdx] = midpoint[dimIdx];
	  (*support_it)[1][dimIdx] = bounds[1][dimIdx];
	}
      } else {
	Real offset = ( (bounds[0][dimIdx] + bounds[1][dimIdx]) / 2.0 - bounds[0][dimIdx]) 
	  / std::pow(2.0,(double)level_index[dimIdx][0]-2.0);
	(*support_it)[0][dimIdx] = point[dimIdx] - offset;
	(*support_it)[1][dimIdx] = point[dimIdx] + offset;
      }
    }
    ++support_it;
  }
}

}
