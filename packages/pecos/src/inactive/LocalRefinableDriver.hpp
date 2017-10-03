/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:	 LocalRefinableDriver
//- Description: 
//- Owner:       Mike Eldred
//- Revised by:  
//- Version:

#ifndef LOCAL_REFINABLE_DRIVER_HPP
#define LOCAL_REFINABLE_DRIVER_HPP

#include "IntegrationDriver.hpp"


namespace Pecos {

///Utility class for RefinableGrids
class CollocationPoint
{
public:
    
  ///Default Constructor
  CollocationPoint();

  ///Destructor
  virtual ~CollocationPoint();

  ///Standard Constructor
  CollocationPoint( const RealVector& point_,
		    const Int2DArray& levelIndex_= Int2DArray(0));
    
  ///Get Point location.
  const RealVector& get_point() const;
    
  ///Get index of point.
  const Int2DArray& get_level_index() const;
    
  ///Get level of point.
  const unsigned int get_level() const;
  
  ///Equality operator.
  bool operator==(const CollocationPoint& col_point) const;
    
  ///Less than operator.
  ///Ordering is lexiographic in the point's coordinates
  ///Used for sorting the grid prior to weeding out duplicates.
  bool operator<(const CollocationPoint& col_point) const;

private:
    
  ///Point coordinates
  RealVector point;

  ///Point index
  Int2DArray levelIndex;

  ///Point level
  mutable unsigned int level;
};


inline CollocationPoint::CollocationPoint()
{}

inline CollocationPoint::~CollocationPoint()
{}

inline CollocationPoint::CollocationPoint( const RealVector& point_,
					   const Int2DArray& levelIndex_):
  point(point_), levelIndex(levelIndex_), level(0)
{}

inline const RealVector& CollocationPoint::get_point() const
{ return point; }

inline const Int2DArray& CollocationPoint::get_level_index() const
{ return levelIndex; }

inline const unsigned int CollocationPoint::get_level() const
{
  if (!level) {
    level = 1;
    for ( Int2DArray::const_iterator it = levelIndex.begin();
	  it != levelIndex.end(); ++it)
      level += (*it)[0] - 1;
  }	
  return level;
}

inline bool CollocationPoint::
operator==(const CollocationPoint& col_point) const
{ return ( point==col_point.get_point() ) ? true : false; }
  
inline bool CollocationPoint::operator<(const CollocationPoint& col_point) const
{
  const RealVector& point2 = col_point.get_point();
  for ( size_t idx = 0; idx < point.length(); idx++ ) {
    if ( point[idx] < point2[idx] ) return true;
    if ( point[idx] > point2[idx] ) return false;
  }
  return false;
}


typedef std::vector<CollocationPoint>::iterator ColPtIterator;
typedef std::vector<CollocationPoint>::const_iterator const_ColPtIterator;


/// Derived nondeterministic class that generates N-dimensional
/// hierarchical adaptive sparse grids for numerical evaluation of expectation
/// integrals over independent standard random variables.

class LocalRefinableDriver: public IntegrationDriver
{
public:
    
  ///Default constructor
  LocalRefinableDriver();
    
  ///Destructor
  virtual ~LocalRefinableDriver();

  ///Unpacks the grid into a RealMatrix
  virtual void compute_grid(RealMatrix& variable_sets);

  ///Returns the grid size
  virtual int grid_size() const;

  ///Returns the size of the highest level
  virtual int highest_level_size() const;

  ///Returns the current grid level.
  virtual int get_current_level() const;

  ///Assembles the sparse grid.
  virtual void initialize_grid(const RealArray& lower_bounds,
			       const RealArray& upper_bounds,
			       const unsigned int starting_level = 1,
			       const short poly_type = PIECEWISE_LINEAR_INTERP,
			       const bool use_derivs = false);

  ///Refines every point at the highest level.  Calls refine_locally().
  virtual void refine_globally();

  ///Refine the points at the highest level indicated by refinement_selector.
  virtual void refine_locally(const BoolDeque& refinement_selector);

  ///TODO: Deletes the points at the highest level and decrements the level counter.
  virtual void undo_refinement();

  ///Get the set of collocation points
  virtual const std::vector<CollocationPoint>& get_collocation_points();

  ///Get the set of points at the highest level
  virtual const_ColPtIterator get_highest_level_start() const;

  ///Get the supports of the basis functions
  virtual const std::vector<Real2DArray>& get_supports() const;

  /// Prints the hierarchical grid.
  friend std::ostream& 
  operator<<(std::ostream& ostr, const LocalRefinableDriver& pointSet);

protected:

  /// 2 by numVars array containing the upper and lower interpolation bounds.
  Real2DArray bounds;
  /// The centroid of the interpolation domain
  RealVector midpoint;

  /// The current adaptive sparse grid level.
  unsigned int current_level;

  /// The number of points at the highest level.
  unsigned int current_level_size;

  /// The total number of collocation points.
  unsigned int num_interp_points;
    
  /// The set of collocation points.
  std::vector<CollocationPoint> interp_points;

  /// Iterator pointing to the first point in the highest level
  ColPtIterator highest_level_start;

  /// The supports of the basis functions.
  std::vector<Real2DArray> supports;

  short poly_type;

  bool isInitialized;

  /// the set of type1 weights (for integration of value interpolants)
  /// associated with each point in the grid
  RealVector type1WeightSets;
  /// the set of type2 weights (for integration of gradient interpolants)
  /// for each derivative component and for each point in the grid
  RealMatrix type2WeightSets;

  ///Computes the quadrature weights for the points at the highest level.
  virtual void set_highest_level_weights();

  ///Computes the support of the basis functions at the highest level.
  virtual void set_highest_level_supports();

  ///Associates a point with a given level/index pair.
  virtual void get_interp_point(const Int2DArray& index,
				RealVector& point) const;

private:

};

inline LocalRefinableDriver::LocalRefinableDriver():
  IntegrationDriver(BaseConstructor()), bounds(), midpoint(), current_level(0),
  current_level_size(0), num_interp_points(0), interp_points(),
  highest_level_start(interp_points.begin()), supports(), poly_type(0),
  isInitialized(false)
{}

inline LocalRefinableDriver::~LocalRefinableDriver()
{}

inline int LocalRefinableDriver::grid_size() const
{ return num_interp_points; }

inline int LocalRefinableDriver::highest_level_size() const
{ return current_level_size; }

inline int LocalRefinableDriver::get_current_level() const
{ return current_level; }

inline const std::vector<CollocationPoint>& LocalRefinableDriver::
get_collocation_points()
{ return interp_points; }

inline const_ColPtIterator LocalRefinableDriver::get_highest_level_start() const
{ return highest_level_start; }

inline const std::vector<Real2DArray>& LocalRefinableDriver::
get_supports() const
{ return supports; }

}

#endif
