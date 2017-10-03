// $Id: HOPSPACK_GssList.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssList.hpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This file is part of HOPSPACK.
//
// HOPSPACK is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library.  If not, see http://www.gnu.org/licenses/.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov)
//                 or Todd Plantenga (tplante@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*!
  @file HOPSPACK_GssList.hpp
  @brief Class declaration for HOPSPACK::GssList.
*/
#ifndef HOPSPACK_GSSLIST_HPP
#define HOPSPACK_GSSLIST_HPP

#include <list>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_GssPoint.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Manages a list of GssPoint instances for a CitizenGSS.
/*!
  This class is similar to ConveyorList, but does not have enough in common
  to use inheritance.  A GssList can copy DataPoint instances to and from
  a ConveyorList, making the necessary type conversion.  GssList also
  supports the notion of ordering points by "best" objective value.
*/
//----------------------------------------------------------------------
class GssList
{
  public:

    //! Constructor.
    GssList (void);

    //! Destructor.
    ~GssList (void);

    //! Set a step length to use in copyFrom().
    void  setDefaultStepLength (const double  dStepLength);

    //! Return true if the list is empty.
    bool  isEmpty (void) const;
  
    //! Push the given point onto the front of the list.
    /*!
     *  @param[in] pNewPoint  Pointer to new point handed to the list.
     */
    void  push (GssPoint *  pNewPoint);

    //! Return a const pointer to the best point in the list, or NULL if empty.
    /*!
     *  Only the single objective value returned by getBestF() is considered;
     *  in particular, all points are assumed to be linearly feasible.
     *
     *  As a side effect, the best point is moved to the end of the list,
     *  making it more efficient to follow with a call to popBest().
     *
     *  @return  Best point, or NULL if the list is empty.
     */
    const GssPoint *  findBest (void);

    //! Pop the next point from the list.
    /*!
     *  @return  Point from the end of the list, or NULL if the list is empty.
     */
    GssPoint *  pop (void);

    //! Pop the best point off the list and return it.
    /*!
     *  Only the single objective value returned by getBestF() is considered;
     *  in particular, all points are assumed to be feasible.
     *
     *  @return  Best point, or NULL if the list is empty.
     */
    GssPoint *  popBest (void);

    //! Copy points into the list, converting from type DataPoint to GssPoint.
    /*!
     *  All points from cInputList are copied, preserving the order.
     *  Those from another citizen are created with an empty parent tag
     *  and the default step length.
     *
     *  @param[in] cInputList  List of all points.
     *  @param[in] cPenalty    Penalty function for nonlinear constraints.
     *  @param[in] cTagList    List of tags belonging to this citizen.
     */
    void  copyFrom (const list< DataPoint * > &  cInputList,
                    const NonlConstrPenalty   &  cPenalty,
                    const list< int >         &  cTagList);

    //! Copy points out of the list to a target list.
    /*!
     *  Points are copied as subclass GssPoint and added as base class DataPoint.
     *  Points are added in order to the front of cOutputList.
     *  The conveyor will take points from the end of the list, so all points
     *  in this list come behind any existing points in cOutputList.
     *
     *  @param[out] cOutputList  List where new points will be addded.
     */
    void  copyTo (list< DataPoint * > &  cOutputList);

    //! Insert (copy) all points from the source to the front of this list.
    /*!
     *  @param[in] cSource   List of points to be inserted.
     */
    void  insertFromList (const GssList &  cSource);

    //! Copy all point tags in the list to cTags.
    /*!
     *  @param[out] cTags  Modify this by adding tags from all points belonging
     *                     to the current instance.
     */
    void  copyTags (vector< int > &  cTags) const;

    //! Erase all points in the list.
    void  clearList (void);

    //! Delete all points in the list except the most recently added n points.
    /*!
     *  The most recently added points are at the front of the list.
     *  @param[in] n  Number of points at front to preserve.
     */
    void prune(const int n = 0);

    //! Prints contents of the list .
    void print (const string label) const;
    
  private:

    //! By design, there is no copy constructor.
    GssList (const GssList &);
    //! By design, there is no assignment operator.
    GssList & operator= (const GssList &);

    //! Move the best point to the end of the list.  The list cannot be empty.
    void  moveBestToEndOfList_ (void);


    //! A list of pointers to trial points.
    typedef list< GssPoint * > GssPointListType;

    //! List of points.
    GssPointListType  _cGssPtList;

    //! Default step length.
    double  _dDefaultStepLength;
};


}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_GSSLIST_HPP
