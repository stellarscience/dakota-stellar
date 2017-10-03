// $Id: HOPSPACK_ConveyorList.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_ConveyorList.hpp $ 

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
  \file HOPSPACK_ConveyorList.hpp
  \brief Class declaration for HOPSPACK::ConveyorList.
*/
#ifndef HOPSPACK_CONVEYORLIST_HPP
#define HOPSPACK_CONVEYORLIST_HPP

#include <list>

#include "HOPSPACK_CacheManager.hpp"
#include "HOPSPACK_DataPoint.hpp"

namespace HOPSPACK
{

//! Manages a list of DataPoint instances for a Conveyor.
/*!
  The purpose of this class is to support communication between the Conveyor
  and Mediator.  The list of points is ordered with the highest priority
  point at the front.  The entire list is assigned a priority from its
  source citizen.
*/
class ConveyorList
{
public:

  //! Constructor 
  ConveyorList (void);

  //! Destructor: deletes the contents of the current list 
  ~ConveyorList();
  
  //! Push a point onto the front of the list.
  /*!
   *  The DataPoint is pushed onto the front of the list.
   *  Together, push() and pop() provide FIFO behavior.
   *  The list takes ownership for the given pointer. 
  */
  void push (DataPoint* trialPointPtr);

  //! Pop a point from the end of the list.
  /*!
   *  If the list is empty, return NULL. Otherwise, pop the DataPoint off
   *  the end of the list.
   *  Together, push() and pop() provide FIFO behavior.
   *  Ownership of the pointer is passed to the calling object.
  */
  DataPoint* pop();

  //! Pop the point with the given tag.
  /*!
   *  If the point with the given tag is not in the list, return NULL.
   *  Otherwise, pop the DataPoint with the given tag off the list.
   *  Ownership of the pointer is passed to the calling object.
  */
  DataPoint* pop (int tag);

  //! Pop the next point if found in the cache.
  /*!
   *  Search from the end of the list and return the first point found
   *  in the cache.  Return NULL if none found.
   *  Ownership of the pointer is passed to the calling object.
  */
  DataPoint *  popNextCached (CacheManager * const  pCache);

  //! Return size of the list
  int size() const;

  //! Return true if the list is empty 
  bool isEmpty() const;

  //! Return a const reference to the list of points
  const list< DataPoint * > &  getPointList() const;

  //! Return a nonconst reference to the list of points
  list< DataPoint * > &  getMutablePointList();

  //! Returns iterator to the first element of stored list 
  list< DataPoint * >::const_iterator begin() const;  

  //! The end() function returns an iterator just past the end of the list.
  list< DataPoint * >::const_iterator end() const;

  //! Set the priority of the list (from the citizen).
  void  setPriority (const int  nPriority);

  //! Get the priority of the list (from the citizen).
  int  getPriority (void) const;

  //! Copies the tags of all points into tagList
  void getTagList(vector<int>& tagList) const;

  //! Prune the list except for the most recently added n points.    
  /*!
   *  Remove and delete points from the end of the list until only n remain.
   *  The most recently added points are at the front of the list.
   */
  void prune(int n = 0);

  //! Return true and assign the tag if w "equals" a point in this list.
  /*!
   *  @param[in] w     Point whose location is to be tested.
   *  @param[out] tag  Is assigned the tag value of the "equal" point,
   *                   or -1 if none found.
   *  @return true     If a point is found in this list that is "equal" to w,
   *                   as determined by ScaledComparison.
   *
   *  The ScaledComparison test necessarily matches the comparison test used
   *  by the Cache.
  */
  bool contains(const DataPoint *w, int& tag) const;

  //! Prints contents of list 
  void print(const string label) const;
    
private:

  //! By design, there is no copy constructor.
  ConveyorList (const ConveyorList &);
  //! By design, there is no assignment operator.
  ConveyorList & operator= (const ConveyorList &);


  //! A list of pointers to trial points
  typedef list< DataPoint * > PointListType;

  //! List of points.
  PointListType  _cPtList;

  //! Priority of the entire list (from the citizen).
  int  _nPriority;
};


}
#endif
