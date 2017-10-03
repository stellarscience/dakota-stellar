// $Id: HOPSPACK_EvalCounter.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_EvalCounter.hpp $ 

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
  \file HOPSPACK_EvalCounter.hpp
  \brief Class declaration for HOPSPACK::EvalCounter
*/
#ifndef HOPSPACK_PERFCOUNTER_HPP
#define HOPSPACK_PERFCOUNTER_HPP

#include <map>
#include "HOPSPACK_common.hpp"

namespace HOPSPACK
{

/*! \brief Counts the different types of evaluations, using messages
  returned by the user's evaluation routine.

  This object is used by the Conveyor for counting the number of
  cached points and evaluated points, on a per worker basis (worker IDs
  are typically the same as MPI node IDs).

  User evaluation can provide a message string, especially for errors.
  The counter object tracks how many times each error string is returned,
  on a per worker basis.

  In summary, the counter tracks the following.
  <ul>
    <li> The total number of successful cache look-ups.
    <li> The total number of evaluations.
    <li> The number of evaluations per worker.
    <li> The total number of times each user message appeared.
    <li> The number of times each user message appeared, per worker.
  </ul>
*/
class EvalCounter
{
public:

  /*! Constructor */
  EvalCounter();

  /*! Destructor */
  ~EvalCounter();


  //! Returns the total number of evaluations
  int getNumEvaluations() const;

  //! Returns the total number of cached evaluations
  int getNumCached() const;

  //! Copies a string with the current global counts into sResult.
  void  getCountString (string &  sResult) const;


  //! Increment the number of cached function values
  void incrementCached();

  //! Increment the number of cached function values
  void incrementPendingCached();

  //! Increment the number of evaluations 
  /*!
    \param workerId - The id of the worker that performed the evaluation. 
    \param msg - The user-defined message about the function evaluation.
   */
  void incrementEvaluated(int workerId, const string& msg);

  //! Print the current counts
  /*!
   *  \param bDisplayDetails  Details include evaluation breakdowns by
   *                          message type and worker ID.
   */
  void print (const bool  bDisplayDetails) const;

private:

  //! Used to count the number of each type of evaluation
  typedef map<string, int> MsgCnt;

  //! Used to keep per worker counts
  typedef map<int, MsgCnt> WkrMsgCnt;

  //! Iterator for MsgCnt
  typedef MsgCnt::iterator MsgCntIterator;
  //! Iterator for WkrMsgCnt
  typedef WkrMsgCnt::iterator WkrMsgCntIterator;
  //! Const iterator for MsgCnt
  typedef MsgCnt::const_iterator MsgCntConstIterator;
  //! Const iterator for WkrMsgCnt
  typedef WkrMsgCnt::const_iterator WkrMsgCntConstIterator;

  
  //! Increments the number of times that the given msg has appeared in mc.
  void incrementMsgCnt(MsgCnt& mc, const string& msg);

  //! The number of each user-defined type of evaluation
  MsgCnt msgCnt;

  //! The number of each user-defined type of evaluation per worker
  WkrMsgCnt wkrMsgCnt;

  //! The total number of cached function values
  int nCached;

  //! The total number of cached function values from pending list.
  int nPendingCached;

  //! The total number of evaluations
  int nEvaluated;

  //! String used by getCountString().
  string sCountString;
};

}

#endif
