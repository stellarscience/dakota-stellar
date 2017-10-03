// $Id: HOPSPACK_CachePoint.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_CachePoint.hpp $ 

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
  \file HOPSPACK_CachePoint.hpp
  \brief Class declaration for HOPSPACK::CachePoint
  \author H. Alton Patrick, Summer 2000
*/

#ifndef HOPSPACK_CACHEPOINT_HPP
#define HOPSPACK_CACHEPOINT_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK 
{

//! Contains a cached point, comparable based on scaling and tolerance.
/*!
  Contains a cached point that allows comparisons between points for
  storage in a splay tree. See CacheManager.
*/
class CachePoint
{

public:

  //! Constructor
  CachePoint (void);

  //! Constructor
  /*! \note Saves a reference to x_in, but DOES NOT COPY IT. In other
    words, this is only a shallow copy, not a deep copy. We are
    exploiting the way that the splay tree is implemented in order to
    minimize copying of vectors. We know that any CachePoint that is
    created by this constructor will later be copied into a new
    CachePoint (and then we will actually do a deep copy).
  */
  CachePoint (const Vector& x_in);

  //! Constructor
  /*! \note Saves a reference to x_in, but DOES NOT COPY IT. In other
    words, this is only a shallow copy, not a deep copy. We are
    exploiting the way that the splay tree is implemented in order to
    minimize copying of vectors. We know that any CachePoint that is
    created by this constructor will later be copied into a new
    CachePoint (and then we will actually do a deep copy).
  */
  CachePoint (const Vector& x_in,
              const Vector& f_in,
              const Vector& cEqs_in,
              const Vector& cIneqs_in);
    
  //! Copy constructor
  /*! \note Deep copy. */
  CachePoint(const CachePoint& source);

  //! Destructor 
  ~CachePoint();

  //! Copy the relevant data from another cached point
  void copyData(const CachePoint& source);

  //! Extract data
  const Vector& getF();

  //! Extract data
  const Vector& getEqs();

  //! Extract data
  const Vector& getIneqs();

  /*! \brief Compare two CachePoints based on a lexicographic ordering of their
    coordinates.

    Components are compared in order.  The first component of "pt" that is
    clearly greater than or less than its counterpart in "this" determines
    whether the entire point is greater than or less than.
    The definition of "clearly" different components is controlled by
    the scaling and tolerance in ScaledComparison.
  */
  bool operator>(const CachePoint& pt) const;

  //! Reverse of operator>
  bool operator<(const CachePoint& pt) const;

  /*! \brief Return true if the two CachePoints are not "equal".

    Two points are not equal if any of their components are clearly different.
    The definition of "clearly" different components is controlled by
    the scaling and tolerance in ScaledComparison.
  */
  bool operator!=(const CachePoint& pt) const;

private:

  //! pointer to data if this data is actually owned by the cached point
  Vector* xPtr;

  //! Reference to the data, either internally or externally owned
  const Vector& x;

  //! Vector of function values
  Vector f;

  //! Vector of nonlinear equality constraint values.
  Vector cEqs;

  //! Vector of nonlinear inequality constraint values.
  Vector cIneqs;

};

}
#endif
