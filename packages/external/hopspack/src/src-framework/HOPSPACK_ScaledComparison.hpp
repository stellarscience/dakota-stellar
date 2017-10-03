// $Id: HOPSPACK_ScaledComparison.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_ScaledComparison.hpp $ 

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
  @file HOPSPACK_ScaledComparison.hpp
  @brief Class declaration for HOPSPACK::ScaledComparison.
*/

#ifndef HOPSPACK_SCALEDCOMPARISON_HPP
#define HOPSPACK_SCALEDCOMPARISON_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Performs scaled comparison tests between two n-dimensional vectors.
/*! The class provides static methods so that CacheSplayTree can be coded
    with templates.

    The test between points is based on a user-supplied scaling vector and
    comparison tolerance.
    The comparisons are used to define whether two points are "equal"
    for purposes of the conveyor cache, and if "not equal", then which one
    is lexicographically "greater than" the other.
    See the methods for the numerical definitions.

    If no user parameter is supplied, then default scale factors = 1,
    and the default tolerance = 0.
 */
//----------------------------------------------------------------------
class ScaledComparison
{
  public:

    //! Set the scale factors.
    /*!
     *  @param[in] cScaling  Vector of scale factors, one per component.
     *                       Factors must be positive.
     *  @throws FATAL_ERROR  If scale factors are not positive.
     */
    static void  setScaling (const Vector &  cScaling);

    //! Set the test tolerance.
    /*!
     *  @param[in] dTolerance  Tolerance for declaring two components "equal".
     *  @throws FATAL_ERROR    If tolerance is negative.
     */
    static void  setTolerance (const double  dTolerance);

    //! Return true if PointA = PointB.
    /*!
     *  The comparison is based on scaled components of the vectors.
     *  Let the two input points be denoted as \f$a\f$ and \f$b\f$.
     *  Let \f$s\f$ be the scaling vector and \f$\tau\f$ the tolerance.
     *  Return true if
     *  \f[
     *    |a_i - b_i| \leq \tau s_i \; \mbox{  for all } i = 1, \dots, n
     *  \f]
     *
     *  @param[in] cPointA  First point.
     *  @param[in] cPointB  Second point.
     *  @return true        If cPointA "equals" cPointB.
     */
    static bool  isEqual (const Vector &  cPointA,
                          const Vector &  cPointB);

    //! Return true if PointA != PointB.
    /*!
     *  The comparison is based on scaled components of the vectors.
     *  Let the two input points be denoted as \f$a\f$ and \f$b\f$.
     *  Let \f$s\f$ be the scaling vector and \f$\tau\f$ the tolerance.
     *  Return true if
     *  \f[
     *    |a_i - b_i| > \tau s_i \; \mbox{  for any } i = 1, \dots, n
     *  \f]
     *
     *  @param[in] cPointA  First point.
     *  @param[in] cPointB  Second point.
     *  @return true        If cPointA does not "equal" cPointB.
     */
    static bool  isNotEqual (const Vector &  cPointA,
                             const Vector &  cPointB);

    //! Return true if PointA > PointB.
    /*!
     *  The comparison is based on lexicographic coordinate ordering, using
     *  scaled components of the vectors.
     *  Let the two input points be denoted as \f$a\f$ and \f$b\f$.
     *  Let \f$s\f$ be the scaling vector and \f$\tau\f$ the tolerance.
     *  The logic is <br>
     *  &nbsp;&nbsp;&nbsp; For i=0, i<n <br>
     *  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; if \f$ |a_i - b_i| > \tau s_i \f$
     *    (components are not equal) <br>
     *  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
     *    then Return \f$ a_i - b_i > \tau s_i \f$ <br>
     *  &nbsp;&nbsp;&nbsp; End For <br>
     *  &nbsp;&nbsp;&nbsp; Return false (points must be equal)
     *
     *  @param[in] cPointA  First point.
     *  @param[in] cPointB  Second point.
     *  @return true        If cPointA "greater than" cPointB.
     */
    static bool  isGreaterThan (const Vector &  cPointA,
                                const Vector &  cPointB);

    //! Return true if PointA < PointB.
    /*!
     *  The comparison is based on lexicographic coordinate ordering, using
     *  scaled components of the vectors.
     *  Let the two input points be denoted as \f$a\f$ and \f$b\f$.
     *  Let \f$s\f$ be the scaling vector and \f$\tau\f$ the tolerance.
     *  The logic is <br>
     *  &nbsp;&nbsp;&nbsp; For i=0, i<n <br>
     *  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; if \f$ |a_i - b_i| > \tau s_i \f$
     *    (components are not equal) <br>
     *  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
     *    then Return \f$ b_i - a_i > \tau s_i \f$ <br>
     *  &nbsp;&nbsp;&nbsp; End For <br>
     *  &nbsp;&nbsp;&nbsp; Return false (points must be equal)
     *
     *  @param[in] cPointA  First point.
     *  @param[in] cPointB  Second point.
     *  @return true        If cPointA "less than" cPointB.
     */
    static bool  isLessThan (const Vector &  cPointA,
                             const Vector &  cPointB);

    //! Print debug information.
    static void  printDebugInfo (void);

private:

    //! The class is a singleton, so default constructor is hidden.
    ScaledComparison (void);
    //! By design, there is no copy constructor.
    ScaledComparison (const ScaledComparison &);
    //! By design, there is no assignment operator.
    ScaledComparison & operator= (const ScaledComparison &);

    //! Throw an exception if two vectors are not the same size.
    static void  checkSizes_ (const Vector &  cPointA,
                              const Vector &  cPointB);

    static bool    _bIsScalingDefined;
    static Vector  _cScalingFactors;
    static double  _dToleranceTau;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_SCALEDCOMPARISON_HPP
