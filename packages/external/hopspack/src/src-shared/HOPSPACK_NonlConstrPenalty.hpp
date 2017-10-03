// $Id: HOPSPACK_NonlConstrPenalty.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_NonlConstrPenalty.hpp $

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
  @file HOPSPACK_NonlConstrPenalty.hpp
  @brief Class declaration for HOPSPACK::NonlConstrPenalty.
*/
#ifndef HOPSPACK_NONLCONSTRPENALTY_HPP
#define HOPSPACK_NONLCONSTRPENALTY_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Defines and computes a penalty term for nonlinear constraints.
/*!
 *  Supported penalty functions are:
 *    None           - returns zero for any constraint values.
 *    L2 squared     - |(c_eq, c_ineq)|_2^2
 *    L1 nonsmooth   - |(c_eq, c_ineq)|_1
 *    L1 smoothed    - Chen and Mangasarian function of (c_eq, c_ineq)
 *    L2 nonsmooth   - |(c_eq, c_ineq)|_2
 *    L2 smoothed    - Positive shift of square root in |(c_eq, c_ineq)|_2
 *    Linf nonsmooth - |(c_eq, c_ineq)|_inf
 *    Linf smoothed  - Liuzi and Lucidi function of (c_eq, c_ineq)
 *  In each case the function is multiplied by the penalty coefficient.
 */
//----------------------------------------------------------------------
class NonlConstrPenalty
{
  public:

    //! Default constructor, penalty term is always zero.
    NonlConstrPenalty (void);

    //! Destructor.
    ~NonlConstrPenalty (void);


    //! Define a penalty function.
    /*!
     *  @param[in] sPenaltyName         Name of the penalty function.
     *  @param[in] dPenaltyCoefficient  Multiply the function by this.
     *  @param[in] dSmoothingFactor     Factor if function is smoothed.
     *  @return                         False if any parameter is incorrect.
     */
    bool  defineFunction (const string &  sPenaltyName,
                          const double    dPenaltyCoefficient,
                          const double    dSmoothingFactor);

    //! Return true if a penalty function is defined.
    bool  isDefined (void) const;

    //! Change the penalty coefficient to a new value.
    void  updateCoefficient (const double  dPenaltyCoefficient);

    //! Change the smoothing factor to a new value.
    void  updateSmoothing (const double  dSmoothingFactor);

    //! Return the name of the penalty function.
    const string &  getPenaltyName (void) const;

    //! Return the penalty coefficient.
    double  getCoefficient (void) const;

    //! Return the smoothing factor, or zero if the function is not smoothed.
    double  getSmoothing (void) const;

    //! Return the computed penalty function, given constraint values.
    /*!
     *  See the comments about the class for a list of supported functions.
     *  If defineFunction() has not been called, then return zero.
     */
    double  computePenalty (const Vector &  cEqs,
                            const Vector &  cIneqs) const;

    //! Print information about the instance.
    void  printDefinition (void) const;
    

  private:

    //! By design, there is no copy constructor.
    NonlConstrPenalty (const NonlConstrPenalty &);
    //! By design, there is no assignment operator.
    NonlConstrPenalty & operator= (const NonlConstrPenalty &);

    //! Return the L2 Squared penalty function.
    double  computeL2Sqrd_ (const Vector &  cEqs,
                            const Vector &  cIneqs) const;

    //! Return the L1 penalty function.
    double  computeL1_ (const Vector &  cEqs,
                        const Vector &  cIneqs) const;

    //! Return the L1 smoothed penalty function.
    double  computeL1Smoothed_ (const Vector &  cEqs,
                                const Vector &  cIneqs) const;

    //! Return the L2 penalty function.
    double  computeL2_ (const Vector &  cEqs,
                        const Vector &  cIneqs) const;

    //! Return the L2 smoothed penalty function.
    double  computeL2Smoothed_ (const Vector &  cEqs,
                                const Vector &  cIneqs) const;

    //! Return the Linf penalty function.
    double  computeLinf_ (const Vector &  cEqs,
                          const Vector &  cIneqs) const;

    //! Return the Linf smoothed penalty function.
    double  computeLinfSmoothed_ (const Vector &  cEqs,
                                  const Vector &  cIneqs) const;

    //! Return |c_eq|_2^2 + |c_ineq|_2^2.
    double  computeSumSqs_ (const Vector &  cEqs,
                            const Vector &  cIneqs) const;


    //! Type of penalty function.
    int     _nPenType;

    //! Coefficient for the entire penalty term.
    double  _dPenCoef;

    //! Smoothing paramter if the penalty function is to be smoothed.
    double  _dSmoothingFactor;
};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_NONLCONSTRPENALTY_HPP
