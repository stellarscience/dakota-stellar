// $Id: HOPSPACK_gssChildReturnCodes.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-common/HOPSPACK_gssChildReturnCodes.hpp $

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
  @file HOPSPACK_gssChildReturnCodes.hpp
  @brief Define return codes from any GSS child citizen to its GSS parent.
*/

#ifndef HOPSPACK_GSSCHILDRETURNCODES_HPP
#define HOPSPACK_GSSCHILDRETURNCODES_HPP


namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Define enumerations
//----------------------------------------------------------------------

enum GssChildReturnCodesType
{
    //! Subproblem experienced a serious error and did not complete.
    REASON_ERROR = 0,

    //! Subproblem converged and is returning the solution point.
    REASON_CONVERGED,

    //! Subproblem used its quota of evaluations; returning best point.
    REASON_NO_MORE_EVALS,

    //! Subproblem was halted by the Mediator; returning best point.
    REASON_HALTED_BY_MEDIATOR
};


//----------------------------------------------------------------------
//  Function declarations in the HOPSPACK namespace
//----------------------------------------------------------------------

void  gssChildPrintReturnCode (const GssChildReturnCodesType  nCode);
    

}


#endif     //-- HOPSPACK_GSSCHILDRETURNCODES_HPP
