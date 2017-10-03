// $Id: HOPSPACK_gssChildReturnCodes.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-common/HOPSPACK_gssChildReturnCodes.cpp $

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
  @file HOPSPACK_gssChildReturnCodes.cpp
  @brief Implement functions declared in HOPSPACK_gssChildReturnCodes.cpp.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"


//----------------------------------------------------------------------
//  Function gssChildPrintReturnCode
//----------------------------------------------------------------------
void  HOPSPACK::gssChildPrintReturnCode (const GssChildReturnCodesType  nCode)
{
    cout << "  Return code = " << (int) nCode;

    if      (nCode == REASON_ERROR)
        cout << " (error)";
    else if (nCode == REASON_CONVERGED)
        cout << " (successful convergence)";
    else if (nCode == REASON_NO_MORE_EVALS)
        cout << " (out of evaluations)";
    else if (nCode == REASON_HALTED_BY_MEDIATOR)
        cout << " (halted by Mediator)";
    else
        cout << " (unknown code!)";

    cout << endl;
    return;
}
