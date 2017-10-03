// $Id: HOPSPACK_MultiStartRepository.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_MultiStartRepository.cpp $

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
  @file HOPSPACK_MultiStartRepository.cpp
  @brief Implement HOPSPACK::MultiStartRepository.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_MultiStartRepository.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
MultiStartRepository::MultiStartRepository (const ProblemDef &  cProbDef,
                                            const LinConstr  &  cLinConstr,
                                            const double        dComparisonTol)
    : _cProbDef (cProbDef),
      _cLinConstr (cLinConstr),
      _dComparisonTol (dComparisonTol)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
MultiStartRepository::~MultiStartRepository (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method addResult
//----------------------------------------------------------------------
void  MultiStartRepository::addResult (const DataPoint &  cNewResult)
{
 cout << "TBD supposed to addResult ";
 cNewResult.leftshift(cout,true);
 cout << endl;
    //TBD
    return;
}


//----------------------------------------------------------------------
//  Method getBestResult
//----------------------------------------------------------------------
bool  MultiStartRepository::getBestResult (DataPoint &  cBestResult) const
{
 cout << "TBD supposed to getBestResult, not working yet\n";
    //TBD
    return( false );
}


//----------------------------------------------------------------------
//  Method getBestResult
//----------------------------------------------------------------------
void  MultiStartRepository::getBestResultList
          (const int                    nMaxNumResults,
                 vector< DataPoint * >  cResults) const
{
    cResults.clear();

    //TBD
 cout << "TBD supposed to getBestResultList, not working yet\n";
    return;
}


}     //-- namespace HOPSPACK
