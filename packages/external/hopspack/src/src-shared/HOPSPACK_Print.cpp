// $Id: HOPSPACK_Print.cpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_Print.cpp $ 

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
  \file HOPSPACK_Print.cpp
  \brief Implement HOPSPACK::Print.
*/

#include "HOPSPACK_float.hpp"
#include "HOPSPACK_Print.hpp"
#include <iomanip>

namespace HOPSPACK
{

//---- STATIC MEMBER DEFAULTS.
Print::PrintType  Print::QUEUE_LISTS = MOST_VERBOSE;
int Print::precision = 3;
Print::PrintType Print::_nDisplayLevel = Print::INPUT_PARAMETERS;

void  Print::setDisplayParameter (const PrintType  nDisplayLevel)
{
    //---- KEEP THE MAXIMUM SYNCHRONIZED WITH PrintType IN THE HEADER FILE.
    if ((nDisplayLevel >= 0) && (nDisplayLevel <= 5))
        _nDisplayLevel = nDisplayLevel;
    return;
}

void  Print::setPrecisionParameter (const int  nPrecision)
{
    precision = nPrecision;
    return;
}

bool  Print::doPrint(const Print::PrintType nPrintLevel)
{
  return (_nDisplayLevel >= nPrintLevel);
}

int  Print::getPrecision (void)
{
    return (precision);
}


}     //-- namespace HOPSPACK

