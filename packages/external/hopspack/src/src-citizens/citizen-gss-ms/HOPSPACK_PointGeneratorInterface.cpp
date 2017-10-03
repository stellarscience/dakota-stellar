// $Id: HOPSPACK_PointGeneratorInterface.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_PointGeneratorInterface.cpp $

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
  @file HOPSPACK_PointGeneratorInterface.cpp
  @brief Partially implement abstract class  HOPSPACK::PointGeneratorInterface.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_GeneratorTBD.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_PointGeneratorInterface.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Construct a point generator instance
//----------------------------------------------------------------------
PointGenerator *  PointGenerator::newInstance
                      (const string     &  sName,
                       const int           nNumToGenerate,
                       const ProblemDef &  cProbDef,
                       const LinConstr  &  cLinConstr)
{
    //---- FIND THE GENERATOR NAME AND CONSTRUCT.
    PointGenerator *  pNewGenerator = NULL;
    if (sName == "TBD")
    {
        try
        {
            pNewGenerator = new GeneratorTBD (nNumToGenerate,
                                              cProbDef,
                                              cLinConstr);
        }
        catch (const char * const)
        {
            pNewGenerator = NULL;
        }
    }
    else
    {
        cerr << "ERROR: Unknown point generator '" << sName << "' for"
             << " GSS-MS" << endl;
        return( NULL );
    }

    return( pNewGenerator );
}


//----------------------------------------------------------------------
//  Constructor for Base class
//----------------------------------------------------------------------
PointGenerator::PointGenerator (void)
{
    return;
}


//----------------------------------------------------------------------
//  Base class needs a destructor even though it's declared pure virtual.
//----------------------------------------------------------------------
PointGenerator::~PointGenerator (void)
{
    return;
}


}     //-- namespace HOPSPACK
