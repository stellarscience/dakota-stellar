// $Id: HOPSPACK_GeneratorTBD.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_GeneratorTBD.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
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
  @file HOPSPACK_GeneratorTBD.cpp
  @brief Implement HOPSPACK::GeneratorTBD.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_GeneratorTBD.hpp"
#include "HOPSPACK_PointGeneratorInterface.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
GeneratorTBD::GeneratorTBD (const int           nNumToGenerate,
                            const ProblemDef &  cProbDef,
                            const LinConstr  &  cLinConstr)
    : PointGenerator(),
      _nMaxNumPoints (nNumToGenerate),
      _cProbDef (cProbDef),
      _cLinConstr (cLinConstr)
{
    _nPointNum = 0;
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
GeneratorTBD::~GeneratorTBD (void)
{
    //TBD
    return;
}


//----------------------------------------------------------------------
//  Method getNextPoint
//----------------------------------------------------------------------
bool  GeneratorTBD::getNextPoint (Vector &                       cStartLocation,
                                  vector< const DataPoint * > &  cEvalList)
{
    //TBD...the test implementation is extremely simple (and stupid)
// cout << "++++++++++++++++++++ TBD enter generatorTBD::getNext " << _nPointNum << endl;

    if (_nPointNum >= _nMaxNumPoints)
        return( false );
    _nPointNum++;

    cStartLocation.resize (_cProbDef.getVarScaling().size());

    //---- IF THE PROBLEM DEFINES AN INITIAL START POINT, USE IT AS THE FIRST
    //---- GENERATED POINT.
    if (_nPointNum == 1)
    {
        Vector  cInitialX = _cProbDef.getInitialX();
        if (cInitialX.empty() == false)
        {
            cStartLocation = cInitialX;
            return( true );
        }
    }

    const Vector &  cLower = _cProbDef.getLowerBnds();
    const Vector &  cUpper = _cProbDef.getUpperBnds();
    for (int  i = 0; i < cStartLocation.size(); i++)
    {
        double  dLowerBnd = cLower[i];
        if (exists (dLowerBnd) == false)
            dLowerBnd = -1.0;
        double  dUpperBnd = cUpper[i];
        if (exists (dUpperBnd) == false)
            dUpperBnd =  1.0;

        double  dRnd = genRandomNumber();
        cStartLocation[i] = dLowerBnd + dRnd * (dUpperBnd - dLowerBnd);
    }
    if (_cLinConstr.isFeasible (cStartLocation) == false)
    {
        _cLinConstr.projectToFeasibility (cStartLocation);
    }
    return( true );

 /*TBD...example of how to do an eval list
 Vector cTbd (2);
 cTbd[0] = 1.0;
 cTbd[1] = 2.0;
 cEvalList.push_back (new DataPoint(ProblemDef::MINIMIZE, cTbd));
 */
}


//----------------------------------------------------------------------
//  Method addResultPair
//----------------------------------------------------------------------
void  GeneratorTBD::addResultPair (const DataPoint &  cStartPoint,
                                   const DataPoint &  cResultPoint,
                                   const ResultState  nResultState)
{
    //TBD...the test implementation does nothing
    return;
}


//----------------------------------------------------------------------
//  Method printDefinition
//----------------------------------------------------------------------
void  GeneratorTBD::printDefinition (void) const
{
    cout << "GeneratorTBD definition\n";
    return;
}


}     //-- namespace HOPSPACK
