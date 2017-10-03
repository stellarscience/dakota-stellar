// $Id: HOPSPACK_Citizen.cpp 183 2010-12-15 18:22:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/HOPSPACK_Citizen.cpp $ 

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
  @file HOPSPACK_Citizen.cpp
  @brief Partially implement abstract class  HOPSPACK::Citizen.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CallbackToMediator.hpp"
#include "HOPSPACK_Citizen.hpp"
#include "HOPSPACK_CitizenGSS.hpp"
#include "HOPSPACK_CitizenGssMS.hpp"
#include "HOPSPACK_CitizenGssNlc.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Construct a citizen instance
//----------------------------------------------------------------------
Citizen *  Citizen::newInstance (const int                         nIdNumber,
                                 const string             &        sName,
                                 const ParameterList      &        cParams,
                                 const ProblemDef         &        cProbDef,
                                 const LinConstr          &        cLinConstr,
                                       CallbackToMediator * const  pCallback,
                                       Citizen            * const  pParent)
{
    //---- NEED A TYPE PARAMETER TO CHOOSE THE CONCRETE IMPLEMENTATION.
    if (cParams.isParameter ("Type") == false)
    {
        cerr << "ERROR found in '" << sName << "' input parameter list:" << endl;
        cerr << "  Required parameter 'Type' is missing." << endl;
        return( NULL );
    }

    string  sCtznType = cParams.getParameter ("Type", "");

    Citizen *  pResult = NULL;
    if (pParent == NULL)
        pResult = makeNewParentInstance_ (sCtznType, nIdNumber, sName,
                                          cParams, cProbDef, cLinConstr,
                                          pCallback);
    else
        pResult = makeNewChildInstance_ (sCtznType, nIdNumber, sName,
                                         cParams, cProbDef, cLinConstr,
                                         pCallback, pParent);
    if (pResult == NULL)
    {
        cerr << "ERROR found in '" << sName << "' input parameter list:" << endl;
        cerr << "  Citizen Type '" << sCtznType << "' not found." << endl;
    }

    return( pResult );
}


//----------------------------------------------------------------------
//  Constructor for Base class
//----------------------------------------------------------------------
Citizen::Citizen (const ParameterList &  cParams,
                  const string        &  sName)
{
    //---- SET THE PRIORITY (SUBCLASSES CANNOT MODIFY THIS).
    _nPriority = cParams.getParameter ("Citizen Priority", 1);
    if (_nPriority < 1)
    {
        cerr << "WARNING: 'Citizen Priority' in '" << sName
             << "' sublist is too small, changing to 1" << endl;
        _nPriority = 1;
    }
    else if (_nPriority > 10)
    {
        cerr << "WARNING: 'Citizen Priority' in '" << sName
             << "' sublist is too large, changing to 10" << endl;
        _nPriority = 10;
    }

    _bShouldIgnoreOtherPoints
        = cParams.getParameter ("Ignore Other Points", false);

    return;
}


//----------------------------------------------------------------------
//  Base class needs a destructor even though it's declared pure virtual.
//----------------------------------------------------------------------
Citizen::~Citizen (void)
{
    return;
}


//----------------------------------------------------------------------
//  Default implementation for Method callbackFromChild
//----------------------------------------------------------------------
void  Citizen::callbackFromChild (const int          nIdNumber,
                                  const int          nReturnCode,
                                  const DataPoint &  cFinalPoint,
                                  const int          nTotalEvals)
{
    //---- DO NOTHING.
    //---- SUBCLASSES THAT ACCEPT A CALLBACK WILL OVERRIDE THIS METHOD.
    return;
}


//----------------------------------------------------------------------
//  Default implementation for Method setEarlyExit
//----------------------------------------------------------------------
void  Citizen::setEarlyExit (void)
{
    //---- IGNORE THE FLAG.
    //---- SUBCLASSES SHOULD OVERRIDE THIS METHOD TO ENSURE ANY SUBSEQUENT
    //---- CALL TO exchange() DOES NOT RETURN NEW TRIAL POINTS.
    return;
}


//----------------------------------------------------------------------
//  Method getPriority
//----------------------------------------------------------------------
int  Citizen::getPriority (void) const
{
    return( _nPriority );
}


//----------------------------------------------------------------------
//  Method shouldIgnoreOtherPoints
//----------------------------------------------------------------------
bool  Citizen::shouldIgnoreOtherPoints (void) const
{
    return( _bShouldIgnoreOtherPoints );
}


//----------------------------------------------------------------------
//  Private Method makeNewParentInstance_
//----------------------------------------------------------------------
Citizen *  Citizen::makeNewParentInstance_
               (const string             &        sCtznType,
                const int                         nIdNumber,
                const string             &        sName,
                const ParameterList      &        cParams,
                const ProblemDef         &        cProbDef,
                const LinConstr          &        cLinConstr,
                      CallbackToMediator * const  pCallback)
{
    Citizen *  pNewCitizen = NULL;

    if (sCtznType == "GSS")
    {
        pNewCitizen = new CitizenGSS (nIdNumber, sName,
                                      cParams, cProbDef, cLinConstr,
                                      NULL);
    }
    else if (sCtznType == "GSS-MS")
    {
        pNewCitizen = new CitizenGssMS (nIdNumber, sName,
                                        cParams, cProbDef, cLinConstr,
                                        pCallback);
    }
    else if (sCtznType == "GSS-NLC")
    {
        pNewCitizen = new CitizenGssNlc (nIdNumber, sName,
                                         cParams, cProbDef, cLinConstr,
                                         pCallback, NULL);
    }

    return( pNewCitizen );
}


//----------------------------------------------------------------------
//  Private Method makeNewChildInstance_
//----------------------------------------------------------------------
Citizen *  Citizen::makeNewChildInstance_
               (const string             &        sCtznType,
                const int                         nIdNumber,
                const string             &        sName,
                const ParameterList      &        cParams,
                const ProblemDef         &        cProbDef,
                const LinConstr          &        cLinConstr,
                      CallbackToMediator * const  pCallback,
                      Citizen            * const  pParent)
{
    Citizen *  pNewCitizen = NULL;

    if (sCtznType == "GSS-child")
    {
        pNewCitizen = new CitizenGSS (nIdNumber, sName,
                                      cParams, cProbDef, cLinConstr,
                                      pParent);
    }
    else if (sCtznType == "GSS-NLC-child")
    {
        pNewCitizen = new CitizenGssNlc (nIdNumber, sName,
                                         cParams, cProbDef, cLinConstr,
                                         pCallback, pParent);
    }

    return( pNewCitizen );
}

    
}     //-- namespace HOPSPACK
