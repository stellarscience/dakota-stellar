// $Id: HOPSPACK_CitizenGssMS.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-ms/HOPSPACK_CitizenGssMS.cpp $

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
  @file HOPSPACK_CitizenGssMS.cpp
  @brief Implement HOPSPACK::CitizenGssMS, subclass of Citizen.
*/

#include <sstream>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CitizenGssMS.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_PointGeneratorInterface.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Declare internal enumerations and constants
//----------------------------------------------------------------------

enum {
    FIN_REASON_CHILD_HALTED = 0,
    FIN_REASON_CHILD_ERROR,
    FIN_REASON_TRIED_ALL_POINTS
};

static const int  nCHILD_IS_IDLE = -999999;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
CitizenGssMS::CitizenGssMS (const int                         nIdNumber,
                            const string             &        sName,
                            const ParameterList      &        cParams,
                            const ProblemDef         &        cProbDef,
                            const LinConstr          &        cLinConstr,
                                  CallbackToMediator * const  pCallback)
    : Citizen (cParams, sName),
      _nID (nIdNumber),
      _sName (sName + " (GSS-MS)"),
      _nState (CONTINUE),
      _cProbDef (cProbDef),
      _cLinConstr (cLinConstr),
      _cParentParams (cParams),
      _pCallbackToMed (pCallback),
      _pGenerator (NULL)
{
    //---- INITIALIZE PARAMETERS FROM THE INPUT LIST.
    if (extractParameters_ (_cParentParams, _cSubprobParams) == false)
    {
        throw "GSS-MS Error";
    }

    //---- CONSTRUCT AN INFORMATION BLOCK FOR EACH POTENTIAL CHILD.
    _cChildren.resize (_nNumConcurrentSubprobs);
    for (int  i = 0; i < (int) _cChildren.size(); i++)
    {
        ChildCtznInfoBlockType *  pChildInfo = new ChildCtznInfoBlockType();
        pChildInfo->nChildID = nCHILD_IS_IDLE;
        pChildInfo->pProbDef = NULL;
        _cChildren[i] = pChildInfo;
    }

    _nCurrentNumStartPoints = 0;
    _nTotalEvals = 0;

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
CitizenGssMS::~CitizenGssMS (void)
{
    if (_pGenerator != NULL)
        delete _pGenerator;

    for (int  i = 0; i < (int) _cChildren.size(); i++)
    {
        ChildCtznInfoBlockType *  pChildInfo = _cChildren[i];
        if (pChildInfo->pProbDef != NULL)
            delete pChildInfo->pProbDef;
        delete pChildInfo;
    }
    _cChildren.erase (_cChildren.begin(), _cChildren.end());

    return;
}


//----------------------------------------------------------------------
//  Method preProcess
//----------------------------------------------------------------------
void  CitizenGssMS::preProcess (void)
{
    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        cout << "###################################################" << endl;
        cout << "###   HOPSPACK GSS-MS Initialization Results    ###" << endl;
        cout << "###   Citizen name: " << getName() << endl;
        cout << endl;
        cout << "  WARNING - HOPSPACK 2.0 default multi-start logic" << endl;
        cout << "            is extremely simple-minded!" << endl;
        cout << "            Please consider writing your own version." << endl;
        cout << endl;

        cout << "Priority = " << getPriority()
             << "  (1=highest, 10=lowest)" << endl;
        cout << endl;

        cout << "*** Parameter List (alphabetical order) ***" << endl;
        _cParentParams.print();
        cout << endl;

        _cProbDef.printDefinition (false);
        _cLinConstr.printDefinition (false);

        _pGenerator->printDefinition();
        cout << endl;

        cout << "### End HOPSPACK GSS-MS Initialization Results  ###" << endl;
        cout << "###################################################" << endl;
    }

    //---- ADD CHILD PARAMETERS APPROPRIATE FOR THE TYPE OF SUBPROBLEM.
    if (_cProbDef.hasNonlinearConstr())
    {
        _cSubprobParams.setParameter ("Type", "GSS-NLC-child");
        _cSubprobParams.setParameter ("Display Subproblem", _nDisplaySubLevel);
        //TBD...have to make this a fraction of nmaxsubprobevals
        _cSubprobParams.setParameter ("Max Subproblem Evaluations",
                                      _nMaxSubprobEvals);
    }
    else
    {
        _cSubprobParams.setParameter ("Type", "GSS-child");
        _cSubprobParams.setParameter ("Maximum Evaluations", _nMaxSubprobEvals);
    }
    _cSubprobParams.setParameter ("Display", _nDisplaySubLevel);
    _cSubprobParams.setParameter ("Ignore Other Points", true);

    //---- ANY ERROR CREATING A CHILD IS FATAL TO THE CITIZEN.
    if (nextIteration_() == false)
    {
        _nState = FINISHED;
        _nFinishedReason = FIN_REASON_CHILD_ERROR;
        return;
    }

    return;
}


//----------------------------------------------------------------------
//  Method postProcess
//----------------------------------------------------------------------
void  CitizenGssMS::postProcess (void)
{
    if (_nState == WAITING)
        _nState = FINISHED;

    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        if (_nState == FINISHED)
        {
            cout << " GSS-MS complete: ";
            if (_nFinishedReason == FIN_REASON_TRIED_ALL_POINTS)
                cout << "Finished all start points";
            else if (_nFinishedReason == FIN_REASON_CHILD_HALTED)
                cout << "Could not proceed after subproblem halted";
            else if (_nFinishedReason == FIN_REASON_CHILD_ERROR)
                cout << "Could not proceed after subproblem error";
            else
                cout << "Unknown reason!";
            cout << endl;
        }
        else
        {
            cout << " GSS-MS did not complete" << endl;
        }

        cout << "  Number of subproblems solved                        = "
             << _nCurrentNumStartPoints++ << endl;
        //TBD...print how many points are in the repository
        cout << "  Evaluated points from this citizen and its children = "
             << _nTotalEvals << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Method exchange
//----------------------------------------------------------------------
void  CitizenGssMS::exchange (const list< DataPoint * > &  cResultList,
                              const list< int >         &  cOwnerTags,
                                    list< DataPoint * > &  cWaitList)
{
    //---- TBD handle evals requested by point generator
    return;
}


//----------------------------------------------------------------------
//  Method getIdNumber
//----------------------------------------------------------------------
int  CitizenGssMS::getIdNumber (void) const
{
    return( _nID );
}


//----------------------------------------------------------------------
//  Method getName
//----------------------------------------------------------------------
const string &  CitizenGssMS::getName (void) const
{
    return( _sName );
}


//----------------------------------------------------------------------
//  Method getState
//----------------------------------------------------------------------
Citizen::State  CitizenGssMS::getState (void) const
{
    return( _nState );
}


//----------------------------------------------------------------------
//  Method callbackFromChild
//----------------------------------------------------------------------
void  CitizenGssMS::callbackFromChild (const int          nIdNumber,
                                       const int          nReturnCode,
                                       const DataPoint &  cFinalPoint,
                                       const int          nTotalEvals)
{
    ChildCtznInfoBlockType *  pChildInfo = NULL;
    for (int  i = 0; i < (int) _cChildren.size(); i++)
    {
        if (_cChildren[i]->nChildID == nIdNumber)
        {
            pChildInfo = _cChildren[i];
            break;
        }
    }
    if (pChildInfo == NULL)
    {
        cerr << "ERROR: Parent citizen '" << getName() << "' received callback"
             << " from unknown child = " << nIdNumber << endl;
        _nState = FINISHED;
        _nFinishedReason = FIN_REASON_CHILD_ERROR;
        return;
    }

    _nCurrentNumStartPoints++;
    _nTotalEvals += nTotalEvals;

    GssChildReturnCodesType  nCode = (GssChildReturnCodesType) nReturnCode;

    if (_nDisplayLevel >= 2)
    {
        cout << endl;
        cout << " " << getName() << " received callback from child "
             << nIdNumber << endl;
        gssChildPrintReturnCode (nCode);
    }
    if (nCode == REASON_ERROR)
    {
        cerr << "WARNING: Child citizen failed to solve subproblem" << endl;
    }
    if (_nDisplayLevel >= 2)
    {
        //TBD...detect invalid point?
        cout << " GSS-MS subproblem solution:" << endl;
        cFinalPoint.leftshift (cout, false);
        cout << endl;
    }

    //TBD...save the point in the respository
    //TBD...mark it if subproblem ran out of evaluations

    //---- RESTORE THIS CHILD TO AN IDLE STATE.
    if (pChildInfo->pProbDef != NULL)
        delete pChildInfo->pProbDef;
    pChildInfo->pProbDef = NULL;
    pChildInfo->nChildID = nCHILD_IS_IDLE;

    if (isTimeToStop_ (nCode, cFinalPoint))
    {
        if (_nDisplayLevel >= 1)
            cout << " GSS-MS '" << getName()
                 << "' is finished, waiting for subproblems to complete"
                 << endl << endl;
        _nState = WAITING;
        return;
    }

    //---- IF THERE ARE NO MORE START POINTS, THEN STATE WILL BE SET TO FINISHED.
    if (nextIteration_() == false)
    {
        _nState = FINISHED;
        _nFinishedReason = FIN_REASON_CHILD_ERROR;
        return;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method extractParameters_
//----------------------------------------------------------------------
bool  CitizenGssMS::extractParameters_ (ParameterList &  cParams,
                                        ParameterList &  cRemainder)
{
    cRemainder = cParams;
    cRemainder.deleteParameter ("Type");

    //---- DISPLAY LEVEL.
    _nDisplayLevel = cParams.getOrSetParameter ("Display", 0);
    if (_nDisplayLevel < 0)
        _nDisplayLevel = 0;
    if (_nDisplayLevel > 2)
        _nDisplayLevel = 2;
    cRemainder.deleteParameter ("Display");

    //---- SUBPROBLEM DISPLAY LEVEL.
    _nDisplaySubLevel
        = cParams.getOrSetParameter ("Display Subproblem", 0);
    if (_nDisplaySubLevel < 0)
        _nDisplaySubLevel = 0;
    if (_nDisplaySubLevel > 3)
        _nDisplaySubLevel = 3;
    cRemainder.deleteParameter ("Display Subproblem");

    _nMaxSubprobEvals
        = cParams.getOrSetParameter ("Max Subproblem Evaluations", -1);
    if (_nMaxSubprobEvals < -1)
        _nMaxSubprobEvals = -1;
    cRemainder.deleteParameter ("Max Subproblem Evaluations");

    _nTotalStartPoints = min (100, _cProbDef.getVarScaling().size() * 5);
    _nTotalStartPoints
        = cParams.getOrSetParameter ("Total Start Points", _nTotalStartPoints);
    if (_nTotalStartPoints <= 0)
    {
        cerr << "ERROR: Invalid nonpositive value for 'Total Start Points'"
             << " in sublist 'GSS-MS'" << endl;
        return( false );
    }
    cRemainder.deleteParameter ("Total Start Points");

    _nNumConcurrentSubprobs
        = cParams.getOrSetParameter ("Concurrent Subproblems", 1);
    if (_nNumConcurrentSubprobs < 0)
        _nNumConcurrentSubprobs = 1;
    if (_nNumConcurrentSubprobs > _nTotalStartPoints)
        _nNumConcurrentSubprobs = _nTotalStartPoints;
    cRemainder.deleteParameter ("Concurrent Subproblems");

    //---- CONSTRUCT A POINT GENERATOR.
    if (cParams.isParameter ("Point Generator") == false)
    {
        cerr << "ERROR: Must specify 'Point Generator' in sublist 'GSS-MS'"
             << endl;
        return( false );
    }
    string  sGen = cParams.getParameter ("Point Generator", "");
    _pGenerator = PointGenerator::newInstance (sGen,
                                               _nTotalStartPoints,
                                               _cProbDef,
                                               _cLinConstr);
    if (_pGenerator == NULL)
    {
        cerr << "ERROR: GSS-MS could not construct point generator" << endl;
        return( false );
    }
    cRemainder.deleteParameter ("Point Generator");

    return( true );
}


//----------------------------------------------------------------------
//  Private Method nextIteration_
//----------------------------------------------------------------------
bool  CitizenGssMS::nextIteration_ (void)
{
    //---- DEFINE CONTAINERS FOR THE POINT GENERATOR TO FILL.
    Vector  cNextLocation;
    vector< const DataPoint * >  cNextEvalsList;

    while (true)
    {
        //---- FIND A CHILD CITIZEN SLOT THAT IS IDLE.
        ChildCtznInfoBlockType *  pChildInfo = NULL;
        for (int  i = 0; i < (int) _cChildren.size(); i++)
        {
            if (_cChildren[i]->nChildID == nCHILD_IS_IDLE)
            {
                pChildInfo = _cChildren[i];
                break;
            }
        }
        if (pChildInfo == NULL)
        {
            //---- MAXIMUM NUMBER OF CONCURRENT CHILDREN ALREADY RUNNING.
            return( true );
        }

        //---- GENERATE THE NEXT START POINT.
        if (_pGenerator->getNextPoint (cNextLocation, cNextEvalsList) == false)
        {
            _nState = FINISHED;
            _nFinishedReason = FIN_REASON_TRIED_ALL_POINTS;
            return( true );
        }

        if (cNextEvalsList.size() > 0)
        {
            //---- THE GENERATOR SUBMITTED A LIST OF POINTS TO BE EVALUATED.
            //TBD...submit points for evaluation, delete ptrs
 cout << "TBD GSS-MS cannot eval genpts yet\n";
 _nState = FINISHED;
 _nFinishedReason = FIN_REASON_CHILD_ERROR; //TBD...phony reason
            for (int  i = 0; i < (int) cNextEvalsList.size(); i++)
                delete cNextEvalsList[i];
            return( true );
        }


        //---- START A NEW CHILD CITIZEN WITH THE NEW START POINT.

        if (isStartPointOK_ (cNextLocation) == false)
            return( false );

        pChildInfo->pProbDef = new ProblemDef (_cProbDef);
        pChildInfo->pProbDef->resetInitialX (cNextLocation);

        if (_nDisplayLevel >= 2)
        {
            cout << " GSS-MS starting new subproblem at the point:" << endl;
            cout << "  x=[";
            cNextLocation.leftshift (cout, true);
            cout << "]" << endl;
        }

        int  nReturnID = -1;
        Citizen *  pChild = NULL;
        try
        {
            nReturnID = _pCallbackToMed->reserveUniqueCitizenID();
            stringstream  tmp;
            tmp << "Citizen " << nReturnID << " (child of " << _nID << ")";
            string  sChildCtznName = tmp.str();

            pChild = Citizen::newInstance (nReturnID,
                                           sChildCtznName,
                                           _cSubprobParams,
                                           *(pChildInfo->pProbDef),
                                           _cLinConstr,
                                           _pCallbackToMed,
                                           this);
        }
        catch (const char * const)
        {
            nReturnID = -1;
        }
        if ((nReturnID == -1) || (pChild == NULL))
        {
            cerr << "ERROR: Failed to create child citizen" << endl;
            return( false );
        }
        pChildInfo->nChildID = nReturnID;

        //---- ADD THE CHILD TO THE MEDIATOR, WHICH WILL START RUNNING IT.
        if (_pCallbackToMed->addChildCitizen (pChild, getIdNumber()) == false)
        {
            cerr << "ERROR: Failed to add child citizen for GSS-MS" << endl;
            return( false );
        }

        if (_nDisplayLevel >= 2)
            cout << " CitizenGssMS started child citizen " << nReturnID
                 << endl << endl;

    }     //-- END while (true)
}


//----------------------------------------------------------------------
//  Private Method isStartPointOK_
//----------------------------------------------------------------------
bool  CitizenGssMS::isStartPointOK_ (const Vector &  cStartPoint) const
{
    if (cStartPoint.size() != _cProbDef.getVarScaling().size())
    {
        cerr << "ERROR: Length of generated start point = "
             << cStartPoint.size()
             << " does not match number of unknowns" << endl;
        return( false );
    }
    if (_cProbDef.isBndsFeasible (cStartPoint) == false)
    {
        cerr << "ERROR: Generated start point violates variable bounds"
             << endl;
        return( false );
    }
    if (_cLinConstr.isFeasible (cStartPoint, true) == false)
    {
        cerr << "ERROR: Generated start point violates linear constraints"
             << endl;
        return( false );
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method isTimeToStop_
//----------------------------------------------------------------------
bool  CitizenGssMS::isTimeToStop_ (const GssChildReturnCodesType  nReturnCode,
                                   const DataPoint &              cTestPoint)
{
    if (nReturnCode == REASON_HALTED_BY_MEDIATOR)
    {
        _nFinishedReason = FIN_REASON_CHILD_HALTED;
        return( true );
    }

    //---- POINT GENERATOR SHOULD NOT EXCEED THE TOTAL REQUESTED START POINTS,
    //---- BUT CHECK JUST IN CASE.
    if (_nCurrentNumStartPoints >= _nTotalStartPoints)
    {
        _nFinishedReason = FIN_REASON_TRIED_ALL_POINTS;
        return( true );
    }

    return( false );
}


}     //-- namespace HOPSPACK
