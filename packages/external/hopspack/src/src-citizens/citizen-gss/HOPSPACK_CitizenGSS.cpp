// $Id: HOPSPACK_CitizenGSS.cpp 183 2010-12-15 18:22:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_CitizenGSS.cpp $ 

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
  @file HOPSPACK_CitizenGSS.cpp
  @brief Implement HOPSPACK::CitizenGSS, subclass of Citizen.
*/

#include <iomanip>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CitizenGSS.hpp"
#include "HOPSPACK_GssIterator.hpp"
#include "HOPSPACK_GssList.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
CitizenGSS::CitizenGSS (const int                    nIdNumber,
                        const string        &        sName,
                        const ParameterList &        cParams,
                        const ProblemDef    &        cProbDef,
                        const LinConstr     &        cLinConstr,
                              Citizen       * const  pParent)
    : Citizen (cParams, sName),
      _nID (nIdNumber),
      _sName (sName + " (GSS)"),
      _cProbDef (cProbDef),
      _cLinConstr (cLinConstr),
      _pParent (pParent),
      _bWillExitEarly (false)
{
    //---- COPY THE PARAMETERS AND ALLOW NEW ONES TO BE APPENDED.
    _cGssParams = cParams;

    //---- CHECK THAT THE PROBLEM TYPE IS APPROPRIATE FOR THIS CITIZEN.
    if (_cProbDef.isDomainContinuous() == false)
    {
        cerr << "ERROR: GSS citizen can only solve problems with"
             << " a continuous domain" << endl;
        throw "GSS Error";
    }
    if (   (_cProbDef.hasNonlinearConstr() == true)
        && (_cGssParams.isParameter ("Penalty Function") == false) )
    {
        cerr << "ERROR: GSS citizen cannot solve problems with"
             << " nonlinear constraints" << endl;
        throw "GSS Error";
    }

    _nMaxToKeepAfterNewBest
        = _cGssParams.getOrSetParameter ("Maximum Queue Size", 0);
    if (_nMaxToKeepAfterNewBest < 0)
    {
        cerr << "WARNING: Invalid negative 'Maximum Queue Size' in GSS sublist,"
             << " changing to zero" << endl;
        _nMaxToKeepAfterNewBest = 0;
    }

    //---- DISPLAY LEVEL.
    _nDisplayLevel = _cGssParams.getOrSetParameter ("Display", 0);
    if (_nDisplayLevel < 0)
        _nDisplayLevel = 0;
    if (_nDisplayLevel > 3)
        _nDisplayLevel = 3;

    //---- DEFINE A PENALTY FUNCTION FOR ANY NONLINEAR CONSTRAINTS.
    _pPenalty = new NonlConstrPenalty();
    if (_cGssParams.isParameter ("Penalty Function"))
    {
        if (_cGssParams.isParameter ("Penalty Parameter") == false)
        {
            cerr << "ERROR: GSS citizen needs 'Penalty Parameter'"
                 << " value for nonlinear constraints" << endl;
            throw "GSS Error";
        }
        const string &  sPenFn
            = _cGssParams.getParameter ("Penalty Function", "");
        double  dPenParm
            = _cGssParams.getDoubleParameter ("Penalty Parameter");
        double  dSmooth
            = _cGssParams.getParameter ("Penalty Smoothing Value", 0.0);
        if (_pPenalty->defineFunction (sPenFn, dPenParm, dSmooth) == false)
        {
            cerr << "ERROR: GSS citizen could not construct penalty"
                 << " function for nonlinear constraints" << endl;
            throw "GSS Error";
        }
    }

    _pGssIterator = new GssIterator (_cProbDef,
                                     _cLinConstr,
                                     _pPenalty,
                                     _cGssParams);
    _cExchangeList.setDefaultStepLength (_pGssIterator->getInitialStep());
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
CitizenGSS::~CitizenGSS (void)
{
    //---- THIS SHOULD ALREADY BE EMPTY, BUT JUST IN CASE PRUNE IT AGAIN.
    _cExchangeList.prune();

    delete _pPenalty;
    delete _pGssIterator;
    return;
}


//----------------------------------------------------------------------
//  Method preProcess
//----------------------------------------------------------------------
void  CitizenGSS::preProcess (void)
{
    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        cout << "##################################################" << endl;
        cout << "###     HOPSPACK GSS Initialization Results    ###" << endl;
        cout << "###     Citizen name: " << getName() << endl;
        cout << endl;

        cout << "Priority = " << getPriority()
             << "  (1=highest, 10=lowest)" << endl;
        cout << endl;

        _pGssIterator->printInitializationInformation();
        cout << endl;
        _cProbDef.printDefinition (false);
        _cLinConstr.printDefinition (false);
        if (_cProbDef.hasNonlinearConstr())
            _pPenalty->printDefinition();

        cout << "### End HOPSPACK GSS Initialization Results    ###" << endl;
        cout << "##################################################" << endl;
    }

    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        cout << " GSS Start Point:" << endl;
        const GssPoint &  cBest = _pGssIterator->getBestPoint();
        cBest.print (cout, false);
    }

    if (_nDisplayLevel >= 3)
        _pGssIterator->printDirections (" Initial Directions");

    return;
}


//----------------------------------------------------------------------
//  Method postProcess
//----------------------------------------------------------------------
void  CitizenGSS::postProcess (void)
{
    // Print the final solution
    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        if (_pGssIterator->isFinished())
        {
            cout << " GSS GssIterator complete: ";
            _pGssIterator->printStopReason();
            cout << endl;
        }
        else
        {
            cout << " GSS GssIterator did not complete" << endl;
        }

        cout << "  Evaluated points from this citizen = "
             << _pGssIterator->getNumGssEvals() << endl;
        cout << endl;

        const GssPoint &  cBest = _pGssIterator->getBestPoint();
        if (cBest.getState() != DataPoint::UNEVALUATED)
        {
            cout << " GSS best point found:" << endl;
            cBest.print (cout, false);
            if (_cProbDef.hasNonlinearConstr())
            {
                cout.setf (ios::scientific);
                cout << "  F + p|C| = " << setprecision (Print::getPrecision())
                     << cBest.getBestF() << endl;
                cout.unsetf (ios::scientific);
            }
        }
    }

    if (_pParent != NULL)
    {
        //---- THE CITIZEN WAS SOLVING A SUBPROBLEM.
        //---- TELL THE PARENT IT IS FINISHED AND RETURN THE BEST POINT.
        //---- (ALL GSS PARENT CITIZENS USE THE SAME RETURN CODES.)
        GssChildReturnCodesType  nReturnCode = REASON_ERROR;
        if (_bWillExitEarly)
            nReturnCode = REASON_HALTED_BY_MEDIATOR;
        else if (_pGssIterator->hasStoppedAndConverged())
            nReturnCode = REASON_CONVERGED;
        else if (_pGssIterator->hasStoppedOutOfEvals())
            nReturnCode = REASON_NO_MORE_EVALS;
        else if (_pGssIterator->isFinished() == false)
            nReturnCode = REASON_HALTED_BY_MEDIATOR;
        _pParent->callbackFromChild (getIdNumber(),
                                     (int) nReturnCode,
                                     _pGssIterator->getBestPoint(),
                                     _pGssIterator->getNumGssEvals());
    }

    return;
}


//----------------------------------------------------------------------
//  Method exchange
//----------------------------------------------------------------------
void  CitizenGSS::exchange (const list< DataPoint * > &  cResultList,
                            const list< int >         &  cOwnerTags,
                                  list< DataPoint * > &  cWaitList)
{
    _cExchangeList.copyFrom (cResultList, *_pPenalty, cOwnerTags);
    printPreDiagnostics_();
    popBestInfeasiblePoints_();

    bool  bFoundNewBest
        = _pGssIterator->pointExchange (_cExchangeList,
                                        shouldIgnoreOtherPoints(),
                                        (_nDisplayLevel >= 3) );
    if (bFoundNewBest)
    {
        //---- GSS THEORY SAYS ALL OLDER POINTS CAN BE ERASED, BUT
        //---- PERFORMANCE MAY IMPROVE IF THEY ARE KEPT.  USER THRESHOLD
        //---- SAYS HOW MANY TO KEEP, AND HERE WE ARBITRARILY ERASE THE
        //---- LOWEST PRIORITY POINTS.
        while ((int) cWaitList.size() > _nMaxToKeepAfterNewBest)
        {
            DataPoint *  pTmp = cWaitList.front();
            delete pTmp;
            cWaitList.pop_front();
        }
    }
    _cExchangeList.copyTo (cWaitList);

    printPostDiagnostics_ (bFoundNewBest);
    _cExchangeList.prune();

    return;
}


//----------------------------------------------------------------------
//  Method setEarlyExit
//----------------------------------------------------------------------
void  CitizenGSS::setEarlyExit (void)
{
    _bWillExitEarly = true;
    return;
}


//----------------------------------------------------------------------
//  Method getIdNumber
//----------------------------------------------------------------------
int  CitizenGSS::getIdNumber (void) const
{
    return( _nID );
}


//----------------------------------------------------------------------
//  Method getName
//----------------------------------------------------------------------
const string &  CitizenGSS::getName (void) const
{
    return( _sName );
}


//----------------------------------------------------------------------
//  Method getState
//----------------------------------------------------------------------
Citizen::State  CitizenGSS::getState (void) const
{
    if (_pGssIterator->isFinished())
    {
        if (_pParent == NULL)
            return( FINISHED );
        else
            return( CHILD_FINISHED );
    }
    return( CONTINUE );
}


//----------------------------------------------------------------------
//  Private Method printPreDiagnostics
//----------------------------------------------------------------------
void  CitizenGSS::printPreDiagnostics_ (void) const
{
    if (_nDisplayLevel >= 2)
        _cExchangeList.print (" GSS result points received from Conveyor");
    return;
}


//----------------------------------------------------------------------
//  Private Method printPostDiagnostics
//----------------------------------------------------------------------
void  CitizenGSS::printPostDiagnostics_ (bool  bFoundNewBest) const
{
    if (bFoundNewBest && (_nDisplayLevel >= 2))
        cout << " GSS shifting to new best point." << endl;

    if (_nDisplayLevel >= 3)
        _pGssIterator->printDirections
            (" Directions after trial point generation");

    if (_nDisplayLevel >= 2)
        _cExchangeList.print (" GSS new trial points returned to Conveyor");  

    if (bFoundNewBest)
    {
        if (_nDisplayLevel >= 1)
        {
            const GssPoint &  cBest = _pGssIterator->getBestPoint();
            cout << " GSS New Best:" << endl;
            cBest.print (cout);
            if (_cProbDef.hasNonlinearConstr())
            {
                cout.setf (ios::scientific);
                cout << "  F + p|C| = " << setprecision (Print::getPrecision())
                     << cBest.getBestF() << endl;
                cout.unsetf (ios::scientific);
            }
        }
        if (_nDisplayLevel >= 3)
            _pGssIterator->printDirections (" New Directions");
    }

    if (_nDisplayLevel >= 1)
    {
        if (getState() == FINISHED)
            cout << " GSS state = FINISHED  - " << getName() << endl;
        else if (getState() == CHILD_FINISHED)
            cout << " GSS state = CHILD_FINISHED  - " << getName() << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method popBestInfeasiblePoints
//----------------------------------------------------------------------
void  CitizenGSS::popBestInfeasiblePoints_ (void)
{
    while (_cExchangeList.isEmpty() == false)
    {
        if (_cLinConstr.isFeasible (_cExchangeList.findBest()->getX()))
            break;

        if (_nDisplayLevel >= 2)
        {
            cout << " Popping off best point because it's linearly infeasible"
                 << ": Tag=" << _cExchangeList.findBest()->getTag() << endl;
        }

        GssPoint *  tmpPtr = _cExchangeList.popBest();
        delete tmpPtr;
    }
    return;
}


}     //-- namespace HOPSPACK
