// $Id: HOPSPACK_CitizenGssNlc.cpp 183 2010-12-15 18:22:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss-nlc/HOPSPACK_CitizenGssNlc.cpp $

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
  @file HOPSPACK_CitizenGssNlc.cpp
  @brief Implement HOPSPACK::CitizenGssNlc, subclass of Citizen.
*/

#include <iomanip>
#include <sstream>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CitizenGssNlc.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_gssChildReturnCodes.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Declare internal enumerations and constants
//----------------------------------------------------------------------

enum {
    FIN_REASON_CHILD_HALTED = 0,
    FIN_REASON_CHILD_ERROR,
    FIN_REASON_STEPTOL_SATISFIED,
    FIN_REASON_INFEAS_CANNOT_IMPROVE,
    FIN_REASON_MAX_EVALS,
    FIN_REASON_UNKNOWN
};


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
CitizenGssNlc::CitizenGssNlc (const int                         nIdNumber,
                              const string             &        sName,
                              const ParameterList      &        cParams,
                              const ProblemDef         &        cProbDef,
                              const LinConstr          &        cLinConstr,
                                    CallbackToMediator * const  pCallback,
                                    Citizen            * const  pParent)
    : Citizen (cParams, sName),
      _nID (nIdNumber),
      _sName (sName + " (GSS-NLC)"),
      _nState (CONTINUE),
      _cProbDef (cProbDef),
      _cLinConstr (cLinConstr),
      _cParentParams (cParams),
      _pCallbackToMed (pCallback),
      _pParent (pParent),
      _bWillExitEarly (false),
      _pChildParams (NULL),
      _pChildProbDef (NULL),
      _pLatestSubprobSol (NULL),
      _nTotalEvals (0),
      _nFinishedReason (FIN_REASON_UNKNOWN)
{
    //---- CHECK THAT THE PROBLEM TYPE IS APPROPRIATE FOR THIS CITIZEN.
    if (_cProbDef.isDomainContinuous() == false)
    {
        cerr << "ERROR: GSS-NLC citizen can only solve problems with"
             << " a continuous domain" << endl;
        throw "GSS-NLC Error";
    }

    //---- INITIALIZE PARAMETERS FROM THE INPUT LIST.
    if (extractParameters_ (_cParentParams, _cSubprobParams) == false)
    {
        throw "GSS-NLC Error";
    }

    _nM = _cProbDef.getNumNonlinearEqs() + _cProbDef.getNumNonlinearIneqs();

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
CitizenGssNlc::~CitizenGssNlc (void)
{
    if (_pChildParams != NULL)
        delete _pChildParams;
    if (_pChildProbDef != NULL)
        delete _pChildProbDef;
    if (_pLatestSubprobSol != NULL)
        delete _pLatestSubprobSol;

    return;
}


//----------------------------------------------------------------------
//  Method preProcess
//----------------------------------------------------------------------
void  CitizenGssNlc::preProcess (void)
{
    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        cout << "###################################################" << endl;
        cout << "###   HOPSPACK GSS-NLC Initialization Results   ###" << endl;
        cout << "###   Citizen name: " << getName() << endl;
        cout << endl;

        cout << "Priority = " << getPriority()
             << "  (1=highest, 10=lowest)" << endl;
        cout << endl;

        cout << "*** Parameter List (alphabetical order) ***" << endl;
        _cParentParams.print();
        cout << endl;

        _cProbDef.printDefinition (false);
        _cLinConstr.printDefinition (false);
        if (_cProbDef.hasNonlinearConstr())
            _cPenalty.printDefinition();

        cout << "### End HOPSPACK GSS-NLC Initialization Results ###" << endl;
        cout << "###################################################" << endl;
    }

    //---- PREPARE PARAMETERS FOR THE FIRST SUBPROBLEM.
    if (_pChildParams != NULL)
        delete _pChildParams;
    _pChildParams = new ParameterList (_cSubprobParams);
    _pChildParams->setParameter ("Type", "GSS-child");
    _pChildParams->setParameter ("Display", _nDisplaySubLevel);
    _pChildParams->setParameter ("Ignore Other Points", _bIgnoreOtherPoints);
    _pChildParams->setParameter ("Step Tolerance", _dCurrentStepTol);
    int  nEvalsCap = _nMaxSubprobEvals;
    if (_nMaxGssNlcEvals != -1)
    {
        nEvalsCap = max (_nMaxGssNlcEvals - _nTotalEvals, 0);
        if (_nMaxSubprobEvals != -1)
            nEvalsCap = min (nEvalsCap, _nMaxSubprobEvals);
    }
    _pChildParams->setParameter ("Maximum Evaluations", nEvalsCap);

    //---- ANY ERROR WITH THE FIRST SUBPROBLEM IS FATAL TO THE CITIZEN.
    _nChildID =  createNewChildCitizen_ (*_pChildParams, _cProbDef, _cPenalty);
    if (_nChildID < 0)
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
void  CitizenGssNlc::postProcess (void)
{
    if (_nState == CHILD_WAITING)
        _nState = FINISHED;

    // Print the final solution
    if (_nDisplayLevel >= 1)
    {
        cout << endl;
        if (_nState == FINISHED)
        {
            cout << " GSS-NLC complete: ";
            if (_nFinishedReason == FIN_REASON_STEPTOL_SATISFIED)
                cout << "Converged - step length smaller than tolerance";
            else if (_nFinishedReason == FIN_REASON_MAX_EVALS)
                cout << "Reached the evaluation limit for this citizen";
            else if (_nFinishedReason == FIN_REASON_CHILD_HALTED)
                cout << "Could not proceed after subproblem halted";
            else if (_nFinishedReason == FIN_REASON_CHILD_ERROR)
                cout << "Could not proceed after subproblem error";
            else if (_nFinishedReason == FIN_REASON_INFEAS_CANNOT_IMPROVE)
            {
                cout << "Best point is infeasible, cannot be improved."
                     << endl;
                cout << "  The problem itself may be infeasible"
                     << " (constraints impossible to satisfy)." << endl;
                cout << "  If the problem is believed to be feasible,"
                     << " then try one of the following:" << endl;
                cout << "  - increase 'Nonlinear Active Tolerance'"
                     << " in sublist 'Problem Definition'" << endl;
                cout << "  - reduce   'Step Tolerance'" << endl;
                cout << "  - increase 'Penalty Parameter Maximum'" << endl;
            }
            else
                cout << "Unknown reason!";
            cout << endl;
        }
        else
        {
            cout << " GSS-NLC did not complete" << endl;
        }

        cout << "  Evaluated points from this citizen and its children = "
             << _nTotalEvals << endl;

        if (   (_pLatestSubprobSol != NULL)
            && (_pLatestSubprobSol->getState() != DataPoint::UNEVALUATED) )
        {
            cout << " GSS-NLC most recent subproblem solution:" << endl;
            printPointWithPen_ (*_pLatestSubprobSol);
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
        else if (_nFinishedReason == FIN_REASON_STEPTOL_SATISFIED)
            nReturnCode = REASON_CONVERGED;
        else if (_nFinishedReason == FIN_REASON_INFEAS_CANNOT_IMPROVE)
            nReturnCode = REASON_CONVERGED;
        else if (_nFinishedReason == FIN_REASON_CHILD_HALTED)
            nReturnCode = REASON_HALTED_BY_MEDIATOR;
        else if (_nFinishedReason == FIN_REASON_MAX_EVALS)
            nReturnCode = REASON_NO_MORE_EVALS;
        _pParent->callbackFromChild (getIdNumber(),
                                     (int) nReturnCode,
                                     *_pLatestSubprobSol,
                                     _nTotalEvals);
    }

    return;
}


//----------------------------------------------------------------------
//  Method exchange
//----------------------------------------------------------------------
void  CitizenGssNlc::exchange (const list< DataPoint * > &  cResultList,
                               const list< int >         &  cOwnerTags,
                                     list< DataPoint * > &  cWaitList)
{
    //---- NOTHING TO DO.
    return;
}


//----------------------------------------------------------------------
//  Method setEarlyExit
//----------------------------------------------------------------------
void  CitizenGssNlc::setEarlyExit (void)
{
    _bWillExitEarly = true;
    return;
}


//----------------------------------------------------------------------
//  Method getIdNumber
//----------------------------------------------------------------------
int  CitizenGssNlc::getIdNumber (void) const
{
    return( _nID );
}


//----------------------------------------------------------------------
//  Method getName
//----------------------------------------------------------------------
const string &  CitizenGssNlc::getName (void) const
{
    return( _sName );
}


//----------------------------------------------------------------------
//  Method getState
//----------------------------------------------------------------------
Citizen::State  CitizenGssNlc::getState (void) const
{
    if ((_nState == FINISHED) && (_pParent != NULL))
        return( CHILD_FINISHED );

    return( _nState );
}


//----------------------------------------------------------------------
//  Method callbackFromChild
//----------------------------------------------------------------------
void  CitizenGssNlc::callbackFromChild (const int          nIdNumber,
                                        const int          nReturnCode,
                                        const DataPoint &  cFinalPoint,
                                        const int          nTotalEvals)
{
    if (nIdNumber != _nChildID)
    {
        cerr << "ERROR: Parent citizen '" << getName() << "' received callback"
             << " from wrong child = " << nIdNumber << endl;
        _nState = FINISHED;
        _nFinishedReason = FIN_REASON_CHILD_ERROR;
        return;
    }
    _nChildID = -1;

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
        cerr << "ERROR: Child GSS citizen failed to solve subproblem" << endl;
        cerr << "       Parent '" << getName() << "' is halting" << endl;
        _nState = FINISHED;
        _nFinishedReason = FIN_REASON_CHILD_ERROR;
        return;
    }
    if (_nDisplayLevel >= 2)
    {
        cout << " GSS-NLC subproblem solution:" << endl;
        printPointWithPen_ (cFinalPoint);
    }

    //---- TEST IF THE PROBLEM IS SOLVED OR CAN MAKE NO FURTHER PROGRESS.
    bool  bStop = isTimeToStop_ (nCode, cFinalPoint);
    if (_pLatestSubprobSol != NULL)
        delete _pLatestSubprobSol;
    _pLatestSubprobSol = new DataPoint (cFinalPoint);
    if (bStop)
    {
        if (_nDisplayLevel >= 1)
            cout << " GSS-NLC '" << getName() << "' is finished" << endl;

        //---- IF THIS IS A CHILD CITIZEN, THEN WAIT FOR THE MEDIATOR TO
        //---- CALL postProcess(), THE RIGHT PLACE TO CALLBACK TO THE PARENT.
        if (_pParent == NULL)
            _nState = FINISHED;
        else
            _nState = CHILD_WAITING;
        return;
    }

    updatePenalty_ (cFinalPoint);

    //---- UPDATE THE SUBPROBLEM STEP TOLERANCE.
    _dCurrentStepTol = max (_dCurrentStepTol * _dStepTolDecrease,
                            _dFinalStepTol);

    //---- PREPARE FOR THE NEXT SUBPROBLEM.
    if (_pChildParams != NULL)
        delete _pChildParams;
    _pChildParams = new ParameterList (_cSubprobParams);
    _pChildParams->setParameter ("Type", "GSS-child");
    _pChildParams->setParameter ("Display", _nDisplaySubLevel);
    _pChildParams->setParameter ("Ignore Other Points", _bIgnoreOtherPoints);
    _pChildParams->setParameter ("Step Tolerance", _dCurrentStepTol);
    int  nEvalsCap = _nMaxSubprobEvals;
    if (_nMaxGssNlcEvals != -1)
    {
        nEvalsCap = max (_nMaxGssNlcEvals - _nTotalEvals, 0);
        if (_nMaxSubprobEvals != -1)
            nEvalsCap = min (nEvalsCap, _nMaxSubprobEvals);
    }
    _pChildParams->setParameter ("Maximum Evaluations", nEvalsCap);

    if (_pChildProbDef != NULL)
        delete _pChildProbDef;
    _pChildProbDef = new ProblemDef (_cProbDef);
    _pChildProbDef->resetInitialX (cFinalPoint.getX(),
                                   cFinalPoint.getVecF(),
                                   cFinalPoint.getEqs(),
                                   cFinalPoint.getIneqs());

    //---- START THE NEXT SUBPROBLEM.
    _nChildID = createNewChildCitizen_ (*_pChildParams,
                                        *_pChildProbDef,
                                        _cPenalty);
    if (_nChildID < 0)
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
bool  CitizenGssNlc::extractParameters_ (ParameterList &  cParams,
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

    //---- DEFINE A PENALTY FUNCTION FOR THE NONLINEAR CONSTRAINTS.
    string  sPenFn
        = cParams.getOrSetParameter ("Penalty Function", "L2 Squared");
    _dInitialPenalty
        = cParams.getOrSetParameter ("Penalty Parameter", 1.0);
    if (_dInitialPenalty < 0.0)
    {
        cerr << "ERROR: Invalid 'Penalty Parameter' value in sublist 'GSS-NLC'"
             << endl;
        return( false );
    }
    _dInitialSmoothing
        = cParams.getOrSetParameter ("Penalty Smoothing Value", 1.0);
    if ((_dInitialSmoothing < 0.0) || (_dInitialSmoothing > 1.0))
    {
        cerr << "ERROR: Invalid 'Penalty Smoothing Value' in sublist 'GSS-NLC'"
             << endl;
        return( false );
    }
    if (_cPenalty.defineFunction (sPenFn,
                                  _dInitialPenalty,
                                  _dInitialSmoothing) == false)
    {
        cerr << "ERROR: GSS-NLC citizen could not construct penalty"
             << " function for nonlinear constraints" << endl;
        return( false );
    }
    cRemainder.deleteParameter ("Penalty Function");
    cRemainder.deleteParameter ("Penalty Parameter");
    cRemainder.deleteParameter ("Penalty Smoothing Value");

    //---- PENALTY PARAMETER ADJUSTMENTS.
    _dMaxPenalty
        = cParams.getOrSetParameter ("Penalty Parameter Maximum", 1.0e+8);
    if (_dMaxPenalty < _dInitialPenalty)
    {
        _dMaxPenalty = 1.0e+8;
        cerr << "WARNING: Cannot have 'Penalty Parameter Maximum' greater"
             << " than 'Penalty Parameter' in sublist 'GSS-NLC'"
             << endl;
        cerr << "         Changing 'Penalty Parameter Maximum' to "
             << _dMaxPenalty << endl;
    }
    _dPenaltyIncrease
        = cParams.getOrSetParameter ("Penalty Parameter Increase", 2.0);
    if (_dPenaltyIncrease <= 1.0)
    {
        _dPenaltyIncrease = 2.0;
        cerr << "WARNING: Must have 'Penalty Parameter Increase' greater"
             << " than 1 in sublist 'GSS-NLC'"
             << endl;
        cerr << "         Changing 'Penalty Parameter Increase' to "
             << _dPenaltyIncrease << endl;
    }
    cRemainder.deleteParameter ("Penalty Parameter Maximum");
    cRemainder.deleteParameter ("Penalty Parameter Increase");

    //---- SMOOTHING UPDATE PARAMETERS.
    if (_cPenalty.getSmoothing() > 0.0)
    {
        _dMinSmoothing
            = cParams.getOrSetParameter ("Smoothing Value Minimum", 1.0e-5);
        if (_dMinSmoothing < 0.0)
        {
            _dMinSmoothing = 0.0;
            cerr << "WARNING: Invalid negative value for"
                 << " 'Smoothing Value Minimum' in sublist 'GSS-NLC'" << endl;
            cerr << "         Changing 'Smoothing Value Minimum' to "
                 << _dMinSmoothing << endl;
        }

        _dSmoothingDecrease
            = cParams.getOrSetParameter ("Smoothing Value Decrease", 0.5);
        if ((_dSmoothingDecrease < 0.0) || (_dSmoothingDecrease > 1.0))
        {
            _dSmoothingDecrease = 0.5;
            cerr << "WARNING: Invalid value for"
                 << " 'Smoothing Value Decrease' in sublist 'GSS-NLC'" << endl;
            cerr << "         Changing 'Smoothing Value Decrease' to "
                 << _dSmoothingDecrease << endl;
        }
    }

    //---- MAXIMUM EVALUATIONS FOR NLC SOLVE, AND EVALUATIONS PER
    //---- LINEAR SUBPROBLEM SOLVER.
    _nMaxGssNlcEvals = cParams.getParameter ("Maximum Evaluations", -1);
    if (_nMaxGssNlcEvals < -1)
        _nMaxGssNlcEvals = -1;

    _nMaxSubprobEvals
        = cParams.getOrSetParameter ("Max Subproblem Evaluations", -1);
    if (_nMaxSubprobEvals < -1)
        _nMaxSubprobEvals = -1;
    if (_nMaxSubprobEvals > _nMaxGssNlcEvals)
    {
        _nMaxSubprobEvals = _nMaxGssNlcEvals;
        cerr << "WARNING: Value of 'Max Subproblem Evaluations'"
             << " exceeds 'Maximum Evaluations' in sublist 'GSS-NLC'" << endl;
        cerr << "         Changing 'Max Subproblem Evaluations' to "
             << _nMaxSubprobEvals << endl;
    }

    _bIgnoreOtherPoints
        = cParams.getOrSetParameter ("Ignore Other Points", true);

    //---- STEP TOLERANCE UPDATE PARAMETERS.
    _dFinalStepTol
        = cParams.getOrSetParameter ("Final Step Tolerance", 0.001);
    if (_dFinalStepTol <= 0.0)
    {
        cerr << "ERROR: 'Final Step Tolerance' must be positive"
             << " in sublist 'GSS-NLC'." << endl;
        throw "GSS-NLC Error";
    }
    _dCurrentStepTol
        = cParams.getOrSetParameter ("Initial Step Tolerance", 0.1);
    if (_dCurrentStepTol < _dFinalStepTol)
    {
        _dCurrentStepTol = _dFinalStepTol;
        cerr << "WARNING: 'Initial Step Tolerance' cannot be smaller than"
             << " 'Final Step Tolerance'" << endl;
        cerr << "         Changing 'Initial Step Tolerance' to "
             << _dCurrentStepTol << endl;
    }
    _dStepTolDecrease
        = cParams.getOrSetParameter ("Step Tolerance Decrease", 0.5);
    if ((_dStepTolDecrease < 0.0) || (_dStepTolDecrease > 1.0))
    {
        _dStepTolDecrease = 0.5;
        cerr << "WARNING: Invalid value for"
             << " 'Step Tolerance Decrease' in sublist 'GSS-NLC'" << endl;
        cerr << "         Changing 'Step Tolerance Decrease' to "
             << _dStepTolDecrease << endl;
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method createNewChildCitizen_
//----------------------------------------------------------------------
int  CitizenGssNlc::createNewChildCitizen_
         (      ParameterList     &  cChildParams,
          const ProblemDef        &  cProbDef,
          const NonlConstrPenalty &  cPenalty)
{
    //---- PENALTY FUNCTION INFORMATION HAS TO BE PASSED TO THE CHILD VIA
    //---- THE INPUT PARAMETERS, BUT A COPY IS NEEDED FOR THIS CITIZEN.
    //---- TO KEEP THEM IN SYNC, GENERATE THE INPUT PARAMETERS FROM THE
    //---- CURRENT INSTANCE.
    cChildParams.setParameter
        ("Penalty Function", cPenalty.getPenaltyName());
    cChildParams.setParameter
        ("Penalty Parameter", cPenalty.getCoefficient());
    cChildParams.setParameter
        ("Penalty Smoothing Value", cPenalty.getSmoothing());

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
                                       cChildParams,
                                       cProbDef,
                                       _cLinConstr,
                                       NULL,
                                       this);
    }
    catch (const char * const)
    {
        nReturnID = -1;
    }
    if ((nReturnID == -1) || (pChild == NULL))
    {
        cerr << "ERROR: Failed to create child GSS citizen" << endl;
        return( -1 );
    }

    //---- ADD THE CHILD TO THE MEDIATOR, WHICH WILL START RUNNING IT.
    if (_pCallbackToMed->addChildCitizen (pChild, getIdNumber()) == false)
    {
        cerr << "ERROR: Failed to add child GSS citizen" << endl;
        return( -2 );
    }

    if (_nDisplayLevel >= 2)
        cout << " CitizenGssNlc started child citizen " << nReturnID
             << endl << endl;

    return( nReturnID );
}


//---------------------------------------------------------------------
//  Private Method isTimeToStop_
//----------------------------------------------------------------------
bool  CitizenGssNlc::isTimeToStop_ (const GssChildReturnCodesType  nReturnCode,
                                    const DataPoint &              cTestPoint)
{
    if (nReturnCode == REASON_HALTED_BY_MEDIATOR)
    {
        _nFinishedReason = FIN_REASON_CHILD_HALTED;
        return( true );
    }

    bool  bIsLinFeas = (   _cProbDef.isBndsFeasible (cTestPoint.getX())
                        && _cLinConstr.isFeasible (cTestPoint.getX()) );
    bool  bIsNonlFeas = _cProbDef.isNonlinearlyFeasible (cTestPoint.getEqs(),
                                                         cTestPoint.getIneqs());

    if (   (nReturnCode == REASON_CONVERGED)
        && (_dCurrentStepTol <= _dFinalStepTol) )
    {
        if (bIsLinFeas && bIsNonlFeas)
        {
            _nFinishedReason = FIN_REASON_STEPTOL_SATISFIED;
            return( true );
        }
    }

    if ((_nMaxGssNlcEvals != -1) && (_nTotalEvals >= _nMaxGssNlcEvals))
    {
        _nFinishedReason = FIN_REASON_MAX_EVALS;
        return( true );
    }

    if (_pLatestSubprobSol != NULL)
    {
        //---- IF THE SUBPPROBLEM MADE NO CHANGE AND THE POINT IS INFEASIBLE
        //---- AND THE PENALTY PARAMETER CANNOT BE INCREASED, THEN WE CANNOT
        //---- MAKE ANY FURTHER IMPROVEMENTS IN FEASIBILITY.
        if (   _pLatestSubprobSol->isSamePoint (cTestPoint, 0.0)
            && (bIsNonlFeas == false)
            && (_cPenalty.getCoefficient() == _dMaxPenalty) )
        {
            _nFinishedReason = FIN_REASON_INFEAS_CANNOT_IMPROVE;
            return( true );
        }
    }

    return( false );
}


//---------------------------------------------------------------------
//  Private Method updatePenalty_
//----------------------------------------------------------------------
void  CitizenGssNlc::updatePenalty_ (const DataPoint &  cSubprobSolution)
{
    //---- DECIDE WHETHER TO INCREASE THE PENALTY PARAMETER.
    double  dThresh1 = _cProbDef.getNonlinearActiveTol();
    double  dThresh2 = _cPenalty.getSmoothing() * ((double) _nM) / 5.0;
    if (cSubprobSolution.getNonlConstrLInfNorm() > max (dThresh1, dThresh2))
    {
        double  dPenCoef = _cPenalty.getCoefficient();
        dPenCoef = min (_dPenaltyIncrease * dPenCoef, _dMaxPenalty);
        _cPenalty.updateCoefficient (dPenCoef);
    }

    //---- UPDATE ANY SMOOTHING OF THE PENALTY FUNCTION.
    double  dSmoothing = _cPenalty.getSmoothing();
    if (dSmoothing > 0.0)
    {
        double  dNew = max (dSmoothing * _dSmoothingDecrease, _dMinSmoothing);
        _cPenalty.updateSmoothing (dNew);
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method printPointWithPen_
//----------------------------------------------------------------------
void  CitizenGssNlc::printPointWithPen_ (const DataPoint &  cPoint) const
{
    cPoint.leftshift (cout, false);
    double  dPen = _cPenalty.computePenalty (cPoint.getEqs(),
                                             cPoint.getIneqs());
    double  dFplusPen = cPoint.getBestF()
                        + (cPoint.getPenaltySign() * dPen);
    cout.setf (ios::scientific);
    cout << ", p|C|=" << setprecision (Print::getPrecision())
         << dPen << endl;
    cout << "  F + p|C| = " << setprecision (Print::getPrecision())
         << dFplusPen << endl;
    cout.unsetf (ios::scientific);

    return;
}


}     //-- namespace HOPSPACK
