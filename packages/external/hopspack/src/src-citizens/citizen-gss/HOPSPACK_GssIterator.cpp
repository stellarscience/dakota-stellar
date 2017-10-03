// $Id: HOPSPACK_GssIterator.cpp 187 2011-04-06 22:39:43Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssIterator.cpp $

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
  \file HOPSPACK_GssIterator.cpp
  \brief Implementation of HOPSPACK::GssIterator
*/

#include <math.h>
#include <algorithm>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_GssIterator.hpp"
#include "HOPSPACK_GssDirections.hpp"
#include "HOPSPACK_GssList.hpp"
#include "HOPSPACK_GssPoint.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"

namespace HOPSPACK
{


GssIterator::GssIterator (const ProblemDef        &        cProbDef,
                          const LinConstr         &        cLinConstr,
                          const NonlConstrPenalty * const  pPenalty,
                                ParameterList     &        cParams)
    :
      _cParams (cParams),
      _cProbDef (cProbDef),
      _cPenalty (*pPenalty),
      _pBestPointPtr (NULL),
      _bIsFinished (false),
      _cDirections (cProbDef, cLinConstr, cParams),
      _cExchangeList(),
      _nState (NOT_CONVERGED),
      _cLinConstr (cLinConstr),
      _nNumGssEvals (0)
{
    if (   (_cParams.isParameterDouble ("Initial Step") == false)
        && (_cProbDef.isAutoScaled() == true) )
    {
        //---- AUTOSCALE MEANS THE PROBLEM HAS VARIABLE BOUNDS BUT NO SCALING.
        //---- IF AUTOSCALE DECIDED TO USE "1", THEN ENLARGE THE INITIAL STEP
        //---- TO ENCOMPASS THE BOUNDS.
        const Vector &  cScaling = _cProbDef.getVarScaling();
        bool  bIsUnity = true;
        for (int  i = 0; i < cScaling.size(); i++)
            if (cScaling[i] != 1.0)
                bIsUnity = false;
        if (bIsUnity == false)
            _cParams.setParameter ("Initial Step", 1.0);
        else
        {
            const Vector &  cUBnds = _cProbDef.getUpperBnds();
            const Vector &  cLBnds = _cProbDef.getLowerBnds();
            double  dMaxRange = 0.0;
            for (int  i = 0; i < cScaling.size(); i++)
                dMaxRange = max (dMaxRange, cUBnds[i] - cLBnds[i]);
            _cParams.setParameter ("Initial Step", dMaxRange / 2.0);
        }
        _dInitialStep = _cParams.getParameter ("Initial Step", 1.0);
    }
    else
    {
        _dInitialStep = _cParams.getOrSetParameter ("Initial Step", 1.0);
    }
    if (_dInitialStep <= 0.0)
    {
        _dInitialStep = _cParams.getDoubleParameter ("Step Tolerance");
        cerr << "WARNING: Invalid negative value for 'Initial Step'"
             << "  <GSS GssIterator>" << endl;
        cerr << "         Changing to " << _dInitialStep << endl;
    }

    //---- DEFAULT IS TRUE, WHICH WORKS BEST WHETHER SERIAL OR PARALLEL.
    _bUseRandomOrder = _cParams.getOrSetParameter("Use Random Order", true);

    //---- EVEN A PROBLEM WITH JUST SIMPLE BOUND CONSTRAINTS NEEDS LAPACK
    //---- ROUTINES TO COMPUTE THE SNAP.
    _bUseSnap = _cParams.getOrSetParameter ("Snap To Boundary", false);
    #if !defined(HAVE_LAPACK)
        if (_bUseSnap)
        {
            cerr << "WARNING: 'Snap To Boundary' cannot be used,"
                 << " need to build with LAPACK" << endl;
            cerr << "         Changing to 'false' and continuing" << endl;
            _bUseSnap = false;
        }
    #endif

    double  dDefault = _cParams.getDoubleParameter ("Step Tolerance") / 2.0;
    _dSnapDistance = _cParams.getOrSetParameter ("Snap Distance", dDefault);
    if (_dSnapDistance < 0.0)
    {
        cerr << "ERROR: Invalid negative value for 'Snap Distance'" << endl;
        cerr << "       <GssIterator::GssIterator>" << endl;
        throw "GSS Error";
    }

    _dSuffImprovementFactorAlpha
        = _cParams.getOrSetParameter ("Sufficient Improvement Factor", 0.01);
    if (_dSuffImprovementFactorAlpha < 0.0)
    {
        cerr << "WARNING: 'Sufficient Improvement Factor' in GSS sublist"
             << " cannot be negative" << endl;
        _dSuffImprovementFactorAlpha = 0.01;
        cerr << "         Changing to " << _dSuffImprovementFactorAlpha << endl;
    }

    _nMaxGssEvals = _cParams.getParameter ("Maximum Evaluations", -1);
    if (_nMaxGssEvals < -1)
        _nMaxGssEvals = -1;

    // Define initial point from user inputs or problem definition
    _pBestPointPtr = initializeBestPointPtr (_cProbDef, _cLinConstr);
    processNewBestPoint();

    // Push a copy of the initial point onto initialList.  If already
    // evaluated, then the cache will return its function value.
    _cInitialList.push (new GssPoint (*_pBestPointPtr));

    return;
}


GssIterator::~GssIterator (void)
{
    if (_pBestPointPtr != NULL)
        delete _pBestPointPtr;
    return;
}


bool  GssIterator::isFinished (void) const
{
    return _bIsFinished;
}


const GssPoint &  GssIterator::getBestPoint (void) const
{
    return( *_pBestPointPtr );
}


int  GssIterator::getNumGssEvals (void) const
{
    return( _nNumGssEvals );
}


double  GssIterator::getInitialStep (void) const
{
    return( _dInitialStep );
}


bool  GssIterator::hasStoppedAndConverged (void) const
{
    if (_bIsFinished == false)
        return( false );
    if ((_nState == STEPLENGTH_CONVERGED) || (_nState == OBJECTIVE_REACHED))
        return( true );
    return( false );
}


bool  GssIterator::hasStoppedOutOfEvals (void) const
{
    if (_bIsFinished == false)
        return( false );
    if (_nState == MAX_EVALS_FOR_CTZN)
        return( true );
    return( false );
}


void GssIterator::printInitializationInformation() const
{
    cout << "*** Parameter List (alphabetical order) ***" << endl;
    _cParams.print();

    return;
}


void GssIterator::printStopReason() const
{
    if (_nState == STEPLENGTH_CONVERGED)
        cout << "Converged - step length smaller than tolerance";
    else if (_nState == OBJECTIVE_REACHED)
        cout << "Converged - objective improved beyond target value";
    else if (_nState == MAX_EVALS_FOR_CTZN)
        cout << "Reached the evaluation limit for this citizen";
    else if (_nState == STOP_FROM_ERROR)
        cout << "Could not proceed after error";
    else
        cout << "Has not stopped yet";
    return;
}


void GssIterator::printDirections(const string label) const
{
    _cDirections.print(label);
    return;
}


bool GssIterator::pointExchange (      GssList &  newList,
                                 const bool       bShouldIgnoreOtherPoints,
                                 const bool       bPrintDetails)
{
    //---- TRANSFER POINTS FROM newList.  ALL POINTS ARE COPIES LOCAL TO
    //---- THIS CITIZEN, SO THEY MUST BE MOVED ON OR DELETED.
    GssPoint *  pNext = NULL;
    while ((pNext = newList.pop()) != NULL)
    {
        //---- CHECK IF THE POINT WAS SUBMITTED FROM THIS CITIZEN.
        vector< int >::iterator  cPos
            = find (_cPointsOwned.begin(), _cPointsOwned.end(), pNext->getTag());

        if (   (cPos != _cPointsOwned.end())
            && (pNext->getState() != DataPoint::EVALUATED_FC_FROM_CACHE) )
        {
            _nNumGssEvals++;
        }

        if (bShouldIgnoreOtherPoints == false)
            _cExchangeList.push (pNext);
        else
        {
            if (cPos != _cPointsOwned.end())
            {
                _cExchangeList.push (pNext);
                _cPointsOwned.erase (cPos);
            }
            else
            {
                delete pNext;
            }
        }
    }

    bool  bFoundNewBestPoint = false;
    if (_cExchangeList.isEmpty() == false)
        bFoundNewBestPoint = processEvaluatedTrialPoints();

    //---- ON THE FIRST ITERATION, A NEW BEST IS DECLARED IF THE INITIAL
    //---- POINT HAS A FUNCTION VALUE ATTACHED.
    if (_cInitialList.isEmpty() == false)
    {
        const GssPoint *  pBest = _cInitialList.findBest();
        if ((pBest != NULL) && (pBest->getBestF() != dne()))
        {
            bFoundNewBestPoint = true;
        }
    }

    if ((_nMaxGssEvals != -1) && (_nNumGssEvals > _nMaxGssEvals))
    {
        _nState = MAX_EVALS_FOR_CTZN;
        _bIsFinished = true;
        _cExchangeList.clearList();
    }
    else
    {
        generateTrialPoints (bPrintDetails);

        //---- TERMINATE THE SEARCH IF TRIAL POINTS ARE NOT GENERATED.
        //---- (BUG TICKET #3, 04/06/2011).
        if (   _cExchangeList.isEmpty()
            && _cDirections.isStepConverged()
            && (_pBestPointPtr != NULL) )
        {
            _nState = STEPLENGTH_CONVERGED;
            _bIsFinished = true;
        }

        newList.insertFromList (_cExchangeList);
        _cExchangeList.clearList();
        newList.copyTags (_cPointsOwned);
    }

    return bFoundNewBestPoint;
}


//----------------------------------------------------------------------
//  Private Methods
//----------------------------------------------------------------------

GssPoint *  GssIterator::initializeBestPointPtr
    (const ProblemDef &  cProbDef,
     const LinConstr  &  cLinConstr) const
{
    //---- DEFINE AN INITIAL X.
    Vector  initialF;
    Vector  initialEqs;
    Vector  initialIneqs;
    Vector  initialX = cProbDef.getInitialX();
    if (initialX.empty())
    {
        const Vector &  bndsLower = cProbDef.getLowerBnds();
        const Vector &  bndsUpper = cProbDef.getUpperBnds();
        initialX.resize (bndsLower.size());
        for (int  i = 0; i < initialX.size(); i++)
        {
            if (exists (bndsUpper[i]) && exists (bndsLower[i]))
                initialX[i] = (bndsUpper[i] + bndsLower[i]) / 2.0;
            else if (exists (bndsUpper[i]))
                initialX[i] = bndsUpper[i];
            else if (exists (bndsLower[i]))
                initialX[i] = bndsLower[i];
            else
                initialX[i] = 0;
        }

        if (cLinConstr.projectToFeasibility (initialX) == false)
        {
            cerr << "ERROR: Cannot generate initial point"  << endl;
            cerr << "       Cannot start GSS solver without"
                 << " a feasible initial point" << endl;
            cerr << "       <GssIterator::initializeBestPointPtr()>" << endl;
            throw "GSS Error";
        }
    }
    else
    {
        initialF = cProbDef.getInitialF();
        initialEqs = cProbDef.getInitialEqs();
        initialIneqs = cProbDef.getInitialIneqs();
    }

    // Assert that x is feasible.
    if (   (cProbDef.isBndsFeasible (initialX) == false)
        || (cLinConstr.isFeasible (initialX) == false) )
    {
        cerr << "ERROR: Infeasible initial point after correcting" << endl;
        cerr << "       Cannot start GSS solver without"
             << " a feasible initial point" << endl;
        cerr << "       <GssIterator::initializeBestPointPtr()>" << endl;
        throw "GSS Error";
    }

    // Pretend the point was generated with an initial step.
    GssPoint *  pResult = new GssPoint (cProbDef.getObjType(),
                                        _cPenalty,
                                        initialX,
                                        _dInitialStep,
                                        GssPoint::NO_PARENT_TAG,
                                        0.0, 0.0,
                                        0.0,
                                        GssPoint::NO_DIR_INDEX);
    if (   (initialF.empty() == false)
        || (initialEqs.empty() == false)
        || (initialIneqs.empty() == false) )
    {
        pResult->setEvalFC (initialF, initialEqs, initialIneqs,
                            "(User Initial Point)");
    }
    return pResult;
}


void  GssIterator::processNewBestPoint (GssPoint *  newBestPointPtr)
{
    // Update the best point
    if (newBestPointPtr != NULL)
    {
        delete _pBestPointPtr;
        _pBestPointPtr = newBestPointPtr;
    }

    // Note that the Mediator checks if the objective target is reached.

    // Update the search directions.
    _cDirections.computeNewDirections (*_pBestPointPtr);

    return;
}


bool GssIterator::processEvaluatedTrialPoints()
{
    bool  bFoundNewBestPoint = false;

    // Check for a new best point
    const GssPoint * const  pNewBest = _cExchangeList.findBest();
    if (   pNewBest->hasSufficientImprovement()
        && pNewBest->isBetterObjThan (*_pBestPointPtr) )
    {
        processNewBestPoint (_cExchangeList.popBest());
        bFoundNewBestPoint = true;

        //---- DESTROY ANY REMAINING POINTS.
        while (_cExchangeList.isEmpty() == false)
        {
            GssPoint *  ptr = _cExchangeList.pop();
            delete ptr;
        }
        _cExchangeList.clearList();
        _cPointsOwned.clear();
    }
    else
    {
        // Otherwise, just process the list
        GssPoint *  ptr;
        bool stepReduced = false;
        while (_cExchangeList.isEmpty() == false)
        {
            ptr = _cExchangeList.pop();
            if (ptr->getParentTag() == _pBestPointPtr->getTag())
            {
                _cDirections.reduceStep (ptr->getDirIndex());
                stepReduced = true;
            }
            delete ptr;
        }
        // Check for step length convergence.  If not converged, append
        // new directions to list.
        if (_cDirections.isStepConverged())
        {
            _nState = STEPLENGTH_CONVERGED;
            _bIsFinished = true;
        }
        else if (stepReduced)
            _cDirections.appendNewDirections();
    }

    return bFoundNewBestPoint;
}


void GssIterator::generateTrialPoints (const bool  bPrintDetails)
{
    // Add initial points if provided.  Happens only on first iteration.
    // This may redundantly insert a user-provided initial point, but
    // there is no harm -- the cache will return it without evaluating.
    if (_cInitialList.isEmpty() == false)
    {
        _cExchangeList.insertFromList (_cInitialList);
        _cInitialList.clearList();
    }

    // Local references
    const Vector& parentX = _pBestPointPtr->getX();

    Vector& x = _tmpReusableVector;
    int n = parentX.size();
    x.resize(n);

    vector< int > dirIndices;
    _cDirections.getDirectionIndices (dirIndices);

    if (_bUseRandomOrder == true)
    {
        // Randomize access order.
        random_shuffle (dirIndices.begin(), dirIndices.end());
    }

    for (int i = 0; i < (int) dirIndices.size(); i++)
    {
        int  idx = dirIndices[i];
        const Vector &  nextDirSet = _cDirections.getDirection (idx);
        const double  dNextStep = _cDirections.getStep (idx);

        double  dMaxStep = _cLinConstr.maxStep (parentX, nextDirSet, dNextStep);

        // Only add nontrivial points to the exchange list.
        if (dMaxStep <= 0)
        {
            // Cannot move along this direction and remain feasible.
            // This can happen if there are linear constraints and Active Tol
            // is too small, as then maxStep() returns zero.
            _cDirections.setStepConverged (idx);
        }
        else
        {
            for (int j = 0; j < n; j++)
                x[j] = parentX[j] + (dMaxStep * nextDirSet[j]);

            //---- SNAP TO THE BOUNDARY IF SO DIRECTED.
            bool  bNeedPrettySnapDetails = false;
            if (_bUseSnap)
            {
                Vector xsnap = x;
                _cLinConstr.snapToBoundary (xsnap, _dSnapDistance);
                if (bPrintDetails)
                {
                    double  dMaxDiff = 0.0;
                    for (int  k = 0; k < x.size(); k++)
                    {
                        double  dDiff = fabs (x[k] - xsnap[k]);
                        if (dDiff > dMaxDiff)
                            dMaxDiff = dDiff;
                    }
                    //---- PRINT THRESHOLD IS ARBITRARY, BUT snapToBoundary
                    //---- ROUTINELY ALTERS x BY A FEW MACHINE EPSILON
                    //---- EVEN WHEN ALREADY ON THE BOUNDARY.
                    if (dMaxDiff >= 1.0e-14)
                    {
                        cout << "    Snap moved point, |diff|_inf = "
                             << dMaxDiff;
                        bNeedPrettySnapDetails = true;
                    }
                }
                if (_cLinConstr.isFeasible (xsnap))
                    x = xsnap;
                else if (bPrintDetails)
                {
                    if (bNeedPrettySnapDetails)
                        cout << endl;
                    cout << "    Snap point ignored; linearly infeasible!"
                         << endl;
                }

                //! Josh: put in maxStep search here (note from Tammy)
            }

            if (_cProbDef.isBndsFeasible (x) == false)
            {
                //---- LINEAR CONSTRAINTS CAN CAUSE ROUNDOFF ERRORS THAT PUSH
                //---- THE TRIAL POINT OUTSIDE THE VARIABLE BOUNDS.
                //---- TRY MAKING A SMALL CORRECTION LIMITED BY THE
                //---- CONSTRAINT TOLERANCE.
                //----
                //---- AS AN EXAMPLE, CONSIDER AN EQUALITY CONSTRAINT ALMOST
                //---- PARALLEL TO A BOUND.  THE CONSTRAINT TOLERANCE ALLOWS
                //---- A CERTAIN AMOUNT OF DEVIATION FROM THE EQUALITY, AND
                //---- THUS VIOLATION OF THE BOUND.  THE CORRECTION IN THIS
                //---- CASE ALLOWS GSS TO CONTINUE; HOWEVER, STEPS CAN NOW
                //---- SLIDE AWAY FROM THE TRUE POINT BECAUSE OF THE ALLOWED
                //---- CONSTRAINT ERROR.  IF THE CONSTRAINT TOLERANCE IS
                //---- TIGHTENED, THEN AN INFEASIBLE STEP IS NEVER TAKEN
                //---- (AND THIS POINT IN CODE NEVER REACHED) BECAUSE maxStep()
                //---- WILL NOT MOVE AWAY FROM THE EQUALITY.
                double  dTol = _cLinConstr.getActiveTol();
                dTol = max (2.0 * dTol, getMachineEpsilon());
                _cProbDef.makeBndsFeasible (dTol, x);
            }

            if (   (_cProbDef.isBndsFeasible (x) == false)
                || (_cLinConstr.isFeasible (x, true) == false) )
            {
                cerr << "WARNING: GSS generated a point infeasible wi/re"
                     << " linear constraints" << endl;
                cerr << "         Cannot continue in this direction" << endl;
                _nState = STOP_FROM_ERROR;
                _bIsFinished = true;
            }
            else
            {
                // Define sufficient improvement amount.
                double  dSuffImpr = _dSuffImprovementFactorAlpha
                                    * (dNextStep * dNextStep);
                double  dPenaltyTerm
                    = _cPenalty.computePenalty (_pBestPointPtr->getEqs(),
                                                _pBestPointPtr->getIneqs());

                // Create a new trial point.
                GssPoint *  newPointPtr
                    = new GssPoint (_cProbDef.getObjType(),
                                    _cPenalty,
                                    x,
                                    dNextStep,
                                    _pBestPointPtr->getTag(),
                                    _pBestPointPtr->getBestF(),
                                    dPenaltyTerm,
                                    dSuffImpr,
                                    idx);
                if (bNeedPrettySnapDetails)
                    cout << ", created with tag = "
                         << newPointPtr->getTag() << endl;

                // Save off trial point information.
                _cDirections.setTrueStepAndTag (idx,
                                                dMaxStep,
                                                newPointPtr->getTag());
                // Push this trial point onto the new trial point list.
                // Ownership of the pointer transfers to the list.
                _cExchangeList.push (newPointPtr);
            }
        }
    }
    return;
}


}     //-- namespace HOPSPACK
