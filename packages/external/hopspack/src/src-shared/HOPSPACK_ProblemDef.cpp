// $Id: HOPSPACK_ProblemDef.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ProblemDef.cpp $

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
  @file HOPSPACK_ProblemDef.cpp
  @brief Implement HOPSPACK::ProblemDef
*/

#include <iomanip>
#include <math.h>      //-- FOR fabs

#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_ParameterList.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Static variables
//----------------------------------------------------------------------
static const string  sPROBDEF = "Problem Definition";


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ProblemDef::ProblemDef (void)
{
    _nNumObjs = 0;
    _nNumVars = 0;
    _nDisplayFlag = 0;
    _nNumNonlEqs = 0;
    _nNumNonlIneqs = 0;

    //---- THIS VALUE MATCHES THE DEFAULT FOR LINEAR CONSTRAINTS.
    _dNonlActTol = 1.0e-7;

    return;
}


//----------------------------------------------------------------------
//  Copy Constructor
//----------------------------------------------------------------------
ProblemDef::ProblemDef (const ProblemDef &  cArg)
{
    _nNumObjs = cArg._nNumObjs;
    _nObjType = cArg._nObjType;
    _dObjTarget = cArg._dObjTarget;
    _dObjTgtPercent = cArg._dObjTgtPercent;

    _nNumVars = cArg._nNumVars;

    _naVarTypes = cArg._naVarTypes;
    _cScaling = cArg._cScaling;
    _bIsAutoScaled = cArg._bIsAutoScaled;
    _cLoBnds = cArg._cLoBnds;
    _cUpBnds = cArg._cUpBnds;

    _cInitialX = cArg._cInitialX;
    _cInitialF = cArg._cInitialF;
    _cInitialEqs = cArg._cInitialEqs;
    _cInitialIneqs = cArg._cInitialIneqs;

    _nNumNonlEqs = cArg._nNumNonlEqs;
    _nNumNonlIneqs = cArg._nNumNonlIneqs;
    _dNonlActTol = cArg._dNonlActTol;

    _nDisplayFlag = cArg._nDisplayFlag;

    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ProblemDef::~ProblemDef (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method initialize
//----------------------------------------------------------------------
bool  ProblemDef::initialize (const ParameterList &  cProbDefParams)
{
    if (setupObj_ (cProbDefParams) == false)
        return( false );
    if (setupVars_ (cProbDefParams) == false)
        return( false );
    if (setupVarBndsAndScaling_ (cProbDefParams) == false)
        return( false );
    if (setupMisc_ (cProbDefParams) == false)
        return( false );
    if (setupInitialPoint_ (cProbDefParams) == false)
        return( false );

    return( true );
}


//----------------------------------------------------------------------
//  Accessors
//----------------------------------------------------------------------

int  ProblemDef::getNumObjs (void) const
{
    return( _nNumObjs );
}

ProblemDef::ObjectiveType  ProblemDef::getObjType (void) const
{
    return( _nObjType );
}

double  ProblemDef::getObjTarget (void) const
{
    return( _dObjTarget );
}

double  ProblemDef::getObjPercentErrorThreshold (void) const
{
    return( _dObjTgtPercent );
}

const vector< ProblemDef::VariableType >  ProblemDef::getVarTypes (void) const
{
    return( _naVarTypes );
}

bool  ProblemDef::isDomainContinuous (void) const
{
    for (int  i = 0; i < _nNumVars; i++)
    {
        if ((_naVarTypes[i] == INTEGER) || (_naVarTypes[i] == ORDINAL))
        {
            if (_cLoBnds[i] != _cUpBnds[i])
                return( false );
        }
        else if (_naVarTypes[i] != CONTINUOUS)
            return( false );
    }
    return( true );
}

const Vector &  ProblemDef::getVarScaling (void) const
{
    return( _cScaling );
}

bool  ProblemDef::isAutoScaled (void) const
{
    return( _bIsAutoScaled );
}

const Vector &  ProblemDef::getLowerBnds (void) const
{
    return( _cLoBnds );
}

const Vector &  ProblemDef::getUpperBnds (void) const
{
    return( _cUpBnds );
}

const Vector &  ProblemDef::getInitialX (void) const
{
    return( _cInitialX );
}

const Vector &  ProblemDef::getInitialF (void) const
{
    return( _cInitialF );
}

const Vector &  ProblemDef::getInitialEqs (void) const
{
    return( _cInitialEqs );
}

const Vector &  ProblemDef::getInitialIneqs (void) const
{
    return( _cInitialIneqs );
}

bool  ProblemDef::hasNonlinearConstr (void) const
{
    return( (_nNumNonlEqs > 0) || (_nNumNonlIneqs > 0) );
}

double  ProblemDef::getNonlinearActiveTol (void) const
{
    return( _dNonlActTol );
}

int  ProblemDef::getNumNonlinearEqs (void) const
{
    return( _nNumNonlEqs );
}

int  ProblemDef::getNumNonlinearIneqs (void) const
{
    return( _nNumNonlIneqs );
}


//----------------------------------------------------------------------
//  Method isBndsFeasible
//----------------------------------------------------------------------
bool  ProblemDef::isBndsFeasible (const Vector &  cX) const
{
    if (cX.size() != _nNumVars)
    {
        cerr << "ERROR: Bad argument length"
             << "  <ProblemDef::isBndsFeasible()>" << endl;
        throw INTERNAL_ERROR;
    }

    for (int  i = 0; i < _nNumVars; i++)
    {
        if (exists (_cLoBnds[i]) && (cX[i] < _cLoBnds[i]))
            return( false );
        if (exists (_cUpBnds[i]) && (cX[i] > _cUpBnds[i]))
            return( false );
    }

    return( true );
}


//----------------------------------------------------------------------
//  Method isNonlinearlyFeasible
//----------------------------------------------------------------------
bool  ProblemDef::isNonlinearlyFeasible (const Vector &  cEqs,
                                         const Vector &  cIneqs) const
{
    if (hasNonlinearConstr() == false)
        return( true );

    if ((cEqs.size() != _nNumNonlEqs) || (cIneqs.size() != _nNumNonlIneqs))
    {
        cerr << "ERROR: Bad argument length"
             << "  <ProblemDef::isNonlinearlyFeasible()>" << endl;
        throw INTERNAL_ERROR;
    }

    for (int  i = 0; i < cEqs.size(); i++)
        if (fabs (cEqs[i]) > _dNonlActTol)
            return( false );

    //---- INEQUALITY VALUES ARE NEGATIVE IF IN VIOLATION.
    for (int  i = 0; i < cIneqs.size(); i++)
        if (cIneqs[i] < (- _dNonlActTol))
            return( false );

    return( true );
}


//----------------------------------------------------------------------
//  Method makeBndsFeasible
//----------------------------------------------------------------------
bool  ProblemDef::makeBndsFeasible (const double    dTol,
                                          Vector &  cX) const
{
    if (cX.size() != _nNumVars)
    {
        cerr << "ERROR: Bad argument length"
             << "  <ProblemDef::makeBndsFeasible()>" << endl;
        throw INTERNAL_ERROR;
    }

    for (int  i = 0; i < _nNumVars; i++)
    {
        if (exists (_cLoBnds[i]) && (cX[i] < _cLoBnds[i]))
        {
            if ((dTol < 0.0) || (_cLoBnds[i] - cX[i] <= dTol))
                cX[i] = _cLoBnds[i];
            else
                return( false );
        }
        if (exists (_cUpBnds[i]) && (cX[i] > _cUpBnds[i]))
        {
            if ((dTol < 0.0) || (cX[i] - _cUpBnds[i] <= dTol))
                cX[i] = _cUpBnds[i];
            else
                return( false );
        }
    }

    return( true );
}


//----------------------------------------------------------------------
//  Method resetInitialX
//----------------------------------------------------------------------
void  ProblemDef::resetInitialX (const Vector &  newX)
{
    if ((newX.empty() == false) && (newX.size() != _nNumVars))
    {
        cerr << "ERROR: Bad argument length for newX"
             << "  <ProblemDef::resetInitialX()>" << endl;
        throw INTERNAL_ERROR;
    }

    _cInitialX = newX;
    _cInitialF.resize (0);
    _cInitialEqs.resize (0);
    _cInitialIneqs.resize (0);
    return;
}


//----------------------------------------------------------------------
//  Method resetInitialX (with function and constraint values)
//----------------------------------------------------------------------
void  ProblemDef::resetInitialX (const Vector &  newX,
                                 const Vector &  newF,
                                 const Vector &  newEqs,
                                 const Vector &  newIneqs)
{
    resetInitialX (newX);

    if (   (newF.size() != _nNumObjs)
        || (newEqs.size() != _nNumNonlEqs)
        || (newIneqs.size() != _nNumNonlIneqs) )
    {
        cerr << "ERROR: Bad argument length"
             << "  <ProblemDef::resetInitialX()>" << endl;
        throw INTERNAL_ERROR;
    }

    _cInitialF = newF;
    _cInitialEqs = newEqs;
    _cInitialIneqs = newIneqs;

    return;
}


//----------------------------------------------------------------------
//  Method isObjTargetReached
//----------------------------------------------------------------------
bool  ProblemDef::isObjTargetReached (const double    dObjValue,
                                            double &  dPercent) const
{
    if (exists (dObjValue) == false)
        return( false );
    if (exists (_dObjTarget) == false)
        return( false );

    if (_nObjType == FIND_FEASIBLE_PT)
    {
        dPercent = 0.0;
        return( true );
    }

    //---- CALCULATE THE DIFFERENCE BETWEEN dObjValue AND THE TARGET
    //---- AS A PERCENTAGE OF THE TARGET, POSITIVE IF TARGET NOT REACHED.
    double  dDiff = dObjValue - _dObjTarget;
    if (_nObjType == MAXIMIZE)
        dDiff = -dDiff;
    if (dDiff <= 0.0)
        dPercent = 0.0;
    else
    {
        double  dMag = fabs (_dObjTarget);
        dMag = max (1.0e-4, dMag);
        dPercent = 100.0 * (dDiff / dMag);
    }

    //---- ABSOLUTE TEST.
    if ((_nObjType == MINIMIZE) && (dObjValue <= _dObjTarget))
        return( true );
    if ((_nObjType == MAXIMIZE) && (dObjValue >= _dObjTarget))
        return( true );

    //---- PERCENTAGE TEST.
    if (exists (_dObjTgtPercent) && (dPercent <= _dObjTgtPercent))
        return( true );

    return( false );
}


//----------------------------------------------------------------------
//  Private Method setupObj_
//----------------------------------------------------------------------
bool  ProblemDef::setupObj_ (const ParameterList &  cParams)
{
    //---- OBJECTIVE USER PARAMETERS ARE OPTIONAL.
    //---- THE DEFAULT IS TO MINIMIZE ONE OBJECTIVE FUNCTION.

    _nNumObjs = cParams.getParameter ("Number Objectives", 1);
    if (_nNumObjs < 0)
    {
        cerr << "ERROR: Bad 'Number Objectives' value " << _nNumObjs
             << " in '" << sPROBDEF << "' sublist" << endl;
        return( false );
    }
    if (_nNumObjs == 0)
    {
        cerr << "ERROR: Currently do not support 'Number Objectives' = 0"
             << " in '" << sPROBDEF << "' sublist" << endl;
        return( false );
    }
    if (_nNumObjs > 1)
    {
        cerr << "ERROR: Currently do not support 'Number Objectives' > 1"
             << " in '" << sPROBDEF << "' sublist" << endl;
        return( false );
    }

    if (_nNumObjs == 1)
    {
        string  cType = cParams.getParameter ("Objective Type", "Minimize");
        if (cType == "Minimize")
            _nObjType = MINIMIZE;
        else if (cType == "Maximize")
            _nObjType = MAXIMIZE;
        else
        {
            cerr << "ERROR: Unknown 'Objective Type' " << cType
                 << " in '" << sPROBDEF << "' sublist" << endl;
            return( false );
        }
    }

    _dObjTarget = cParams.getParameter ("Objective Target", dne());

    _dObjTgtPercent = cParams.getParameter ("Objective Percent Error", dne());
    if (exists (_dObjTgtPercent) && (exists (_dObjTarget) == false))
    {
        cerr << "WARNING: Cannot define 'Objective Percent Error' without"
             << " also defining 'Objective Target'" << endl;
        cerr << "         Ignoring 'Objective Percent Error'"
             << " in '" << sPROBDEF << "' sublist" << endl;
        _dObjTgtPercent = dne();
    }
    if (exists (_dObjTgtPercent) && (_dObjTgtPercent < 0.0))
    {
        cerr << "WARNING: Cannot make 'Objective Percent Error' less than zero"
             << endl;
        cerr << "         Changing 'Objective Percent Error' to zero"
             << " in '" << sPROBDEF << "' sublist" << endl;
        _dObjTgtPercent = 0.0;
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method setupVars_
//----------------------------------------------------------------------
bool  ProblemDef::setupVars_ (const ParameterList &  cParams)
{
    _nNumVars = 0;
    if (cParams.isParameterInt ("Number Unknowns") == false)
    {
        cerr << "ERROR: Need 'Number Unknowns'"
             << " in '" << sPROBDEF << "' sublist" << endl;
        return( false );
    }
    _nNumVars = cParams.getParameter ("Number Unknowns", 0);
    if (_nNumVars <= 0)
    {
        cerr << "ERROR: Bad 'Number Unknowns' = " << _nNumVars
             << " in '" << sPROBDEF << "' sublist" << endl;
        return( false );
    }

    if (cParams.isParameterCharVec ("Variable Types"))
    {
        vector< char >  caTmp = cParams.getCharVecParameter ("Variable Types");
        if ((int) caTmp.size() != _nNumVars)
        {
            cerr << "ERROR: Length of 'Variable Types' = " << caTmp.size()
                 << " does not match 'Number Unknowns' = " << _nNumVars << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            return( false );
        }
        _naVarTypes.resize (_nNumVars);
        for (int  i = 0; i < _nNumVars; i++)
        {
            if ((caTmp[i] == 'C') || (caTmp[i] == 'c'))
                _naVarTypes[i] = CONTINUOUS;
            else if ((caTmp[i] == 'I') || (caTmp[i] == 'i'))
                _naVarTypes[i] = INTEGER;
            else if ((caTmp[i] == 'O') || (caTmp[i] == 'o'))
                _naVarTypes[i] = ORDINAL;
            else
            {
                cerr << "ERROR: Unknown variable type '" << caTmp[i]
                     << "' for element [" << (i + 1) << "]" << endl;
                cerr << "       See 'Variable Types'"
                     << " in '" << sPROBDEF << "' sublist" << endl;
                return( false );
            }
        }
    }
    else
    {
        _naVarTypes.assign (_nNumVars, CONTINUOUS);
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method setupVarBndsAndScaling_
//----------------------------------------------------------------------
bool  ProblemDef::setupVarBndsAndScaling_ (const ParameterList &  cParams)
{
    if (cParams.isParameterVector ("Lower Bounds"))
    {
        _cLoBnds = cParams.getVectorParameter ("Lower Bounds");
        if (_cLoBnds.size() != _nNumVars)
        {
            cerr << "ERROR: Length of 'Lower Bounds' = " << _cLoBnds.size()
                 << " does not match 'Number Unknowns' = " << _nNumVars << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            return( false );
        }
    }
    if (cParams.isParameterVector ("Upper Bounds"))
    {
        _cUpBnds = cParams.getVectorParameter ("Upper Bounds");
        if (_cUpBnds.size() != _nNumVars)
        {
            cerr << "ERROR: Length of 'Upper Bounds' = " << _cLoBnds.size()
                 << " does not match 'Number Unknowns' = " << _nNumVars << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            return( false );
        }
    }
    if ((_cLoBnds.empty() == false) && (_cUpBnds.empty() == false))
    {
        for (int  i = 0; i < _nNumVars; i++)
        {
            if (   exists (_cLoBnds[i]) && exists (_cUpBnds[i])
                && (_cUpBnds[i] < _cLoBnds[i]) )
            {
                cerr << "ERROR: Variable bounds are inconsistent for element ["
                     << (i + 1) << "]" << endl;
                cerr << "       See 'Lower Bounds' and 'Upper Bounds'"
                     << " in '" << sPROBDEF << "' sublist" << endl;
                return( false );
            }
        }
    }

    _bIsAutoScaled = false;

    if (cParams.isParameterVector ("Scaling"))
    {
        _cScaling = cParams.getVectorParameter ("Scaling");
        if (_cScaling.size() != _nNumVars)
        {
            cerr << "ERROR: Length of 'Scaling' = " << _cScaling.size()
                 << " does not match 'Number Unknowns' = " << _nNumVars << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            return( false );
        }
    }

    //---- EITHER SCALING OR VARIABLE BOUNDS MUST BE PROVIDED SO THAT CITIZENS
    //---- HAVE SOME WAY TO SCALE STEP LENGTHS.
    if (_cScaling.empty())
    {
        if (_cLoBnds.empty() || _cUpBnds.empty())
        {
            cerr << "ERROR: Must define 'Scaling' or all variable bounds"
                 << " in '" << sPROBDEF << "' sublist" << endl;
            return( false );
        }
        for (int  i = 0; i < _nNumVars; i++)
        {
            if (   (exists (_cLoBnds[i]) == false)
                || (exists (_cUpBnds[i]) == false) )
            {
                cerr << "ERROR: Must define 'Scaling' or"
                     << " all variable bounds"
                     << " in '" << sPROBDEF << "' sublist" << endl;
                cerr << "  One or both bounds are undefined for element ["
                     << (i + 1) << "]" << endl;
                return( false );
            }
        }

        //---- AUTOMATIC SCALING COULD BE THE WIDTH OF THE VARIABLE BOUNDS,
        //---- BUT PERFORMS POORLY UNLESS BOUNDS ARE MEANINGFUL.
        for (int  i = 0; i < _nNumVars; i++)
        {
            _cScaling.push_back (1.0);
        }
        _bIsAutoScaled = true;
    }
    for (int  i = 0; i < _nNumVars; i++)
    {
        if (_cScaling[i] <= 0.0)
        {
            cerr << "ERROR: 'Scaling' must be positive, element ["
                 << (i + 1) << "] is not" << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            return( false );
        }
    }

    //---- IF VARIABLE BOUNDS ARE NOT PROVIDED, THEN SET THEM AS UNBOUNDED.
    if (_cLoBnds.empty())
        _cLoBnds.assign (_nNumVars, dne());

    if (_cUpBnds.empty())
        _cUpBnds.assign (_nNumVars, dne());

    return( true );
}


//----------------------------------------------------------------------
//  Private Method setupInitialPoint_
//----------------------------------------------------------------------
bool  ProblemDef::setupInitialPoint_ (const ParameterList &  cParams)
{
    //---- INITIAL POINT USER PARAMETERS ARE OPTIONAL.
    //---- THE DEFAULT IS NONE.

    bool  bModifiedInitialX = false;
    if (cParams.isParameterVector ("Initial X"))
    {
        _cInitialX = cParams.getVectorParameter ("Initial X");
        if (_cInitialX.size() != _nNumVars)
        {
            cerr << "ERROR: Length of 'Initial X' = " << _cInitialX.size()
                 << " does not match 'Number Unknowns' = " << _nNumVars << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            _cInitialX.resize (0);
            return( false );
        }
        for (int  i = 0; i < _nNumVars; i++)
        {
            if (exists (_cInitialX[i]) == false)
            {
                cerr << "ERROR: Element [" << (i + 1) << "]"
                     << " of 'Initial X' is undefined"
                     << " in '" << sPROBDEF << "' sublist" << endl;
                _cInitialX.resize (0);
                return( false );
            }
        }
        if (isBndsFeasible (_cInitialX) == false)
        {
            cerr << "WARNING: The point 'Initial X' violates"
                 << " an upper or lower variable bound" << endl;
            cerr << "         Modifying 'Initial X' to be feasible"
                 << " in '" << sPROBDEF << "' sublist" << endl;
            makeBndsFeasible (-1.0, _cInitialX);
            bModifiedInitialX = true;
        }
    }

    if (cParams.isParameterVector ("Initial F"))
    {
        _cInitialF = cParams.getVectorParameter ("Initial F");
        if (_cInitialF.size() != _nNumObjs)
        {
            cerr << "ERROR: Length of 'Initial F' = " << _cInitialF.size()
                 << " does not match 'Number Objectives' = " << _nNumObjs
                 << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            _cInitialF.resize (0);
            return( false );
        }
        if (bModifiedInitialX == true)
        {
            cerr << "WARNING: Ignoring 'Initial F' because 'Initial X'"
                 << " was modified" << endl;
            _cInitialF.resize (0);
        }
    }
    else if (cParams.isParameterDouble ("Initial F"))
    {
        cerr << "WARNING: Parameter 'Initial F' in '" << sPROBDEF << "' sublist"
             << " should be a vector" << endl;
        cerr << "         Ignoring 'Initial F'" << endl;
    }

    if (cParams.isParameterVector ("Initial Nonlinear Eqs"))
    {
        _cInitialEqs = cParams.getVectorParameter ("Initial Nonlinear Eqs");
        if (_cInitialEqs.size() != _nNumNonlEqs)
        {
            cerr << "ERROR: Length of 'Initial Nonlinear Eqs' = "
                 << _cInitialEqs.size()
                 << " does not match 'Number Nonlinear Eqs' = "
                 << _nNumNonlEqs << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            _cInitialEqs.resize (0);
            return( false );
        }
        if (bModifiedInitialX == true)
        {
            cerr << "WARNING: Ignoring 'Initial Nonlinear Eqs'"
                 << " because 'Initial X' was modified" << endl;
            _cInitialEqs.resize (0);
        }
    }

    if (cParams.isParameterVector ("Initial Nonlinear Ineqs"))
    {
        _cInitialIneqs = cParams.getVectorParameter ("Initial Nonlinear Ineqs");
        if (_cInitialIneqs.size() != _nNumNonlIneqs)
        {
            cerr << "ERROR: Length of 'Initial Nonlinear Ineqs' = "
                 << _cInitialIneqs.size()
                 << " does not match 'Number Nonlinear Ineqs' = "
                 << _nNumNonlIneqs << endl;
            cerr << "       See sublist '" << sPROBDEF << "'" << endl;
            _cInitialIneqs.resize (0);
            return( false );
        }
        if (bModifiedInitialX == true)
        {
            cerr << "WARNING: Ignoring 'Initial Nonlinear Ineqs'"
                 << " because 'Initial X' was modified" << endl;
            _cInitialIneqs.resize (0);
        }
    }

    //---- USER CANNOT HAVE EVALUATION VALUES WITH NO INITIAL POINT.
    //---- (BUT IT IS OK TO HAVE AN INITIAL POINT WITH NO EVALUATIONS.)
    if (_cInitialX.empty() && (_cInitialF.empty() == false))
    {
        cerr << "WARNING: Ignoring 'Initial F' in '" << sPROBDEF << "' sublist"
             << "; need 'Initial X'" << endl;
        _cInitialF.resize (0);
    }
    if (_cInitialX.empty() && (_cInitialEqs.empty() == false))
    {
        cerr << "WARNING: Ignoring 'Initial Nonlinear Eqs' in '"
             << sPROBDEF << "' sublist" << "; need 'Initial X'" << endl;
        _cInitialEqs.resize (0);
    }
    if (_cInitialX.empty() && (_cInitialIneqs.empty() == false))
    {
        cerr << "WARNING: Ignoring 'Initial Nonlinear Ineqs' in '"
             << sPROBDEF << "' sublist" << "; need 'Initial X'" << endl;
        _cInitialEqs.resize (0);
    }

    //---- CHECK THAT ALL INITIAL EVALUATION DATA IS PROVIDED, NOT JUST
    //---- PART OF IT.
    if (   (_cInitialF.empty() == false)
        || (_cInitialEqs.empty() == false)
        || (_cInitialIneqs.empty() == false) )
    {
        bool  bIsOK = true;
        if (_cInitialF.size() != _nNumObjs)
        {
            cerr << "WARNING: Ignoring initial point data in '"
                 << sPROBDEF << "' sublist; need 'Initial F'" << endl;
            bIsOK = false;
        }
        if (_cInitialEqs.size() != _nNumNonlEqs)
        {
            cerr << "WARNING: Ignoring initial point data in '"
                 << sPROBDEF << "' sublist; need 'Initial Nonlinear Eqs'"
                 << endl;
            bIsOK = false;
        }
        if (_cInitialIneqs.size() != _nNumNonlIneqs)
        {
            cerr << "WARNING: Ignoring initial point data in '"
                 << sPROBDEF << "' sublist; need 'Initial Nonlinear Ineqs'"
                 << endl;
            bIsOK = false;
        }
        if (bIsOK == false)
        {
            _cInitialF.resize (0);
            _cInitialEqs.resize (0);
            _cInitialIneqs.resize (0);
        }
    }

    return( true );
}


//----------------------------------------------------------------------
//  Private Method setupMisc_
//----------------------------------------------------------------------
bool  ProblemDef::setupMisc_ (const ParameterList &  cParams)
{
    _nDisplayFlag = cParams.getParameter ("Display", _nDisplayFlag);
    if (_nDisplayFlag < 0)
        _nDisplayFlag = 0;
    if (_nDisplayFlag > 2)
        _nDisplayFlag = 2;

    _nNumNonlEqs = cParams.getParameter ("Number Nonlinear Eqs", 0);
    if (_nNumNonlEqs < 0)
    {
        cerr << "WARNING: Cannot have negative 'Number Nonlinear Eqs'"
             << " in '" << sPROBDEF << "' sublist" << endl;
        cerr << "         Changing 'Number Nonlinear Eqs' to zero" << endl;
        _nNumNonlEqs = 0;
    }

    _nNumNonlIneqs = cParams.getParameter ("Number Nonlinear Ineqs", 0);
    if (_nNumNonlIneqs < 0)
    {
        cerr << "WARNING: Cannot have negative 'Number Nonlinear Ineqs'"
             << " in '" << sPROBDEF << "' sublist" << endl;
        cerr << "         Changing 'Number Nonlinear Ineqs' to zero" << endl;
        _nNumNonlIneqs = 0;
    }

    _dNonlActTol
        = cParams.getParameter ("Nonlinear Active Tolerance", _dNonlActTol);

    return( true );
}


//----------------------------------------------------------------------
//  Method printDefinition
//----------------------------------------------------------------------
void  ProblemDef::printDefinition (const bool  bDisplayFull) const
{
    if (_nDisplayFlag <= 0)
        return;

    if ((_nDisplayFlag < 2) || (bDisplayFull == false))
    {
        //---- DISPLAY SUMMARY ONLY.

        cout << "Problem Definition" << endl;
        printObjDefinition_();
        printVarSummary_();
        printInitPointSummary_();
        cout << endl;
    }
    else
    {
        //---- DISPLAY ALL DETAILS.

        cout << "Problem Definition (full display)" << endl;
        printObjDefinition_();
        printVarSummary_();

        //---- PRINT A LIST WITH EVERY VARIABLE, ITS BOUNDS, AND TYPE.
        cout << "  Variable bounds and scaling:" << endl;
        for (int  i = 0; i < _nNumVars; i++)
        {
            cout << "  ";

            if (_naVarTypes[i] == CONTINUOUS)
                cout << " (cont)     ";
            else if (_naVarTypes[i] == INTEGER)
                cout << " (integer)  ";
            else if (_naVarTypes[i] == ORDINAL)
                cout << " (ordinal)  ";
            else
                cout << " (unknown)  ";

            if (exists (_cLoBnds[i]))
            {
                cout << setw(14) << setprecision(6)
                     << setiosflags (ios::scientific) << _cLoBnds[i];
                cout << " <= ";
            }
            else
                cout << "              " << "    ";

            printVarName_ (i);

            if (exists (_cUpBnds[i]))
            {
                cout << " <= ";
                cout << setw(14) << setprecision(6)
                     << setiosflags (ios::scientific) << _cUpBnds[i];
            }
            else
                cout << "    "  << "              ";

            cout << "  scale=" << setw(11) << setprecision(4)
                 << setiosflags (ios::scientific) << _cScaling[i];
            cout << endl;
        }

        printInitPointSummary_();
        if (_cInitialX.empty() == false)
        {
            for (int  i = 0; i < _nNumVars; i++)
            {
                cout << "    Initial ";
                printVarName_ (i);
                cout << " = " << setw(14) << setprecision(6)
                     << setiosflags (ios::scientific) << _cInitialX[i] << endl;
            }
            if (_cInitialEqs.empty() == false)
            {
                for (int  i = 0; i < _nNumNonlEqs; i++)
                {
                    cout << "    Initial c_e[" << setw(5) << i << "]";
                    cout << " = " << setw(14) << setprecision(6)
                         << setiosflags (ios::scientific) << _cInitialEqs[i]
                         << "     (nonlinear eq)" << endl;
                }
            }
            if (_cInitialIneqs.empty() == false)
            {
                for (int  i = 0; i < _nNumNonlIneqs; i++)
                {
                    cout << "    Initial c_i[" << setw(5) << i << "]";
                    cout << " = " << setw(14) << setprecision(6)
                         << setiosflags (ios::scientific) << _cInitialIneqs[i]
                         << "     (nonlinear ineq)" << endl;
                }
            }
        }

        cout << "End of Problem Definition (full display)" << endl;
        cout << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method printObjDefinition_
//----------------------------------------------------------------------
void ProblemDef::printObjDefinition_ (void) const
{
    if (_nNumObjs == 1)
    {
        if (_nObjType == MINIMIZE)
        {
            cout << "  Minimize 1 objective";
            if (_dObjTarget != dne())
                cout << ", objective target = " << _dObjTarget;
            cout << endl;
        }
        else if (_nObjType == MAXIMIZE)
        {
            cout << "  Maximize 1 objective";
            if (_dObjTarget != dne())
                cout << ", objective target = " << _dObjTarget;
            cout << endl;
        }
        else if (_nObjType == FIND_FEASIBLE_PT)
        {
            cout << "  Find any feasible point (no objective target)" << endl;
        }
    }
    else
    {
        cout << "  " << _nNumObjs << " objectives" << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method printVarSummary_
//----------------------------------------------------------------------
void ProblemDef::printVarSummary_ (void) const
{
    int  nNumCont    = 0;
    int  nNumInteger = 0;
    int  nNumOrdinal = 0;
    for (int  i = 0; i < _nNumVars; i++)
    {
        if (_naVarTypes[i] == CONTINUOUS)
            nNumCont++;
        else if (_naVarTypes[i] == INTEGER)
            nNumInteger++;
        else if (_naVarTypes[i] == ORDINAL)
            nNumOrdinal++;
    }
    cout << "  "   << setw (5) << _nNumVars << " variables" << endl;
    cout << "    " << setw (5) << nNumCont << " continuous variables" << endl;
    cout << "    " << setw (5) << nNumInteger << " integer variables" << endl;
    cout << "    " << setw (5) << nNumOrdinal << " ordinal variables" << endl;

    int  nNumLoBndOnly = 0;
    int  nNumUpBndOnly = 0;
    int  nNumBothBnds  = 0;
    int  nNumFree      = 0;
    for (int  i = 0; i < _nNumVars; i++)
    {
        if (exists (_cLoBnds[i]) && exists (_cUpBnds[i]))
            nNumBothBnds++;
        else if (exists (_cLoBnds[i]))
            nNumLoBndOnly++;
        else if (exists (_cUpBnds[i]))
            nNumUpBndOnly++;
        else
            nNumFree++;
    }
    cout << "  " << setw (5) << nNumBothBnds
         << " vars with bounds above and below" << endl;
    cout << "  " << setw (5) << nNumUpBndOnly
         << " vars with upper bound only" << endl;
    cout << "  " << setw (5) << nNumLoBndOnly
         << " vars with lower bound only" << endl;
    cout << "  " << setw (5) << nNumFree
         << " vars with no upper or lower bound" << endl;

    if ((_nNumNonlEqs > 0) || (_nNumNonlIneqs > 0))
    {
        cout << "  " << setw (5) << _nNumNonlEqs
             << " nonlinear equality constraints" << endl;
        cout << "  " << setw (5) << _nNumNonlIneqs
             << " nonlinear inequality constraints" << endl;
        cout << "  Tolerance for nonlinear constraint feasibility = "
             << setw(14) << setprecision(6)
             << setiosflags (ios::scientific) << _dNonlActTol << endl;
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method printVarName_
//----------------------------------------------------------------------
void ProblemDef::printVarName_ (const int  nVarNum) const
{
    cout << "x[" << setw (5) << nVarNum << "]";
    return;
}


//----------------------------------------------------------------------
//  Private Method printInitPointSummary_
//----------------------------------------------------------------------
void ProblemDef::printInitPointSummary_ (void) const
{
    if (_cInitialX.empty())
        cout << "  Initial point not defined" << endl;
    else
    {
        cout << "  Initial point defined";
        if (_cInitialF.empty() == false)
        {
            if (_nNumObjs == 1)
            {
                cout << ", with objective value = ";
                if (exists (_cInitialF[0]))
                    cout << setw(19) << setprecision(11)
                         << setiosflags (ios::scientific) << _cInitialF[0];
                else
                    cout << "DNE";
            }
            else if (_nNumObjs > 1)
            {
                cout << ", with objective values = [ ";
                for (int  i = 0; i < _nNumObjs; i++)
                {
                    if (exists (_cInitialF[i]))
                        cout << setw(19) << setprecision(11)
                             << setiosflags (ios::scientific) << _cInitialF[i]
                             << " ";
                    else
                        cout << "DNE ";
                }
                cout << "]";
            }
        }
        else
            cout << ", but no objective value";
        cout << endl;
    }

    return;
}


}     //-- namespace HOPSPACK
