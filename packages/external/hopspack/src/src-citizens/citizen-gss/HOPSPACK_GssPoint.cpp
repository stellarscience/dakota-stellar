// $Id: HOPSPACK_GssPoint.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssPoint.cpp $

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
  @file HOPSPACK_GssPoint.cpp
  @brief Implement HOPSPACK::GssPoint, subclass of DataPoint.
*/

#include <iomanip>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_GssPoint.hpp"
#include "HOPSPACK_Print.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructors
//----------------------------------------------------------------------

GssPoint::GssPoint (const ProblemDef::ObjectiveType  nObjGoal,
                    const NonlConstrPenalty &        cPenalty,
                    const Vector &                   cX,
                    const double                     dStep)
    : DataPoint (nObjGoal, cX),
      _nParentTag (NO_PARENT_TAG),
      _nDirIndex (NO_DIR_INDEX),
      _dStep (dStep),
      _dParentObjective (dne()),
      _dSufficientImprovementAmount (dne()),
      _cPenalty (cPenalty)
{
    return;
}

GssPoint::GssPoint (const ProblemDef::ObjectiveType  nObjGoal,
                    const NonlConstrPenalty &        cPenalty,
                    const Vector &                   cX,
                    const double                     dStep,
                    const int                        nParentTag,
                    const double                     dParentObj,
                    const double                     dPenaltyTerm,
                    const double                     dSuffImprvAmt,
                    const int                        nDirIndex)
    : DataPoint (nObjGoal, cX),
      _nParentTag (nParentTag),
      _nDirIndex (nDirIndex),
      _dStep (dStep),
      _dSufficientImprovementAmount (dSuffImprvAmt),
      _cPenalty (cPenalty)
{
    if (_nObjGoal == ProblemDef::MAXIMIZE)
        _dParentObjective = dParentObj - dPenaltyTerm;
    else     //-- MINIMIZE OR FIND_FEASIBLE_PT
        _dParentObjective = dParentObj + dPenaltyTerm;
    return;
}

GssPoint::GssPoint (const DataPoint         &  cArg,
                    const NonlConstrPenalty &  cPenalty,
                    const double               dStep)
    : DataPoint (cArg),
      _nParentTag (NO_PARENT_TAG),
      _nDirIndex (NO_DIR_INDEX),
      _dStep (dStep),
      _dParentObjective (dne()),
      _dSufficientImprovementAmount (dne()),
      _cPenalty (cPenalty)
{
    return;
}



//----------------------------------------------------------------------
//  Copy Constructor
//----------------------------------------------------------------------
GssPoint::GssPoint (const GssPoint &  cArg)
    : DataPoint (cArg),
      _nParentTag (cArg._nParentTag),
      _nDirIndex (cArg._nDirIndex),
      _dStep (cArg._dStep),
      _dParentObjective (cArg._dParentObjective),
      _dSufficientImprovementAmount (cArg._dSufficientImprovementAmount),
      _cPenalty (cArg._cPenalty)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
GssPoint::~GssPoint (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method getParentTag
//----------------------------------------------------------------------
int  GssPoint::getParentTag (void) const
{
    return( _nParentTag );
}


//----------------------------------------------------------------------
//  Method getStepLength
//----------------------------------------------------------------------
double  GssPoint::getStepLength (void) const
{
    return( _dStep );
}


//----------------------------------------------------------------------
//  Method getDirIndex
//----------------------------------------------------------------------
int  GssPoint::getDirIndex (void) const
{
    return( _nDirIndex );
}


//----------------------------------------------------------------------
//  Method getBestF
//----------------------------------------------------------------------
double  GssPoint::getBestF (void) const
{
    double  dF = DataPoint::getBestF();

    if (_cPenalty.isDefined() == false)
        return( dF );

    double  dPen = _cPenalty.computePenalty (DataPoint::getEqs(),
                                             DataPoint::getIneqs());
    return( dF + (getPenaltySign() * dPen) );
}


//----------------------------------------------------------------------
//  Method isBetterObjThan
//----------------------------------------------------------------------
bool  GssPoint::isBetterObjThan (const GssPoint &  other) const
{
    bool  bAreObjsComparable;

    if (_cPenalty.isDefined() == false)
        return( DataPoint::isBetterObjThan (other, bAreObjsComparable) );


    //---- INEQUALITIES ASSUME THAT PENALTIES ARE ALWAYS NONNEGATIVE.

    const DataPoint &  cOtherDP = (const DataPoint &) other;
    double  dMyPen = _cPenalty.computePenalty (DataPoint::getEqs(),
                                               DataPoint::getIneqs());
    double  dOtherPen = _cPenalty.computePenalty (cOtherDP.getEqs(),
                                                  cOtherDP.getIneqs());

    bool  bIsFBetterThanOther
        = DataPoint::isBetterObjThan (other, bAreObjsComparable);
    if (bAreObjsComparable == false)
    {
        if (_nObjGoal == ProblemDef::FIND_FEASIBLE_PT)
            return( dMyPen < dOtherPen );
        else
            return( bIsFBetterThanOther );
    }

    //---- BOTH OBJECTIVES EXIST.
    double  dMyF = DataPoint::getBestF();
    double  dOtherF = cOtherDP.getBestF();
    if (_nObjGoal == ProblemDef::MINIMIZE)
        return( (dMyF + dMyPen) < (dOtherF + dOtherPen) );
    else
        return( (dMyF - dMyPen) > (dOtherF - dOtherPen) );
}


//----------------------------------------------------------------------
//  Method hasSufficientImprovement
//----------------------------------------------------------------------
bool  GssPoint::hasSufficientImprovement (void) const
{
    //---- SUFFICIENT IMPROVEMENT ONLY APPLIES TO POINTS GENERATED BY GSS.
    if (getParentTag() == NO_PARENT_TAG)
        return( true );

    double  dMyF = DataPoint::getBestF();

    if (_cPenalty.isDefined() == false)
    {
        //---- ARBITRARILY TRUE IF THE OBJECTIVE DOES NOT MATTER.
        if (_nObjGoal == ProblemDef::FIND_FEASIBLE_PT)
            return( true );

        if (exists (_dParentObjective) == false)
            return( true );
        if (exists (dMyF) == false)
            return( false );

        if (_nObjGoal == ProblemDef::MINIMIZE)
            return( dMyF < (_dParentObjective - _dSufficientImprovementAmount) );
        else
            return( dMyF > (_dParentObjective + _dSufficientImprovementAmount) );
    }

    //---- THE PARENT OBJECTIVE WILL INCLUDE THE PENALTY TERM.
    double  dMyPen = _cPenalty.computePenalty (DataPoint::getEqs(),
                                               DataPoint::getIneqs());
    if ( (exists (dMyF) == false) || (exists (dMyPen) == false) )
        return( false );
    if (exists (_dParentObjective) == false)
        return( true );
    if (_nObjGoal == ProblemDef::MAXIMIZE)
    {
        return( (dMyF - dMyPen)
                    > (_dParentObjective + _dSufficientImprovementAmount) );
    }
    else     //-- MINIMIZE OR FIND_FEASIBLE_PT
    {
        return( (dMyF + dMyPen)
                    < (_dParentObjective - _dSufficientImprovementAmount) );
    }
}


//----------------------------------------------------------------------
//  Method print
//----------------------------------------------------------------------
void  GssPoint::print (ostream &   stream,
                       const bool  bIncludeMsg) const
{
    DataPoint::leftshift (stream, bIncludeMsg, false);

    if (_cPenalty.isDefined())
    {
        cout.setf (ios::scientific);
        cout << ", p|C|="
             << setprecision (Print::getPrecision())
             << _cPenalty.computePenalty (DataPoint::getEqs(),
                                          DataPoint::getIneqs());
        cout.unsetf (ios::scientific);
    }

    cout << ", Step=" << _dStep;

    int  i = getParentTag();
    if (i == NO_PARENT_TAG)
        cout << ", ParentTag=(none)" << endl;
    else
        cout << ", ParentTag=" << getParentTag()
             << ", DirIx=" << _nDirIndex << endl;

    return;
}


}     //-- namespace HOPSPACK
