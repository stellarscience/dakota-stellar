// $Id: HOPSPACK_DataPoint.cpp 202 2012-03-07 22:26:01Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_DataPoint.cpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2012 Sandia Corporation
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
  @file HOPSPACK_DataPoint.cpp
  @brief Implement HOPSPACK::DataPoint.
*/

#include <math.h>                 //-- FOR fabs, sqrt, log10
#include <algorithm>              //-- FOR set_difference
#include <iterator>               //-- FOR insert_iterator

#include "HOPSPACK_DataPoint.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_Print.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Define static variables
//----------------------------------------------------------------------

int  DataPoint::_nTagGlobalCounter = 0;


//---- CHANGE THIS TO ENABLE MEMORY LEAK INFORMATION.
//---- DO NOT SET IT true FOR PRODUCTION CODE AS LISTS CAN GET QUITE LARGE.
bool  DataPoint::_bDebuggingLeaks = false;

int  DataPoint::_nDebugTagCounter = 0;
vector< int >  DataPoint::_cDebugCreateList;
vector< int >  DataPoint::_cDebugDeleteList;


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
DataPoint::DataPoint (const ProblemDef::ObjectiveType  nObjGoal,
                      const Vector &                   cX)
    :
    _nObjGoal (nObjGoal),
    _cX (cX),
    _nTag (_nTagGlobalCounter),
    _nState (UNEVALUATED)
{
    _nTagGlobalCounter++;

    if (_bDebuggingLeaks)
    {
        _nDebugTag = (_nDebugTagCounter++);
        cout << "=====Debug DataPoint constructor, tag = "
             << _nDebugTag << endl;
        _cDebugCreateList.push_back (_nDebugTag);
    }

    return;
}


//----------------------------------------------------------------------
//  Copy Constructor
//----------------------------------------------------------------------
DataPoint::DataPoint (const DataPoint &  cArg)
    :
    _nObjGoal (cArg._nObjGoal),
    _cX (cArg._cX),
    _cFns (cArg._cFns),
    _cEqs (cArg._cEqs),
    _cIneqs (cArg._cIneqs),
    _nTag (cArg._nTag),
    _sMsg (cArg._sMsg),
    _nState (cArg._nState)
{
    if (_bDebuggingLeaks)
    {
        _nDebugTag = (_nDebugTagCounter++);
        cout << "=====Debug DataPoint copy constructor, tag = "
             << _nDebugTag << endl;
        _cDebugCreateList.push_back (_nDebugTag);
    }

    return;
}


//----------------------------------------------------------------------
//  Assignment Operator
//----------------------------------------------------------------------
DataPoint &  DataPoint::operator= (const DataPoint &  cArg)
{
    _nObjGoal = cArg._nObjGoal;
    _cX = cArg._cX;
    _cFns = cArg._cFns;
    _cEqs = cArg._cEqs;
    _cIneqs = cArg._cIneqs;
    _nTag = cArg._nTag;
    _sMsg = cArg._sMsg;
    _nState = cArg._nState;

    if (_bDebuggingLeaks)
    {
        _nDebugTag = (_nDebugTagCounter++);
        cout << "=====Debug DataPoint operator=, tag = "
             << _nDebugTag << endl;
        _cDebugCreateList.push_back (_nDebugTag);
    }

    return( *this );
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
DataPoint::~DataPoint (void)
{
    if (_bDebuggingLeaks)
    {
        cout << "=====Debug DataPoint destructor, tag = "
             << _nDebugTag << endl;
        _cDebugDeleteList.push_back (_nDebugTag);
    }

    return;
}


//----------------------------------------------------------------------
//  Accessors
//----------------------------------------------------------------------

const Vector &  DataPoint::getX (void) const
{
    return( _cX );
}

const Vector &  DataPoint::getVecF (void) const
{
    return( _cFns );
}

double  DataPoint::getBestF (void) const
{
    if (_nState == UNEVALUATED)
        return( dne() );
    if (_nObjGoal == ProblemDef::FIND_FEASIBLE_PT)
        return( dne() );
    if (_cFns.empty())
        return( dne() );

    if (_cFns.size() == 1)
        return( _cFns[0] );

    double  dResult = dne();
    for (int  i = 0; i < _cFns.size(); i++)
    {
        if (exists (_cFns[i]))
        {
            if (exists (dResult) == false)
                dResult = _cFns[i];
            else
            {
                if ((_nObjGoal == ProblemDef::MAXIMIZE) && (_cFns[i] > dResult))
                    dResult = _cFns[i];
                if ((_nObjGoal == ProblemDef::MINIMIZE) && (_cFns[i] < dResult))
                    dResult = _cFns[i];
            }
        }
    }
    return( dResult );
}

const Vector &  DataPoint::getEqs (void) const
{
    return( _cEqs );
}

const Vector &  DataPoint::getIneqs (void) const
{
    return( _cIneqs );
}

double  DataPoint::getNonlConstrL2Norm (void) const
{
    double  dEq = _cEqs.norm();

    //---- INEQUALITY VALUES ARE NEGATIVE IF IN VIOLATION.
    double  dIneqSqrd = 0.0;
    for (int  i = 0; i < (int) _cIneqs.size(); i++)
        if (_cIneqs[i] < 0.0)
            dIneqSqrd += _cIneqs[i] * _cIneqs[i];

    double  dResultSqrd = (dEq * dEq) + dIneqSqrd;
    return( sqrt (dResultSqrd) );
}

double  DataPoint::getNonlConstrLInfNorm (void) const
{
    double  dResult = 0.0;
    for (int  i = 0; i < _cEqs.size(); i++)
    {
        double  dTmp = fabs (_cEqs[i]);
        if (dTmp > dResult)
            dResult = dTmp;
    }
    //---- INEQUALITY VALUES ARE NEGATIVE IF IN VIOLATION.
    for (int  i = 0; i < _cIneqs.size(); i++)
    {
        if ((- _cIneqs[i]) > dResult)
            dResult = - _cIneqs[i];
    }

    return( dResult );
}

double  DataPoint::getPenaltySign (void) const
{
    if (_nObjGoal == ProblemDef::MAXIMIZE)
        return( -1.0 );
    return( 1.0 );
}

int  DataPoint::getTag (void) const
{
    return( _nTag );
}

DataPoint::State  DataPoint::getState (void) const
{
    return( _nState );
}


//----------------------------------------------------------------------
//  Manipulators
//----------------------------------------------------------------------

void  DataPoint::setEvalFC (const Vector &  cFns,
                            const Vector &  cEqs,
                            const Vector &  cIneqs,
                            const string &  sMsg)
{
    _cFns = cFns;
    _cEqs = cEqs;
    _cIneqs = cIneqs;
    _sMsg = sMsg;
    _nState = EVALUATED_FC;
    return;
}

void  DataPoint::setCachedFC (const Vector &  cFns,
                              const Vector &  cEqs,
                              const Vector &  cIneqs,
                              const string &  sMsg)
{
    _cFns = cFns;
    _cEqs = cEqs;
    _cIneqs = cIneqs;
    _sMsg = sMsg;
    _nState = EVALUATED_FC_FROM_CACHE;
    return;
}


//----------------------------------------------------------------------
//  Method isBetterObjThan (another point)
//----------------------------------------------------------------------
bool  DataPoint::isBetterObjThan (const DataPoint &  cOther,
                                        bool      &  bAreObjsComparable) const
{
    bAreObjsComparable = false;

    //---- ARBITRARILY TRUE IF THE OBJECTIVE DOES NOT MATTER.
    if (_nObjGoal == ProblemDef::FIND_FEASIBLE_PT)
        return( true );

    bool  bIsThisEval = (_nState == EVALUATED_FC)
                        || (_nState == EVALUATED_FC_FROM_CACHE);

    bool  bIsOtherEval = (cOther._nState == EVALUATED_FC)
                         || (cOther._nState == EVALUATED_FC_FROM_CACHE);

    //---- AN EVALUATED POINT IS BETTER THAN AN UNEVALUATED POINT.
    if (bIsThisEval && (bIsOtherEval == false))
        return( true );
    if (bIsOtherEval && (bIsThisEval == false))
        return( false );

    //---- IF NEITHER POINT HAS BEEN EVALUATED, USE THE TAG TIE-BREAKER.
    if ((bIsThisEval == false) && (bIsOtherEval == false))
        return( this->getTag() < cOther.getTag() );

    //---- FALSE IF POINTS HAVE BEEN EVALUATED AND ARE IDENTICAL.
    if (this->getTag() == cOther.getTag())
        return( false );

    double  dThisF  = this->getBestF();
    double  dOtherF = cOther.getBestF();
    
    //---- A POINT WHOSE OBJECTIVE EXISTS IS BETTER THAN A POINT WITH NONE.
    if (exists (dThisF) && (exists (dOtherF) == false))
        return( true );
    if (exists (dOtherF) && (exists (dThisF) == false))
        return( false );

    //---- IF NEITHER OBJECTIVE EXISTS, USE THE TAG TIE-BREAKER.
    if ((exists (dThisF) == false) && (exists (dOtherF) == false))
        return( this->getTag() < cOther.getTag() );

    bAreObjsComparable = true;

    //---- COMPARE THE OBJECTIVES.
    if (dThisF < dOtherF)
    {
        if (_nObjGoal == ProblemDef::MINIMIZE)
            return( true );
        else
            return( false );
    }
    if (dThisF > dOtherF)
    {
        if (_nObjGoal == ProblemDef::MINIMIZE)
            return( false );
        else
            return( true );
    }

    //---- TIE-BREAKER IF OBJECTIVES ARE EQUAL.
    return( this->getTag() < cOther.getTag() );
}


//----------------------------------------------------------------------
//  Method isSamePoint
//----------------------------------------------------------------------
bool  DataPoint::isSamePoint (const DataPoint &  cOther,
                              const double       dTolerance) const
{
    if (_cX.size() != cOther._cX.size())
    {
        cerr << "ERROR: Bad argument length"
             << "  <DataPoint::isSamePoint()>" << endl;
        throw INTERNAL_ERROR;
    }

    for (int  i = 0; i < (int) _cX.size(); i++)
    {
        if (fabs (_cX[i] - cOther._cX[i]) > dTolerance)
            return( false );
    }
    return( true );
}


//----------------------------------------------------------------------
//  Method leftshift
//----------------------------------------------------------------------
void  DataPoint::leftshift (      ostream &  stream,
                            const bool       bIncludeMsg,
                            const bool       bPrintAllX) const
{
    stream << "  ";
    stream << "Tag=" << _nTag;
    if (_bDebuggingLeaks)
        stream << " (dbgTag=" << _nDebugTag << ")";

    if ((bPrintAllX == false) && (_cX.size() >= 10))
    {
        stream << ", Size of x=" << _cX.size();
    }
    else
    {
        stream << ", x=[";
        _cX.leftshift (stream);
        stream << "]";
    }

    stream << ", State=";
    if (_nState == UNEVALUATED)
        stream << "UNEVL";
    else if (_nState == EVALUATED_FC)
        stream << "EVL-F";
    else if (_nState == EVALUATED_FC_FROM_CACHE)
        stream << "CACHE";
    else
        stream << "???";

    if (_nState != UNEVALUATED)
    {
        if (bIncludeMsg)
            stream << ", '" << _sMsg << "'";

        if (_cFns.size() > 1)
            stream << ", F=[";
        else
            stream << ", F=";
        _cFns.leftshift (stream);
        if (_cFns.size() > 1)
            stream << "]";

        if ((_cEqs.size() > 0) || (_cIneqs.size() > 0))
        {
            stream << endl;
            stream << "      ";
            int  kDigits = 1 + ((int) log10 ((double) _nTag));
            for (int  i = 0; i < kDigits; i++)
                stream << " ";

            stream << ", c_e=[";
            _cEqs.leftshift (stream);
            stream << "]";

            stream << ", c_i=[";
            _cIneqs.leftshift (stream);
            stream << "]";
        }
    }

    return;
}


//----------------------------------------------------------------------
//  Method debugPrintMemoryLists
//----------------------------------------------------------------------
void  DataPoint::debugPrintMemoryLists (void)
{
    if (_bDebuggingLeaks == false)
        return;

    cout << "=====Debug DataPoint memory lists (begin)" << endl;
    cout << "  Total  instances = " << _nDebugTagCounter << endl;
    cout << "  Unique instances = " << _nTagGlobalCounter << endl;

    for (int  i = 0; i < (int) _cDebugCreateList.size(); i++)
    {
        cout << "  Created " << _cDebugCreateList[i];
        if (i < (int) _cDebugDeleteList.size())
            cout << "  Deleted " << _cDebugDeleteList[i] << endl;
    }

    sort (_cDebugCreateList.begin(), _cDebugCreateList.end());
    sort (_cDebugDeleteList.begin(), _cDebugDeleteList.end());
    vector< int > cDifferences;
    set_difference (_cDebugCreateList.begin(), _cDebugCreateList.end(),
                    _cDebugDeleteList.begin(), _cDebugDeleteList.end(),
                    insert_iterator< vector< int > > (cDifferences,
                                                      cDifferences.begin()) );
    cout << endl;
    cout << "  Tag(s) of points not deleted:";
    vector< int >::iterator  it;
    for (it = cDifferences.begin(); it != cDifferences.end(); it++)
        cout << " " << (*it);
    cout << endl;

    cout << "=====Debug DataPoint memory lists (end)" << endl;
    return;
}


}     //-- namespace HOPSPACK
