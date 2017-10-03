// $Id: HOPSPACK_ExecutorSerial.cpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-executor/HOPSPACK_ExecutorSerial.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2014 Sandia Corporation
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
  @file HOPSPACK_ExecutorSerial.cpp
  @brief Implement HOPSPACK::ExecutorSerial.
*/

#include <iomanip>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_ExecutorSerial.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_SystemTimer.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
ExecutorSerial::ExecutorSerial (Evaluator * const  pEvaluator)
    :
    _pEvaluator (pEvaluator)
{
    _pTimers = new SystemTimer (2);
    _pTimers->start (0);

    _bIsAvailable = true;
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
ExecutorSerial::~ExecutorSerial (void)
{
    if (_pTimers != NULL)
        delete _pTimers;
    return;
}


//----------------------------------------------------------------------
//  Method submit
//----------------------------------------------------------------------
bool  ExecutorSerial::submit (const int              nTag,
                              const Vector&          cX,
                              const EvalRequestType  nRequest)
{
    if (_bIsAvailable == false)
        return( false );

    if (Print::doPrint (Print::MOST_VERBOSE))
        cout << "ExecutorSerial calling Evaluator for tag " << nTag << endl;

    _nResultTag = nTag;

    //---- EVALUATE IMMEDIATELY, SAVING ALL RESULTS.
    _pTimers->start (1);
    if (nRequest == EVALREQTYPE_F)
    {
        _cResultF.resize (0);
        _pEvaluator->evalF (nTag, cX,
                            _cResultF, _sResultMsg);
    }
    else if (nRequest == EVALREQTYPE_FC)
    {
        _cResultF.resize (0);
        _cResultEqs.resize (0);
        _cResultIneqs.resize (0);
        _pEvaluator->evalFC (nTag, cX,
                             _cResultF, _cResultEqs, _cResultIneqs, _sResultMsg);
    }
    else
    {
        cerr << "ERROR: Evaluator request type " << nRequest
             << " not implemented <ExecutorSerial::submit>" << endl;
        throw INTERNAL_ERROR;
    }
    _pTimers->stop (1);

    _bIsAvailable = false;
    return( true );
}


//----------------------------------------------------------------------
//  Method isReadyForWork
//----------------------------------------------------------------------
bool  ExecutorSerial::isReadyForWork (void) const
{
    return( _bIsAvailable );
}


//----------------------------------------------------------------------
//  Method recv FC
//----------------------------------------------------------------------
int  ExecutorSerial::recv (int    &  nTag,
                           Vector &  cF,
                           Vector &  cEqs,
                           Vector &  cIneqs,
                           string &  sMsg)
{
    if (_bIsAvailable)
        //---- NO RESULT IS STORED, SO RETURN NOTHING.
        return( 0 );

    nTag = _nResultTag;
    cF = _cResultF;
    cEqs = _cResultEqs;
    cIneqs = _cResultIneqs;
    sMsg = _sResultMsg;

    _bIsAvailable = true;
    return( 1 );
}


//----------------------------------------------------------------------
//  Method getEvaluatorType
//----------------------------------------------------------------------
string  ExecutorSerial::getEvaluatorType (void) const
{
    return( _pEvaluator->getEvaluatorType() );
}


//----------------------------------------------------------------------
//  Method printDebugInfo
//----------------------------------------------------------------------
void  ExecutorSerial::printDebugInfo (void) const
{
    cout << "  HOPSPACK_ExecutorSerial -- does not use workers" << endl;
    cout << "    isReadyForWork() returns = " << isReadyForWork() << endl;
    _pEvaluator->printDebugInfo();

    return;
}


//----------------------------------------------------------------------
//  Method printTimingInfo
//----------------------------------------------------------------------
void  ExecutorSerial::printTimingInfo (void) const
{
    _pTimers->stop (0);
    cout << "Total wall clock time in serial Executor: "
         << _pTimers->getTotalTime (0) << " secs" << endl;

    int  nWidth = 8;
    streamsize  nCurrPrec  = cout.precision (3);
    cout.setf (ios::fixed | ios::right);
    cout << "  Serial evaluations      "
         << setw (nWidth) << _pTimers->getTotalTime (1) << endl;
    cout.precision (nCurrPrec);

    return;
}

    
}     //-- namespace HOPSPACK
