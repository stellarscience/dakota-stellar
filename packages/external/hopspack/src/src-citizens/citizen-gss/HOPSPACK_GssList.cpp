// $Id: HOPSPACK_GssList.cpp 217 2013-11-25 21:59:49Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssList.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2013 Sandia Corporation
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
  @file HOPSPACK_GssList.cpp
  @brief Implement HOPSPACK::GssList
*/

#include <algorithm>     //-- FOR find

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_GssList.hpp"
#include "HOPSPACK_NonlConstrPenalty.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
GssList::GssList (void)
    : _dDefaultStepLength (0.0)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
GssList::~GssList()
{
    prune();
    return;
}


//----------------------------------------------------------------------
//  Method setDefaultStepLength
//----------------------------------------------------------------------
void  GssList::setDefaultStepLength (const double  dStepLength)
{
    _dDefaultStepLength = dStepLength;
    return;
}


//----------------------------------------------------------------------
//  Method isEmpty
//----------------------------------------------------------------------
bool  GssList::isEmpty (void) const
{
    return( _cGssPtList.empty() );
}


//----------------------------------------------------------------------
//  Method push
//----------------------------------------------------------------------
void  GssList::push (GssPoint *  pNewPoint)
{
    _cGssPtList.push_front (pNewPoint);
    return;
}


//----------------------------------------------------------------------
//  Method findBest
//----------------------------------------------------------------------
const GssPoint *  GssList::findBest (void)
{
    moveBestToEndOfList_();
    return( _cGssPtList.back() );
}


//----------------------------------------------------------------------
//  Method pop
//----------------------------------------------------------------------
GssPoint *  GssList::pop (void)
{
    if (_cGssPtList.empty())
        return( NULL );

    GssPoint *  pResult = _cGssPtList.back();
    _cGssPtList.pop_back();
    return( pResult );
}


//----------------------------------------------------------------------
//  Method popBest
//----------------------------------------------------------------------
GssPoint *  GssList::popBest (void)
{
    moveBestToEndOfList_();
    return( pop() );
}


//----------------------------------------------------------------------
//  Method copyFrom
//----------------------------------------------------------------------
void  GssList::copyFrom (const list< DataPoint * > &  cInputList,
                         const NonlConstrPenalty   &  cPenalty,
                         const list< int >         &  cTagList)
{
    list< DataPoint * >::const_iterator  it;
    for (it = cInputList.begin(); it != cInputList.end(); it++)
    {
        int  nTag = (*it)->getTag();
        if (find (cTagList.begin(), cTagList.end(), nTag) != cTagList.end())
        {
            //---- THE POINT IS OURS; CAST IT AND MAKE A NEW COPY.
            const GssPoint *  pNext = (GssPoint *) (*it);
            _cGssPtList.push_back (new GssPoint (*pNext));
        }
        else
        {
            //---- THE POINT IS NOT OURS; COPY IT BUT SPECIFY THERE IS
            //---- NO PARENT AND USE THE DEFAULT STEP LENGTH.
            GssPoint *  pNext = new GssPoint (**it,
                                              cPenalty,
                                              _dDefaultStepLength);
            _cGssPtList.push_back (pNext);
        }
    }
    return;
}


//----------------------------------------------------------------------
//  Method copyTo
//----------------------------------------------------------------------
void  GssList::copyTo (list< DataPoint * > &  cOutputList)
{
    GssPointListType::reverse_iterator  it;
    for (it = _cGssPtList.rbegin(); it != _cGssPtList.rend(); it++)
    {
        cOutputList.push_front (new GssPoint (**it));
    }
    return;
}


//----------------------------------------------------------------------
//  Method insertFromList
//----------------------------------------------------------------------
void  GssList::insertFromList (const GssList &  cSource)
{
    _cGssPtList.insert (_cGssPtList.begin(),
                        cSource._cGssPtList.begin(),
                        cSource._cGssPtList.end());
    return;
}


//----------------------------------------------------------------------
//  Method copyTags
//----------------------------------------------------------------------
void  GssList::copyTags (vector< int > &  cTags) const
{
    GssPointListType::const_iterator  it;
    for (it = _cGssPtList.begin(); it != _cGssPtList.end(); it++)
        cTags.push_back ((*it)->getTag());
    return;
}


//----------------------------------------------------------------------
//  Method clearList
//----------------------------------------------------------------------
void  GssList::clearList (void)
{
    _cGssPtList.clear();
    return;
}


//----------------------------------------------------------------------
//  Method prune
//----------------------------------------------------------------------
void  GssList::prune (const int  n)
{
    if (n <= 0)
    {
        GssPointListType::iterator  it;
        for (it = _cGssPtList.begin(); it != _cGssPtList.end(); it ++)
            delete *it;
        _cGssPtList.clear();
    }
    else
    {
        int p = ((int) _cGssPtList.size()) - n;
        for (int i = 0; i < p; i++)
        {
            GssPoint *  popped = pop();
            delete popped;
        }
    }
    return;
}



//----------------------------------------------------------------------
//  Method print
//----------------------------------------------------------------------
void  GssList::print (const string label) const
{
    cout << label << ":" << endl;

    if (_cGssPtList.empty())
    {
        cout << "  <empty>" << endl;
        return;
    }

    GssPointListType::const_reverse_iterator  it;
    for (it = _cGssPtList.rbegin(); it != _cGssPtList.rend(); it++)
    {
        GssPoint *  pTmp = *it;
        pTmp->print (cout);
    }

    return;
}


//----------------------------------------------------------------------
//  Private Method moveBestToEndOfList_
//----------------------------------------------------------------------
void  GssList::moveBestToEndOfList_ (void)
{
    if (_cGssPtList.empty())
    {
        cerr << "ERROR: List is empty"
             << "       <GssList::moveBestToEndOfList()>." << endl;
        throw "GSS Error";
    }

    if (_cGssPtList.size() == 1)
        return;

    // Find the index of the best point
    GssPointListType::iterator  bestIterator = _cGssPtList.begin();
    GssPointListType::iterator  tpi = _cGssPtList.begin();

    for (tpi++; tpi != _cGssPtList.end(); tpi++)
        if ((*tpi)->isBetterObjThan (**bestIterator))
            bestIterator = tpi;

    // Swap the best point to the end of the list
    GssPoint *  tmp = *bestIterator;
    *bestIterator = _cGssPtList.back();
    _cGssPtList.back() = tmp;

    return;
}


}     //-- namespace HOPSPACK
