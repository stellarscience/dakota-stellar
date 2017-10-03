// $Id: HOPSPACK_ConveyorList.cpp 217 2013-11-25 21:59:49Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_ConveyorList.cpp $ 

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
  \file HOPSPACK_ConveyorList.cpp
  \brief Implement HOPSPACK::ConveyorList
*/

#include <list>

#include "HOPSPACK_CacheManager.hpp"
#include "HOPSPACK_ConveyorList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_ScaledComparison.hpp"


HOPSPACK::ConveyorList::ConveyorList (void)
{
}

HOPSPACK::ConveyorList::~ConveyorList()
{
  prune();
}

void HOPSPACK::ConveyorList::push(DataPoint* pushed)
{
  _cPtList.push_front(pushed);
}

HOPSPACK::DataPoint* HOPSPACK::ConveyorList::pop()
{
  if (_cPtList.empty())
    return NULL;

  DataPoint* popped = _cPtList.back();
  _cPtList.pop_back();
  return popped;
}

HOPSPACK::DataPoint* HOPSPACK::ConveyorList::pop(int tag)
{
  if (_cPtList.empty())
    return NULL;

  PointListType::iterator it = _cPtList.begin();
  for (; it != _cPtList.end(); it ++)
    if (((*it)->getTag()) == tag)
      break;

  if (it == _cPtList.end())
    return NULL;

  DataPoint* popped = *it;
  _cPtList.erase(it);

  return popped;
}

HOPSPACK::DataPoint *  HOPSPACK::ConveyorList::popNextCached
    (CacheManager * const  pCache)
{
    if (_cPtList.empty())
        return( NULL );

    PointListType::reverse_iterator  it;
    for (it = _cPtList.rbegin(); it != _cPtList.rend(); it++)
    {
        Vector  cDummyF;
        Vector  cDummyE;
        Vector  cDummyI;
        if (pCache->isCached ((*it)->getX(), cDummyF, cDummyE, cDummyI))
        {
            DataPoint *  pResult = *it;

            //---- TO ERASE, THE REVERSE ITERATOR HAS TO BE CONVERTED
            //---- TO A FORWARD ITERATOR USING base(), AND IT MUST BE MOVED.
            //---- SEE http://www.ddj.com/cpp/184401406.
            _cPtList.erase ((++it).base());

            return( pResult );
        }
    }

    return( NULL );
}

int HOPSPACK::ConveyorList::size() const
{
    return( (int) _cPtList.size() );
}

bool HOPSPACK::ConveyorList::isEmpty() const
{
  return _cPtList.empty();
}

const list< HOPSPACK::DataPoint * > &  HOPSPACK::ConveyorList::getPointList() const
{
    return( (const list< DataPoint *> &) _cPtList );
}

list< HOPSPACK::DataPoint * > &  HOPSPACK::ConveyorList::getMutablePointList()
{
    return( _cPtList );
}

list< HOPSPACK::DataPoint * >::const_iterator HOPSPACK::ConveyorList::begin() const
{
  return _cPtList.begin();
}

list< HOPSPACK::DataPoint * >::const_iterator HOPSPACK::ConveyorList::end() const
{
  return _cPtList.end();
}

void HOPSPACK::ConveyorList::setPriority (const int  nPriority)
{
  _nPriority = nPriority;

  if (_nPriority < 1)
    _nPriority = 1;
  if (_nPriority > 10)
    _nPriority = 10;

  return;
}

int HOPSPACK::ConveyorList::getPriority (void) const
{
  return( _nPriority );
}

void HOPSPACK::ConveyorList::getTagList(vector<int>& t) const
{
  t.clear();
  PointListType::const_iterator it;
  for (it = _cPtList.begin(); it != _cPtList.end(); it++)
  {
    int tag = (*it)->getTag();
    t.insert(t.begin(), tag);
  }
}

void HOPSPACK::ConveyorList::prune(int n)
{
  if (n <= 0)
  {
    PointListType::iterator it;
    for (it = _cPtList.begin(); it != _cPtList.end(); it++)
      delete *it;
    
    _cPtList.clear();
  }
  else
  {
    int p = size() - n;
    DataPoint* popped;
    for (int i = 0; i < p; i ++)
    {
      popped = pop();
      delete popped;
    }
  }
}

bool HOPSPACK::ConveyorList::contains(const DataPoint* w, int& tag) const
{
  tag = -1;
  if (_cPtList.empty())
    return false;
  
  const Vector& wx = w->getX();
  PointListType::const_iterator tpi;
  for (tpi = _cPtList.begin(); tpi != _cPtList.end(); tpi ++)
  {
    const Vector& zx = (*tpi)->getX();
    if (ScaledComparison::isEqual (wx, zx) == true)
    {
      tag = (*tpi)->getTag();
      return true;
    }
  }
  
  return false;
}

void HOPSPACK::ConveyorList::print(const string label) const
{
  cout << label << ":" << endl;
  
  if (_cPtList.empty())
  {
    cout << "  <empty>" << endl;
    return;
  }

  PointListType::const_iterator it;
  for (it = _cPtList.begin(); it != _cPtList.end(); it++)
  {
    DataPoint *  pTmp = *it;
    pTmp->leftshift (cout, true);
    cout << endl;
  }

  return;
}
