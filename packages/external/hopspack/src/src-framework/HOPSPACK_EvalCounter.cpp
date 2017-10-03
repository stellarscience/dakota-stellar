// $Id: HOPSPACK_EvalCounter.cpp 217 2013-11-25 21:59:49Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_EvalCounter.cpp $ 

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
  \file HOPSPACK_EvalCounter.cpp
  \brief Implement HOPSPACK::EvalCounter
*/

#include <iomanip>
#include <sstream>

#include "HOPSPACK_EvalCounter.hpp"


HOPSPACK::EvalCounter::EvalCounter() :
  nCached(0),
  nPendingCached(0),
  nEvaluated(0)
{
}

HOPSPACK::EvalCounter::~EvalCounter()
{
}

int HOPSPACK::EvalCounter::getNumEvaluations() const
{
  return nEvaluated;
}

int HOPSPACK::EvalCounter::getNumCached() const
{
  return nCached;
}

void  HOPSPACK::EvalCounter::getCountString (string &  sResult) const
{
    stringstream  sTmp;
    for (MsgCntConstIterator i = msgCnt.begin(); i != msgCnt.end(); i++)
        sTmp << " " << (*i).first << ": " << (*i).second;
    sResult = sTmp.str();
    return;
}

void HOPSPACK::EvalCounter::incrementCached()
{
  nCached++;
}

void HOPSPACK::EvalCounter::incrementPendingCached()
{
  nPendingCached++;
}

void HOPSPACK::EvalCounter::incrementEvaluated(int workerId, const string& msg)
{
  nEvaluated++;
  incrementMsgCnt(msgCnt, msg);
  incrementMsgCnt(wkrMsgCnt[workerId], msg);
}

// PRIVATE
void HOPSPACK::EvalCounter::incrementMsgCnt(MsgCnt& mc, const string& msg)
{
  MsgCntIterator i = mc.find(msg);
  if (i == mc.end())
    mc[msg] = 1;
  else
    (*i).second ++;
  return;
}

void HOPSPACK::EvalCounter::print (const bool  bDisplayDetails) const
{
  int  nWidth = 4;
  cout.setf (ios::fixed | ios::right);
  cout << "Evaluation breakdown:" << endl;
  cout << "  Number executed by workers = "
       << setw (nWidth) << nEvaluated << endl;
  cout << "  Number from cache          = "
       << setw (nWidth) << nCached << endl;
  cout << "  Number from pending list   = "
       << setw (nWidth) << nPendingCached << endl;

  if (bDisplayDetails == false)
    return;

  if (nEvaluated == 0)
    return;

  cout << "Evaluation executions grouped by message type:" << endl;
  int  nMaxWidth = 0;
  for (MsgCntConstIterator i = msgCnt.begin(); i != msgCnt.end(); i ++)
  {
      int  k = (int) ((*i).first).size();
      if (k > nMaxWidth)
          nMaxWidth = k;
  }
  nMaxWidth += 3;
  for (MsgCntConstIterator i = msgCnt.begin(); i != msgCnt.end(); i ++)
  {
    cout << "  '" << (*i).first << "'"
         << setw (nMaxWidth - ((*i).first).size()) << ": "
         << setw (3) << (*i).second << endl;
  }

  cout << "Evaluation executions grouped by worker ID then message type:"
       << endl;
  for (WkrMsgCntConstIterator j = wkrMsgCnt.begin(); j != wkrMsgCnt.end(); j++)
  {
    cout << "  Worker #" << (*j).first << endl;
    for (MsgCntConstIterator i = (*j).second.begin(); i != (*j).second.end(); i++)
    {
      cout << "    '" << (*i).first << "'"
           << setw (nMaxWidth - ((*i).first).size()) << ": "
           << setw (3) << (*i).second << endl;
    }
  }
  cout.unsetf(std::ios::fixed);
  cout.unsetf(std::ios::right);
  return;
}
