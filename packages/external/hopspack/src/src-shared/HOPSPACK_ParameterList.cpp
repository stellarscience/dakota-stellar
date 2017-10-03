// $Id: HOPSPACK_ParameterList.cpp 203 2012-05-14 22:27:30Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ParameterList.cpp $ 

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
  \file HOPSPACK_ParameterList.cpp
  \brief Implement HOPSPACK::ParameterList.
*/

#include "HOPSPACK_ParameterEntry.hpp"
#include "HOPSPACK_ParameterList.hpp"


//----------------------------------------------------------------------
// Static local data
//----------------------------------------------------------------------
static HOPSPACK::ParameterList  cEMPTY_SUBLIST;


//----------------------------------------------------------------------
// Constructors and Destructor
//----------------------------------------------------------------------

HOPSPACK::ParameterList::ParameterList()
{}

HOPSPACK::ParameterList::ParameterList(const ParameterList& source) 
{
  params = source.params;
}

HOPSPACK::ParameterList& HOPSPACK::ParameterList::operator=(const ParameterList& source)
{
  if (&source == this)
    return *this;

  params = source.params;
  return *this;
}

HOPSPACK::ParameterList::~ParameterList()
{
}


//----------------------------------------------------------------------
// Sublist methods
//----------------------------------------------------------------------

const HOPSPACK::ParameterList& HOPSPACK::ParameterList::sublist(const string& name) const
{
  // Find name in list, if it exists.
  ConstIterator i = params.find(name);

  // If it does not exist, return an empty list.
  if (i == params.end())
  {
      return( cEMPTY_SUBLIST );
  }

  // If it does exist and is a list, return the list value.
  if (entry(i).isList()) 
    return (entry(i).getListValue());

  // Otherwise, the parameter exists but is not a list. Throw an error.
  cerr << "ERROR: Parameter " << name << " is not a list." << endl;
  throw INTERNAL_ERROR;
}

HOPSPACK::ParameterList& HOPSPACK::ParameterList::getOrSetSublist(const string& name)
{
  // Find name in list, if it exists.
  Iterator i = params.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params.end()) {
    if (entry(i).isList()) 
      return (entry(i).getListValue());
    else
    {
      cerr << "ERROR: Parameter " << name << " is not a list." << endl;
      throw INTERNAL_ERROR;
    }
  }

  // If it does not exist, create a new empty list and return a reference
  return params[name].setList(true);
}


//----------------------------------------------------------------------
// Set parameter methods
//----------------------------------------------------------------------

void HOPSPACK::ParameterList::setParameter(const string& name, bool value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isBool() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, int value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isInt() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, double value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isDouble() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, const char* value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isString() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, const string& value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isString() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, const vector< char >  value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isVector() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, const Vector& value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isVector() == false))
        return;

    params[name].setValue(value);
    return;
}

void HOPSPACK::ParameterList::setParameter(const string& name, const Matrix& value)
{
    ConstIterator i = params.find(name);
    if ((i != params.end()) && (entry(i).isMatrix() == false))
        return;

    params[name].setValue(value);
    return;
}


//----------------------------------------------------------------------
// Get parameter methods
//----------------------------------------------------------------------
  
bool HOPSPACK::ParameterList::getParameter(const string& name, bool nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isBool()))
    return entry(i).getBoolValue();
  return nominal;
}

int HOPSPACK::ParameterList::getParameter(const string& name, int nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isInt()))
    return entry(i).getIntValue();
  return nominal;
}

double HOPSPACK::ParameterList::getParameter(const string& name, double nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isDouble()))
    return entry(i).getDoubleValue();
  return nominal;
}

const string& HOPSPACK::ParameterList::getParameter(const string& name, const char* nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  // Save nominal char* value as a string, and return the string value.
  tmpstrings.push_back(nominal);
  return tmpstrings.back();
}

const string& HOPSPACK::ParameterList::getParameter(const string& name, const string& nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();
  return nominal;
}
  
const HOPSPACK::Vector& HOPSPACK::ParameterList::getParameter(const string& name, const Vector& nominal) const
{
  ConstIterator i = params.find(name);
  if ((i != params.end()) && (entry(i).isVector()))
    return entry(i).getVectorValue();
  return nominal;
}
  

double HOPSPACK::ParameterList::getDoubleParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isDouble()))
    return entry(i).getDoubleValue();

  cerr << "HOPSPACK::ParameterList::getValueParameter - no such parameter (" << name << ")"<< endl;
  throw INTERNAL_ERROR;
}


//----------------------------------------------------------------------
// GetOrSet parameter methods
//----------------------------------------------------------------------

bool HOPSPACK::ParameterList::getOrSetParameter(const string& name,
                                                bool nominal)
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isBool()))
    return entry(i).getBoolValue();

  cerr << "HOPSPACK::ParameterList::getParameter - get error for bool" << endl;
  throw INTERNAL_ERROR;
}

int HOPSPACK::ParameterList::getOrSetParameter(const string& name,
                                               int nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isInt()))
    return entry(i).getIntValue();

  cerr << "HOPSPACK::ParameterList::getParameter - get error for int" << endl;
  throw INTERNAL_ERROR;
}

double HOPSPACK::ParameterList::getOrSetParameter(const string& name,
                                                  double nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isDouble()))
    return entry(i).getDoubleValue();

  cerr << "HOPSPACK::ParameterList::getParameter - get error for double" << endl;
  throw INTERNAL_ERROR;
}

const string& HOPSPACK::ParameterList::getOrSetParameter(const string& name,
                                                         const char* nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  cerr << "HOPSPACK::ParameterList::getParameter - get error for string" << endl;
  throw INTERNAL_ERROR;
}

const string& HOPSPACK::ParameterList::getOrSetParameter(const string& name,
                                                         const string& nominal) 
{
  ConstIterator i = params.find(name);

  if (i == params.end()) {
    params[name].setValue(nominal, true);
    i = params.find(name);
  }

  if ((i != params.end()) && (entry(i).isString()))
    return entry(i).getStringValue();

  cerr << "HOPSPACK::ParameterList::getParameter - get error for string" << endl;
  throw INTERNAL_ERROR;
}


//----------------------------------------------------------------------
// Get parameter methods without a nominal value
//----------------------------------------------------------------------

const vector< char > &  HOPSPACK::ParameterList::getCharVecParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isCharVec()))
    return entry(i).getCharVecValue();

  cerr << "HOPSPACK::ParameterList::getCharVecParameter - no such parameter (" << name << ")" << endl;
  throw INTERNAL_ERROR;
}

const HOPSPACK::Vector& HOPSPACK::ParameterList::getVectorParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isVector()))
    return entry(i).getVectorValue();

  cerr << "HOPSPACK::ParameterList::getVectorParameter - no such parameter (" << name << ")" << endl;
  throw INTERNAL_ERROR;
}

const HOPSPACK::Matrix& HOPSPACK::ParameterList::getMatrixParameter(const string& name) const
{
  ConstIterator i = params.find(name);

  if ((i != params.end()) && (entry(i).isMatrix()))
    return entry(i).getMatrixValue();

  cerr << "HOPSPACK::ParameterList::getMatrixParameter - no such parameter (" << name << ")" << endl;
  throw INTERNAL_ERROR;
}


//----------------------------------------------------------------------
// Parameter existence methods
//----------------------------------------------------------------------

bool HOPSPACK::ParameterList::isParameterBool(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isBool());

  return false;
}

bool HOPSPACK::ParameterList::isParameterInt(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isInt());

  return false;
}

bool HOPSPACK::ParameterList::isParameterDouble(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isDouble());

  return false;
}

bool HOPSPACK::ParameterList::isParameterString(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isString());

  return false;
}

bool HOPSPACK::ParameterList::isParameterCharVec(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isCharVec());

  return false;
}

bool HOPSPACK::ParameterList::isParameterSublist(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isList());

  return false;
}

bool HOPSPACK::ParameterList::isParameterVector(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isVector());

  return false;
}

bool HOPSPACK::ParameterList::isParameterMatrix(const string& name) const
{
  ConstIterator i = params.find(name);

  if (i != params.end())
    return (entry(i).isMatrix());

  return false;
}

bool HOPSPACK::ParameterList::isParameter(const string& name) const
{
  return (params.find(name) != params.end());
}


//----------------------------------------------------------------------
// Method deleteParameter
//----------------------------------------------------------------------
void HOPSPACK::ParameterList::deleteParameter(const string& name)
{
    Iterator  it = params.find (name);
    if (it != params.end())
        params.erase (it);

    return;
}


//----------------------------------------------------------------------
// Method print
//----------------------------------------------------------------------

ostream& HOPSPACK::ParameterList::print(ostream& stream, int indent) const
{
  if (params.begin() == params.end()) 
  {
    for (int j = 0; j < indent; j ++)
      stream << ' ';
    stream << "[empty list]" << endl;
  }
  else 
    for (ConstIterator i = params.begin(); i != params.end(); ++i) 
    {
      for (int j = 0; j < indent; j ++)
        stream << ' ';
      if (entry(i).isList()) 
      {
        stream << name(i) << " -> " << endl;
        entry(i).getListValue().print(stream, indent + 2);
      }
      else
        stream << name(i) << " = " << entry(i) << endl;
    }
  return stream;
}


//----------------------------------------------------------------------
// MPI pack and unpack methods
//----------------------------------------------------------------------

#if defined(HAVE_MPI)

    //---- DEFINE MESSAGE IDS FOR MPI/PVM BUFFERS.
    static const int  nMPIMSG_NEW_ENTRY   = 1;
    static const int  nMPIMSG_END_OF_LIST = 2;

void HOPSPACK::ParameterList::pack (GenProcComm &  cGPC) const
{
    for (ConstIterator i = params.begin(); i != params.end(); ++i)
    {
        cGPC.pack (nMPIMSG_NEW_ENTRY);
        cGPC.pack (name(i));
        entry(i).pack (cGPC);
    }
    cGPC.pack (nMPIMSG_END_OF_LIST);
    return;
}

void HOPSPACK::ParameterList::unpack (GenProcComm &  cGPC)
{
    int     code;
    string  name;
    cGPC.unpack (code);

    while (code != nMPIMSG_END_OF_LIST)
    {
        cGPC.unpack (name);
        params[name].unpack (cGPC);
        cGPC.unpack (code);
    }
    return;
}

#endif     //-- HAVE_MPI


//----------------------------------------------------------------------
// Private methods
//----------------------------------------------------------------------

#ifdef SNL_TFLOPS_ENV

const string& HOPSPACK::ParameterList::name(ConstIterator i) const
{
  return ((*i).first);
}

HOPSPACK::ParameterEntry& HOPSPACK::ParameterList::entry(Iterator i)
{
  return ((*i).second);
}

const HOPSPACK::ParameterEntry& HOPSPACK::ParameterList::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else

const string& HOPSPACK::ParameterList::name(ConstIterator i) const
{
  return (i->first);
}

HOPSPACK::ParameterEntry& HOPSPACK::ParameterList::entry(Iterator i)
{
  return (i->second);
}

const HOPSPACK::ParameterEntry& HOPSPACK::ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}

#endif
