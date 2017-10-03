// $Id: HOPSPACK_ParameterEntry.cpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ParameterEntry.cpp $ 

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
  \file HOPSPACK_ParameterEntry.cpp
  \brief Implement HOPSPACK::ParameterEntry.
*/

#include "HOPSPACK_ParameterEntry.hpp"
#include "HOPSPACK_ParameterList.hpp"

HOPSPACK::ParameterEntry::ParameterEntry() : 
  type(HOPSPACK_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(false)
{
}

HOPSPACK::ParameterEntry::ParameterEntry(const ParameterEntry& source) :
  type(HOPSPACK_NONE),
  bval(false),
  ival(0),
  dval(0),
  sval(""), 
  cvval(), 
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(false)
{
  operator=(source);
}

HOPSPACK::ParameterEntry& HOPSPACK::ParameterEntry::operator=(const ParameterEntry& source)
{
  if (&source == this)
    return *this;

  reset();

  type = source.type;
  bval = source.bval;
  ival = source.ival;
  dval = source.dval;
  sval = source.sval;
  cvval = source.cvval;
  
  if ((type == HOPSPACK_LIST) && (source.lval != NULL)) 
  {
    lval = new ParameterList(*source.lval);
  }
  
  vectorval = source.vectorval;
  matrixval = source.matrixval;

  isGotten = source.isGotten;
  isSetByGet = source.isSetByGet;

  return *this;
}

HOPSPACK::ParameterEntry::ParameterEntry(bool value, bool isCreatedByGet) : 
  type(HOPSPACK_BOOL),
  bval(value),
  ival(0),
  dval(0),
  sval(""),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(isCreatedByGet)
{
}

HOPSPACK::ParameterEntry::ParameterEntry(int value, bool isCreatedByGet) : 
  type(HOPSPACK_INT),
  bval(false),
  ival(value),
  dval(0),
  sval(""),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

HOPSPACK::ParameterEntry::ParameterEntry(double value, bool isCreatedByGet) : 
  type(HOPSPACK_DOUBLE),
  bval(false),
  ival(0),
  dval(value),
  sval(""),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

HOPSPACK::ParameterEntry::ParameterEntry(const string& value, bool isCreatedByGet) : 
  type(HOPSPACK_STRING),
  bval(false),
  ival(0),
  dval(0),
  sval(value),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

HOPSPACK::ParameterEntry::ParameterEntry(const Vector& value, bool isCreatedByGet) : 
  type(HOPSPACK_VECTOR),
  bval(false),
  ival(0),
  dval(0),
  sval("" ),
  cvval(),
  lval(NULL),
  vectorval(value),
  matrixval(),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

HOPSPACK::ParameterEntry::ParameterEntry(const Matrix& value, bool isCreatedByGet) : 
  type(HOPSPACK_VECTOR),
  bval(false),
  ival(0),
  dval(0),
  sval("" ),
  cvval(),
  lval(NULL),
  vectorval(),
  matrixval(value),
  isGotten(false),
  isSetByGet(isCreatedByGet) 
{
}

HOPSPACK::ParameterEntry::~ParameterEntry() 
{
  reset();
}

void HOPSPACK::ParameterEntry::reset()
{
  type = HOPSPACK_NONE;

  delete lval;
  lval = NULL;

  isGotten = false;
  isSetByGet = false;
}

void HOPSPACK::ParameterEntry::setValue(bool value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_BOOL;
  bval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(int value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_INT;
  ival = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(double value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_DOUBLE;
  dval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(const char* value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(const string& value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_STRING;
  sval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(const vector< char > value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_CHARVECTOR;
  cvval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(const Vector& value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_VECTOR;
  vectorval = value;
  isSetByGet = isCreatedByGet;
}

void HOPSPACK::ParameterEntry::setValue(const Matrix& value, bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_MATRIX;
  matrixval = value;
  isSetByGet = isCreatedByGet;
}

HOPSPACK::ParameterList& HOPSPACK::ParameterEntry::setList(bool isCreatedByGet)
{
  reset();
  type = HOPSPACK_LIST;
  lval = new ParameterList();
  isSetByGet = isCreatedByGet;
  isGotten = true;
  return *lval;
}


bool HOPSPACK::ParameterEntry::isBool() const
{
  return (type == HOPSPACK_BOOL);
}

bool HOPSPACK::ParameterEntry::isInt() const
{
  return (type == HOPSPACK_INT);
}

bool HOPSPACK::ParameterEntry::isDouble() const
{
  return (type == HOPSPACK_DOUBLE);
}

bool HOPSPACK::ParameterEntry::isString() const
{
  return (type == HOPSPACK_STRING);
}

bool HOPSPACK::ParameterEntry::isCharVec() const
{
  return (type == HOPSPACK_CHARVECTOR);
}

bool HOPSPACK::ParameterEntry::isList() const
{
  return (type == HOPSPACK_LIST);
}

bool HOPSPACK::ParameterEntry::isMatrix() const
{
  return (type == HOPSPACK_MATRIX);
}

bool HOPSPACK::ParameterEntry::isVector() const
{
  return (type == HOPSPACK_VECTOR);
}

bool HOPSPACK::ParameterEntry::getBoolValue() const
{
  if (type != HOPSPACK_BOOL)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getBoolValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return bval;
}

int HOPSPACK::ParameterEntry::getIntValue() const
{
  if (type != HOPSPACK_INT)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getIntValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return ival;
}

double HOPSPACK::ParameterEntry::getDoubleValue() const
{
  if (type != HOPSPACK_DOUBLE)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getDoubleValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return dval;
}

const string& HOPSPACK::ParameterEntry::getStringValue() const
{
  if (type != HOPSPACK_STRING)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getStringValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return sval;
}

const vector< char > &  HOPSPACK::ParameterEntry::getCharVecValue() const
{
  if (type != HOPSPACK_CHARVECTOR)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getCharVecValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return cvval;
}

HOPSPACK::ParameterList& HOPSPACK::ParameterEntry::getListValue() 
{
  if (type != HOPSPACK_LIST)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getListValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return *lval;
}

const HOPSPACK::ParameterList& HOPSPACK::ParameterEntry::getListValue() const
{
  if (type != HOPSPACK_LIST)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getListValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return *lval;
}

const HOPSPACK::Vector& HOPSPACK::ParameterEntry::getVectorValue() const
{
  if (type != HOPSPACK_VECTOR)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getVectorValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return vectorval;
}

const HOPSPACK::Matrix& HOPSPACK::ParameterEntry::getMatrixValue() const
{
  if (type != HOPSPACK_MATRIX)
  {
      cerr << "ERROR: Requested wrong parameter type"
           << "  <ParameterEntry::getMatrixValue()>" << endl;
      throw INTERNAL_ERROR;
  }
  isGotten = true;
  return matrixval;
}

ostream& HOPSPACK::ParameterEntry::leftshift(ostream& stream) const
{
  switch(type) {
  case HOPSPACK_BOOL: 
    stream << (bval ? "true" : "false");
    break;
  case HOPSPACK_INT:
    stream << ival;
    break;
  case HOPSPACK_DOUBLE:
    stream << dval;
    break;
  case HOPSPACK_STRING:
    stream << "\"" << sval << "\"";
    break;
  case HOPSPACK_CHARVECTOR:
    for (int  i = 0; i < (int) cvval.size(); i++)
      cout << cvval[i] << ' ';
    break;
  case HOPSPACK_LIST:
    break;
  case HOPSPACK_VECTOR:
    vectorval.leftshift (stream);
    break;
  case HOPSPACK_MATRIX:
    matrixval.formattedPrint ("    ", stream);
    break;
  default:
    stream << "(empty non-typed parameter)";
    break;
  }

  if (isSetByGet)
    stream << "   [default]";
  else if (!isGotten)
    stream << "   [unused]";
  

  return stream;
}


#if defined(HAVE_MPI)

void  HOPSPACK::ParameterEntry::pack (GenProcComm &  cGPC) const
{
  switch(type) {
  case HOPSPACK_BOOL: 
    cGPC.pack(HOPSPACK_BOOL);
    cGPC.pack(bval);
    break;
  case HOPSPACK_INT:
    cGPC.pack(HOPSPACK_INT);
    cGPC.pack(ival);
    break;
  case HOPSPACK_DOUBLE:
    cGPC.pack(HOPSPACK_DOUBLE);
    cGPC.pack(dval);
    break;
  case HOPSPACK_STRING:
    cGPC.pack(HOPSPACK_STRING);
    cGPC.pack(sval);
    break;
  case HOPSPACK_CHARVECTOR:
    //---- SHOULD NOT NEED THIS.
    cerr << "HOPSPACK::ParameterEntry::pack - No charvec implementation";
    throw INTERNAL_ERROR;
  case HOPSPACK_LIST:
    cGPC.pack(HOPSPACK_LIST);
    lval->pack(cGPC);
    break;
  case HOPSPACK_VECTOR:
    cGPC.pack(HOPSPACK_VECTOR);
    cGPC.pack(vectorval);
    break;
  case HOPSPACK_MATRIX:
    //---- SHOULD NOT NEED THIS.
    cerr << "HOPSPACK::ParameterEntry::pack - No matrix implementation";
    throw INTERNAL_ERROR;
  default:
    cerr << "HOPSPACK::ParameterEntry::pack - Empty non-typed parameter";
    throw INTERNAL_ERROR;
    break;
  }

  cGPC.pack(isGotten);
  cGPC.pack(isSetByGet);
}

void  HOPSPACK::ParameterEntry::unpack (GenProcComm &  cGPC)
{
  int itype;
  cGPC.unpack(itype);

  switch(itype) {
  case HOPSPACK_BOOL: 
    type = HOPSPACK_BOOL;
    cGPC.unpack(bval);
    break;
  case HOPSPACK_INT:
    type = HOPSPACK_INT;
    cGPC.unpack(ival);
    break;
  case HOPSPACK_DOUBLE:
    type = HOPSPACK_DOUBLE;
    cGPC.unpack(dval);
    break;
  case HOPSPACK_STRING:
    type = HOPSPACK_STRING;
    cGPC.unpack(sval);
    break;
  case HOPSPACK_CHARVECTOR:
    //---- SHOULD NOT NEED THIS.
    cerr << "HOPSPACK::ParameterEntry::pack - No charvec implementation";
    throw INTERNAL_ERROR;
  case HOPSPACK_LIST:
    type = HOPSPACK_LIST;
    lval = new ParameterList();
    lval->unpack(cGPC);
    break;
  case HOPSPACK_VECTOR:
    type = HOPSPACK_VECTOR;
    cGPC.unpack(vectorval);
    break;
  case HOPSPACK_MATRIX:
    //---- SHOULD NOT NEED THIS.
    cerr << "HOPSPACK::ParameterEntry::pack - No matrix implementation";
    throw INTERNAL_ERROR;
  default:
    cerr << "HOPSPACK::ParameterEntry::pack - Empty non-typed parameter";
    throw INTERNAL_ERROR;
    break;
  }

  cGPC.unpack(isGotten);
  cGPC.unpack(isSetByGet);
}
#endif     //-- HAVE_MPI


ostream& operator<<(ostream& stream, const HOPSPACK::ParameterEntry& e)
{
  return e.leftshift(stream);
}


