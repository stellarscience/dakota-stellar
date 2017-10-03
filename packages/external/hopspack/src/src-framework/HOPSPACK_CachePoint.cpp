// $Id: HOPSPACK_CachePoint.cpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_CachePoint.cpp $ 

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
  \file HOPSPACK_CachePoint.cpp
  \brief Implement HOPSPACK::CachePoint
*/

#include "HOPSPACK_CachePoint.hpp" 
#include "HOPSPACK_ScaledComparison.hpp"


HOPSPACK::CachePoint::CachePoint (void) :
    xPtr(NULL),
    x(*xPtr),
    f(),
    cEqs(),
    cIneqs()
{
}

HOPSPACK::CachePoint::CachePoint (const Vector& x_in) :
    xPtr(NULL),
    x(x_in),
    f(),
    cEqs(),
    cIneqs()
{
}

HOPSPACK::CachePoint::CachePoint (const Vector& x_in,
                                  const Vector& f_in,
                                  const Vector& cEqs_in,
                                  const Vector& cIneqs_in) :
    xPtr(NULL),
    x(x_in),
    f(f_in),
    cEqs(cEqs_in),
    cIneqs(cIneqs_in)
{
}

HOPSPACK::CachePoint::CachePoint(const CachePoint& source) :
  xPtr(new Vector(source.x)),
  x(*xPtr),
  f(source.f),
  cEqs(source.cEqs),
  cIneqs(source.cIneqs)
{
}

HOPSPACK::CachePoint::~CachePoint()
{
  delete xPtr;
}

void HOPSPACK::CachePoint::copyData(const CachePoint& source)
{
  f = source.f;
  cEqs = source.cEqs;
  cIneqs = source.cIneqs;
}

const HOPSPACK::Vector& HOPSPACK::CachePoint::getF()
{
  return f;
}

const HOPSPACK::Vector& HOPSPACK::CachePoint::getEqs()
{
  return cEqs;
}

const HOPSPACK::Vector& HOPSPACK::CachePoint::getIneqs()
{
  return cIneqs;
}

bool HOPSPACK::CachePoint::operator>(const CachePoint& pt) const
{
    return( ScaledComparison::isGreaterThan (x, pt.x) );
}

bool HOPSPACK::CachePoint::operator<(const CachePoint& pt) const
{
    return( ScaledComparison::isLessThan (x, pt.x) );
}

bool HOPSPACK::CachePoint::operator!=(const CachePoint& pt) const
{
    return( ScaledComparison::isNotEqual (x, pt.x) );
}
