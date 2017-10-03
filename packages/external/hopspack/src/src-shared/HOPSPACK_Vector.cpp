// $Id: HOPSPACK_Vector.cpp 217 2013-11-25 21:59:49Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_Vector.cpp $ 

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
  \file HOPSPACK_Vector.cpp
  \brief Implement HOPSPACK::Vector.
*/

#include <iomanip>
#include <math.h>                 //-- FOR sqrt

#include "HOPSPACK_float.hpp"
#include "HOPSPACK_LapackWrappers.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_Vector.hpp"


HOPSPACK::Vector::Vector()
{
}

HOPSPACK::Vector::Vector(int n) : vec(n)
{
}

HOPSPACK::Vector::Vector(int n, double val) : vec(n, val)
{
}

HOPSPACK::Vector::Vector(int n, double *x) : vec(x,x+n)
{
}

HOPSPACK::Vector::Vector(const HOPSPACK::Vector& x) : vec(x.vec)
{
}

HOPSPACK::Vector::Vector(const vector<double>& x) : vec(x)
{
}

HOPSPACK::Vector& HOPSPACK::Vector::operator=(const Vector& x)
{
  vec = x.vec;
  return *this;
}

HOPSPACK::Vector::~Vector()
{
}

void HOPSPACK::Vector::resize(int n)
{
  vec.resize(n);
}

void HOPSPACK::Vector::reserve(int n)
{
  vec.reserve(n);
}

void HOPSPACK::Vector::push_back(double a)
{
  vec.push_back(a);
}

void HOPSPACK::Vector::assign(int n, double alpha)
{
  vec.assign(n,alpha);
}

void HOPSPACK::Vector::append(const Vector& x)
{
  if (x.size() > 0)
    vec.insert(vec.end(), x.vec.begin(), x.vec.end());
}

void HOPSPACK::Vector::append(int n, double alpha)
{
  vec.insert(vec.end(), n, alpha);
}

void HOPSPACK::Vector::erase(int i)
{
  vec.erase(vec.begin()+i);
}


int HOPSPACK::Vector::size() const
{
    return( (int) vec.size() );
}

bool HOPSPACK::Vector::empty() const
{
  return vec.empty();
}

const vector<double>& HOPSPACK::Vector::getStlVector() const
{
  return vec;
}

double HOPSPACK::Vector::operator[](int i) const
{
  return vec[i];
}

bool HOPSPACK::Vector::operator==(const Vector& a) const
{
  return (vec == a.vec);
}

bool HOPSPACK::Vector::operator!=(const Vector& a) const
{
  return (vec != a.vec);
}

void  HOPSPACK::Vector::leftshift(ostream& stream) const
{
  leftshift (stream, Print::getPrecision());
  return;
}

void  HOPSPACK::Vector::leftshift(ostream& stream, int precision) const
{
  if (vec.size() == 0)
  {
    stream << "(empty)";
    return;
  }

  int  nPrec = precision;
  if (nPrec < 0)
    nPrec =  Print::getPrecision();

  stream.setf(ios::scientific);
  stream.precision(nPrec);

  for (int i = 0; i < (int) vec.size(); i++) 
  {
    double  d = vec[i];
    if (exists(d))
      #if defined(WIN32)
        //---- WINDOWS USES FORMAT SX.XXXeSXXX.
        stream << setw(nPrec + 8) << d << " ";
      #else
        //---- UNIX USES FORMAT    SX.XXXeSXX.
        stream << setw(nPrec + 7) << d << " ";
      #endif
    else
    {
      stream << " DNE";
      #if defined(WIN32)
        for (int k = 0; k < nPrec + 5; k++)
          stream << " ";
      #else
        for (int k = 0; k < nPrec + 4; k++)
          stream << " ";
      #endif
    }
  }

  stream.unsetf(ios::scientific);
  return;
}

void HOPSPACK::Vector::zero()
{
  vec.assign(vec.size(), 0);
}

void HOPSPACK::Vector::scale(double alpha)
{
  for (int i = 0; i < (int) vec.size(); i++)
    vec[i] = alpha * vec[i];
}

void HOPSPACK::Vector::scale(const Vector& s)
{
  if (s.vec.size() != vec.size())
  {
    cerr << "ERROR: Vector size mismatch  <HOPSPACK::Vector.scale()>" << endl;
    throw INTERNAL_ERROR;
  }
  
  for (int i = 0; i < (int) vec.size(); i++)
    vec[i] = vec[i] * s.vec[i];
}

double HOPSPACK::Vector::norm() const
{
  double norm=0;
  for (int i = 0; i < (int) vec.size(); i++)
    norm += vec[i]*vec[i];
  return sqrt(norm);
}

double HOPSPACK::Vector::minElement() const
{
  if (vec.empty())
  {
    cerr << "ERROR: Vector is empty  <HOPSPACK::Vector.minElement()>" << endl;
    throw INTERNAL_ERROR;
  }

  double minvi = vec[0];
  for (int i = 1; i < (int) vec.size(); i++)
    minvi = (minvi < vec[i]) ? minvi : vec[i];
  return minvi;
}

double HOPSPACK::Vector::maxElement() const
{
  if (vec.empty())
  {
    cerr << "ERROR: Vector is empty  <HOPSPACK::Vector.maxElement()>" << endl;
    throw INTERNAL_ERROR;
  }

  double maxvi = vec[0];
  for (int i = 1; i < (int) vec.size(); i++)
    maxvi = (maxvi > vec[i]) ? maxvi : vec[i];
  return maxvi;
}

double HOPSPACK::Vector::dot(const Vector& x) const
{
  if (x.vec.size() != vec.size())
  {
    cerr << "ERROR: Vector size mismatch  <HOPSPACK::Vector.dot()>" << endl;
    throw INTERNAL_ERROR;
  }
  
  double  z = LapackWrappers::getTheInstance().ddot (this->vec.size(),
                                                     &(this->vec[0]),
                                                     &(x.vec[0]));
  return z;
}

HOPSPACK::Vector& HOPSPACK::Vector::operator+=(const Vector& x)
{
  if (x.vec.size() != vec.size())
  {
    cerr << "ERROR: Vector size mismatch  <HOPSPACK::Vector.operator+=()>"
         << endl;
    throw INTERNAL_ERROR;
  }
  
  for (int i = 0; i < (int) vec.size(); i++)
    vec[i] += x.vec[i];

  return *this;
}

HOPSPACK::Vector& HOPSPACK::Vector::operator-=(const Vector& x)
{
  if (x.vec.size() != vec.size())
  {
    cerr << "ERROR: Vector size mismatch  <HOPSPACK::Vector.operator-=()>"
         << endl;
    throw INTERNAL_ERROR;
  }
  
  for (int i = 0; i < (int) vec.size(); i++)
    vec[i] -= x.vec[i];

  return *this;
}

double& HOPSPACK::Vector::operator[](int i)
{
  return vec[i];
}

