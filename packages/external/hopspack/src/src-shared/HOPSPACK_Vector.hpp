// $Id: HOPSPACK_Vector.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_Vector.hpp $ 

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
  \file HOPSPACK_Vector.hpp
  \brief Class declaration of HOPSPACK::Vector.
*/

#ifndef HOPSPACK_VECTOR_HPP
#define HOPSPACK_VECTOR_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK
{

//! Extends properties of Standard Template Library (STL) vector<double>.
class Vector
{

public:

  //! Constructor -- creates vector of length 0.
  Vector();

  //! Constructor -- creates vector of length n.
  Vector(int n);

  //! Constructor -- creates vector of length n, each element assigned the value val.
  Vector(int n, double val);

  //! Constructor -- creates a vector by copying the first n elements pointed by x. 
  Vector(int n, double *x);

  //! Copy constructor.
  Vector(const Vector& x);
  
  //! Copy constructor.
  Vector(const vector<double>& x);

  //! Copies x.
  Vector& operator=(const Vector& x);  
  
  //! Destructor.
  ~Vector();

  //@{ \name Memory allocating/altering methods

  //! Resizes vector length.
  void resize(int n);
  
  //! Sets the capacity of the vector to at least n.
  void reserve(int n);

  //! Appends a to the end of the vector, increasing vector length by 1.
  void push_back(double a);
  
  //! Reset vector to size n and all values to value alpha
  void assign(int n, double alpha);
  
  //! Appends x to end of the vector.
  void append(const Vector& x);

  //! Appends n copies of alpha to the end of the vector.
  void append(int n, double alpha);

  //! Deletes the ith element from the vector, decreasing the length by 1.
  void erase(int i);

  //@}

  //@{ \name Data access methods
  //! Returns the length of the vector.
  int size() const;

  //! Return true if vector is size zero, false otherwise.
  bool empty() const;

  //! Returns an equivalent STL vector.
  const vector<double>& getStlVector() const;

  //! Returns the ith element of the vector.
  double operator[](int i) const;
  
  //! Returns true if two vectors are equal, false otherwise.
  bool operator==(const Vector& x) const;

  //! Returns true if two vectors are not equal, false otherwise.
  bool operator!=(const Vector& x) const;

  //@}

  //@{ \name Vector mathematics

  //! Sets all entries of the vector to zero.
  void zero();

  //! Scales all elements of the vector by scalar alpha.
  void scale(double alpha);  
  
  //! On exit, the ith element is scaled by s[i].
  void scale(const Vector& s);  
    
  //! Returns the two norm of the vector.
  double norm() const;

  //! Returns the smallest element of the vector.
  double minElement() const;

  //! Returns the largest element of the vector.
  double maxElement() const;

  //! Compute the dot product between v (this vector) and x.
  double dot(const Vector& x) const;

  //! Adds x to this vector.
  Vector& operator+=(const Vector& x);
  
  //! Substracts x from this vector.
  Vector& operator-=(const Vector& x);

  //! Returns a reference to the ith element.
  double& operator[](int i);

  //@}

  //@{ \name Printing

  //! Prints out vector to specified stream.
  void leftshift(ostream& stream) const;

  //! Prints out vector to specified stream with given precision.
  /*!
   *  @param stream IOstream to send output.
   *  @param precision Number of digits after the decimal
   */
  void leftshift(ostream& stream, int precision) const;

  //@}

private:
  
  //! The vector.
  vector<double> vec;

};

}

#endif
