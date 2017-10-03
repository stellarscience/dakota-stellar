// $Id: HOPSPACK_CacheManager.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_CacheManager.hpp $ 

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

/*! \file HOPSPACK_CacheManager.hpp
    \brief Class declaration for HOPSPACK::CacheManager
*/

#ifndef HOPSPACK_CACHEMANAGER_HPP
#define HOPSPACK_CACHEMANAGER_HPP

#include "HOPSPACK_CachePoint.hpp" 
#include "HOPSPACK_CacheSplayTree.hpp" 
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK 
{

/*! 
  \brief Manages cached points efficiently.

  The cache supports the ability to read and/or write cached function
  values from files.  Each line corresponds to one point and its values.
  A typical single line may look like the following.

  \verbatim
  x=[ 0.0000e+00 -1.2500e-02  ] f=[ 3.1250e-04 ] ceq=[ (empty)] cineq=[ (empty) ]
  \endverbatim

  The line is parsed as follows:
  <ul>
  <li> The line must begin with "x=[".
  <li> This is followed by a minimum of one space.
  <li> Entries of x, separated by spaces, are read until "]" is
       encountered, surrounded by white space.
  <li> This is followed by a minimum of one space.
  <li> The next characters are "f=[".
  <li> This is followed by a minimum of one space.
  <li> Entries of f, separated by spaces, are read until "]" is
       encountered, surrounded by white space.  The single entry "(empty)"
       is allowed, meaning that f has no value.
  <li> The next characters are "ceq=[".
  <li> Another vector of entries follows for nonlinear equality constraints,
       following the same format as "f=[".
  <li> The next characters are "cineq=[".
  <li> Another vector of entries follows for nonlinear inequality constraints,
       following the same format as "f=[".
  </ul>

  Any line that does not conform to the above format is ignored.
  Comments are not allowed, but empty lines can be inserted to delimit sections.

  \author H. Alton Patrick, Summer 2000<br>
  Todd Plantenga, Tamara G. Kolda
*/
class CacheManager 
{

public:

  /*! Constructor */
   /*!
   *  @param[in] params  User input from "Mediator" sublist.
   */
 CacheManager (const ParameterList  &  params);

  /*! Destructor */
  ~CacheManager();

  //! Add the given point to the cache.
  bool insert(const Vector& x,
              const Vector& f,
              const Vector& cEqs,
              const Vector& cIneqs);
  
  //! Return true if x is cached and fill in the function values.
  bool isCached(const Vector& x,
                      Vector& f,
                      Vector& cEqs,
                      Vector& cIneqs);

  //! Print debug information about the instance.
  void  printDebugInfo (void) const;

private:

  //! By design, there is no copy constructor.
  CacheManager (const CacheManager &);
  //! By design, there is no assignment operator.
  CacheManager & operator= (const CacheManager &);

  //! Parse the cache input file
  void parseInputFile(const string &  filename);

  //! Process a single line from the cache input file, extracting one point
  bool processInputLine(string& line);

  bool readVectorFromLine (string            &  line,
                           string::size_type &  line_pos,
                           Vector            &  result);

  //! Open the output file for the cache
  void openOutputFile(const string &  filename);

  //! Write a given cache point to the output file
  void writeToOutputFile(const Vector& x,
                         const Vector& f,
                         const Vector& cEqs,
                         const Vector& cIneqs);

  //! Close the output file
  void closeOutputFile();


  //! Pointer to splay tree containing the cache
  CacheSplayTree<CachePoint>* treeptr;   

  //! Use cache output file?
  bool isFout;

  //! Cache output file
  ofstream fout;

  string outname;

  string inname;
  bool bCanOpenInname;

  //! Precision of file output (digits after the decimal point)
  int precision;
};
}

#endif
