// $Id: HOPSPACK_utils.hpp 203 2012-05-14 22:27:30Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_utils.hpp $ 

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
  \file HOPSPACK_utils.hpp
  \brief No classes - declare utility functions in the HOPSPACK namespace.
*/

#ifndef HOPSPACK_UTILS_HPP
#define HOPSPACK_UTILS_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ParameterList.hpp"

namespace HOPSPACK
{

/*!  Get the next string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the next string (i.e., the
  next contiguous block of non-space characters). This is an empty
  string if no string is found.

  \retval Returns true if the string is successfully found,
  false otherwise.
  
*/
bool getNextString(const string& line, string::size_type& pos, string& value);

/*!  Get the next string on the given line, starting at
  position pos, and convert it to a double. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error in reading the
  next string, this is set to std::string::npos upon return.

  \param value - On output, filled in with the double value constrained
  in the next string (i.e., the next contiguous block of non-space
  characters).
  
  \retval Returns true if the next string contains a double,
  false otherwise.
  
*/
bool getNextDouble(const string& line, string::size_type& pos, double& value);


/*!  Get the next string on the given line, starting at
  position pos, and convert it to a int. 

  \param[in] line - Line of text from which to read

  \param[in,out] pos - On input, the starting position in the line. On output,
  the next position after the string (which may be
  std::string::npos). If there is any sort of error in reading the
  next string, this is set to std::string::npos upon return.

  \param[out] value - On output, filled in with the int value constrained
  in the next string (i.e., the next contguous block of non-space
  characters).
  
  \retval Returns true if the next string contains an integer, false otherwise.
*/
bool getNextInt(const string& line, string::size_type& pos, int& value);


/*!
  \brief Helper function for parseTextInputFile.
 */
bool processTextInputFileLine(const string& line,
                              ParameterList& params,
                              ParameterList*& subPtr,
                              ifstream &fin);


/*!
  \brief Parse an HOPSPACK input file and store the data in the given parameter list.

  \param filename - The file name.
  
  \param params - The parameter list that is to be filled in by this
  function

  \return Returns false if there are any problems parsing the input
  file, true otherwise.


*/
bool parseTextInputFile(const string filename, ParameterList& params);


}

#endif

