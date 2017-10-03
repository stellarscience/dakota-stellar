// $Id: HOPSPACK_ParameterEntry.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ParameterEntry.hpp $ 

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
  \file HOPSPACK_ParameterEntry.hpp
  \brief Class declaration of HOPSPACK::ParameterEntry.
*/

#ifndef HOPSPACK_PARAMETER_ENTRY_HPP
#define HOPSPACK_PARAMETER_ENTRY_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Vector.hpp"
#include "HOPSPACK_Matrix.hpp"

#if defined(HAVE_MPI)
  #include "HOPSPACK_GenProcComm.hpp"
#endif


namespace HOPSPACK {

class ParameterList; // another parameter type (forward declaration)

//! Contains a particular user parameter.
/*! 
  Stores a List of parameters \b or a single parameter, which can be
  a bool, int, double, string, HOPSPACK::Vector.
  For any single parameter, two attributes are tracked: 
  <ol>
  <li> Used (isGotten): True if it has been accessed by a "get" function (mutable)
  <li> Default (isSetByGet): True if it was set by a "get" function
  </ol>

  The "Used" parameter is mutable, which means it can be changed even
  for a const entry.
 */
class ParameterEntry {

public:

  /** \name Constructors / Destructor

  The "Used" parameter is set to false.
  The "Default" parameter is set to false, unless the optional second argument is true.

  */
  //@{ 
  
  //! Create an empty entry that doesn't contain anything
  ParameterEntry();
  
  //! Copy constructor
  ParameterEntry(const ParameterEntry& source);

  //! Create an entry containing a bool. 
  ParameterEntry(bool value, bool isCreatedByGet = false);

  //! Create an entry containing an int. 
  ParameterEntry(int value, bool isCreatedByGet = false);

  //! Create an entry containing a double.
  ParameterEntry(double value, bool isCreatedByGet = false);

  //! Create an entry containing a string.
  ParameterEntry(const string& value, bool isCreatedByGet = false);

  //! Create an entry containing a Vector. 
  ParameterEntry(const Vector& value, bool isCreatedByGet = false);

  //! Create an entry containing a Matrix.
  ParameterEntry(const Matrix& value, bool isCreatedByGet = false);

  //! Destructor
  ~ParameterEntry();

  //@}

  //@{ \name Copy

  //! Copy an entry
  ParameterEntry& operator=(const ParameterEntry& source);

  //@}
  
  /** @name Sublists
   *
   * Functions for handling parameters that are themselves lists.  */
  //@{
  //! Create a sublist. Sets "Default" according to the optional argument. Sets "Used" to false.
  ParameterList& setList(bool isCreatedByGet = false);
  //! Extract the reference to the list that is stored in the entry. Set "Used" to true.
  ParameterList& getListValue();
  //! Extract the reference to the list that is stored in the entry. Set "Used" to true.
  const ParameterList& getListValue() const;
  //@}

  /** @name Set functions. 

  Similar to the corresponding constructors. Completely erases the old entry.
  The "Used" parameter is (re-)set to false.
  The "Default" parameter is set to false, unless the optional second argument is true.

  */
  //@{ 
  //! Set bool value
  void setValue(bool value, bool isCreatedByGet = false);
  //! Set int value
  void setValue(int value, bool isCreatedByGet = false);
  //! Set double value
  void setValue(double value, bool isCreatedByGet = false);
  //! Set string value
  void setValue(const char* value, bool isCreatedByGet = false);
  //! Set string value
  void setValue(const string& value, bool isCreatedByGet = false);
  //! Set a character vector parameter
  void setValue(const vector< char >  value, bool isCreatedByGet = false);
  //! Set Vector value
  void setValue(const Vector& value, bool isCreatedByGet = false);
  //! Set Matrix value
  void setValue(const Matrix& value, bool isCreatedByGet = false);
  //@}

  /** @name Is functions. 
   
    Return true if the parameter is of the specified type; otherwise,
    return false. Has no affect on "Default" or "Used" status.
  */
  //@{ 
  //! Return true for bool entry
  bool isBool() const;
  //! Return true for int entry
  bool isInt() const;
  //! Return true for double entry
  bool isDouble() const;
  //! Return true for string entry
  bool isString() const;
  //! Return true for character vector entry
  bool isCharVec() const;
  //! Return true for list entry (i.e., a sublist)
  bool isList() const;
  //! Return true for Vector entry
  bool isVector() const;
  //! Return true for Matrix entry
  bool isMatrix() const;
  //@}
  
  /** @name Get functions. 
   
    Returns the value of the entry. Sets "Used" to true. Does not
    affect "Default" flag.
  */
  //@{ 
  //! Get bool value
  bool getBoolValue() const;
  //! Get int value
  int getIntValue() const;
  //! Get double value
  double getDoubleValue() const;
  //! Get string value
  const string& getStringValue() const;
  //! Get character vector value
  const vector< char > &  getCharVecValue() const;
  //! Get Vector value
  const Vector& getVectorValue() const;
  //! Get Matrix value
  const Matrix& getMatrixValue() const;
  //@}


#if defined(HAVE_MPI)
  /** @name Pack/Unpack for Inter-Process Communication

  These methods work in conjunction with MPI or PVM objects.
  */

  //@{

  //! Pack the entry for MPI transmission
  void  pack (GenProcComm &) const;

  //! Unpack the entry from MPI transmission
  void  unpack (GenProcComm &);
  
  //@}
#endif     //-- HAVE_MPI


  //@{ \name Printing

  //! Output the parameter to the given stream. 
  /*! 
    Formats the output as "<type,value>", except in the case of a
    list which just outputs "\<sublist\>". If the parameter has not yet
    been set, it outputs "\<NONE\>". This is the function called by the
    ostream operator<<. 
  */
  ostream& leftshift(ostream& stream) const;

  //@}

private:

  //! Reset the entry
  void reset();
  
  //! All possible parameter types that this class can store
  enum EntryType { 
    //! No entry type set yet (will be set later by setValue()
    HOPSPACK_NONE, 
    //! Boolean
    HOPSPACK_BOOL, 
    //! Integer
    HOPSPACK_INT, 
    //! Double
    HOPSPACK_DOUBLE, 
    //! String
    HOPSPACK_STRING,
    //! Character vector
    HOPSPACK_CHARVECTOR,
    //! Sublist (HOPSPACK::ParameterList)
    HOPSPACK_LIST,
    //! HOPSPACK::Vector
    HOPSPACK_VECTOR,
    //! HOPSPACK::Matrix
    HOPSPACK_MATRIX
  };

  //! Type of parameter stored in this object.
  EntryType type;

  //! Boolean value, if this is of type BOOL
  bool bval;

  //! Integer value, if this is of type INT
  int ival;

  //! Double value, if this is of type DOUBLE
  double dval;

  //! String value, if this is of type STRING
  string sval;

  //! Character vector value, if this is of type CHARVECTOR
  vector< char > cvval;

  //! Pointer to list, if this is of type LIST
  ParameterList* lval;		

  //! Vector value
  Vector vectorval;

  //! Matrix value
  Matrix matrixval;

  //! True if this parameter been accessed by a "get" function
  mutable bool isGotten;

  //! True if this parameter is a nominal value assigned by a "get" function
  mutable bool isSetByGet;

};

} // namespace HOPSPACK

//! Output the parameter. Relies of leftshift operator defined in the class.
ostream& operator<<(ostream& stream, const HOPSPACK::ParameterEntry& e);

#endif
