// $Id: HOPSPACK_ParameterList.hpp 203 2012-05-14 22:27:30Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_ParameterList.hpp $ 

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
  \file HOPSPACK_ParameterList.hpp
  \brief Class declaration of HOPSPACK::ParameterList.
*/

#ifndef HOPSPACK_PARAMETER_LIST_HPP
#define HOPSPACK_PARAMETER_LIST_HPP

#include <map>

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ParameterEntry.hpp"
#include "HOPSPACK_Matrix.hpp"

#if defined(HAVE_MPI)
  #include "HOPSPACK_GenProcComm.hpp"
#endif


namespace HOPSPACK {

//! Contains a set of user parameters.
/*! A ParameterList stores a set of parameters, optionally organized in
  sublists. Each parameter has a name and a ParameterEntry. The name is a string,
  case-sensitive.  The ParameterEntry is itself a ParameterList (a sublist)
  or a bool, int, double, string, Vector, or Matrix.
  Each ParameterEntry keeps track of whether it was set to be the default value
  or set by the input, and whether or not it has been used (i.e., accessed by
  a get method).
*/

class ParameterList {

public:


    //! Constructor
  ParameterList();

  //! Copy constructor (deep copy)
  ParameterList(const ParameterList& source);

  //! Copy (deep copy)
  ParameterList& operator=(const ParameterList& source);

  //! Destructor
  ~ParameterList();

  //-------------------------
  /** @name Sublists 

  The entries in a parameter list can be organized into sublists.

  */
  //@{

  //! Create a new sublist \b or return the sublist if it already exists.
  /*! 
    Creates an empty sublist and returns a reference to the
    sublist. If the list already exists, returns a reference to that
    sublist. If the name exists but is not a sublist, throws an error.
  */
  ParameterList& getOrSetSublist(const string& name);

  //! Returns a const reference to the sublist, or throws an error.
  /*! 
    If the name is not found, then an empty list is returned.
    If the name exists but is not a sublist, throws an error.
  */
  const ParameterList& sublist(const string& name) const;
  //@}
  

  //-------------------------
  /** @name Setting Parameters 

    Methods to set parameters in a ParameterList. Each parameter has a name and
    a value which is either bool, int, double, string, Vector, or Matrix.
    Any previous value is overwritten.

    If necessary, use static_cast<type>() to resolve the type when its
    ambiguous. For example,

    \code
    ParameterList list;
    list.setParameter("some name", static_cast<double>(0))
    \endcode

    creates a double parameter named "some name" with a value of zero.

    Both char* and string are stored as strings internally. 

    Sets "Default" to false.
  */
  //@{

  //! Set a bool parameter
  void setParameter(const string& name, bool value);

  //! Set an int parameter
  void setParameter(const string& name, int value);

  //! Set a double parameter
  void setParameter(const string& name, double value);

  //! Set a string parameter
  void setParameter(const string& name, const char* value);

  //! Set a string parameter
  void setParameter(const string& name, const string& value);

  //! Set a character vector parameter
  void setParameter(const string& name, const vector< char >  value);

  //! Set a Vector parameter
  void setParameter(const string& name, const Vector& value);

  //! Set a Matrix parameter
  void setParameter(const string& name, const Matrix& value);
  //@}


  /** @name Getting Parameters 
   
    Get the value of a parameter from the list. Returns the nominal
    value if the parameter is not already specified.

    Caller should use static_cast<type>() when the type is ambiguous. 

    Both char* and string map return string values.
  */
  //@{

  //! Get a bool parameter - no change to the list
  bool getParameter(const string& name, bool nominal) const;

  //! Get an int parameter - no change to the list
  int getParameter(const string& name, int nominal) const;

  //! Get a double parameter - no change to the list
  double getParameter(const string& name, double nominal) const;

  //! Get a string parameter - no change to the list
  const string& getParameter(const string& name, const char* nominal) const;

  //! Get a string parameter - no change to the list
  const string& getParameter(const string& name, const string& nominal) const;

  //! Get a Vector parameter - no change to the list
  const Vector& getParameter(const string& name, const Vector& nominal) const;

  //@}

  /** @name GetOrSet Parameters 
   
    Get the value of a parameter from the list, or add the nominal value
    to the list if the parameter is not already specified (also setting "Default"
    to true).

    Use static_cast<type>() when the type is ambiguous. 

    Both char* and string map return string values.
  */
  //@{

  //! Get a bool parameter
  bool getOrSetParameter(const string& name, bool nominal);

  //! Get an int parameter
  int getOrSetParameter(const string& name, int nominal);

  //! Get a double parameter
  double getOrSetParameter(const string& name, double nominal);

  //! Get a string parameter
  const string& getOrSetParameter(const string& name, const char* nominal);

  //! Get a string parameter
  const string& getOrSetParameter(const string& name, const string& nominal);

  //@}

  /** @name Getting Parameters Without Nominal Value

  Get the value of a parameter from the list. Throws an error if the parameter
  does not exist or is not the correct type.
  */
  //@{

  //! Get a double parameter - no change to the list
  double getDoubleParameter(const string& name) const;

  //! Get a character vector parameter - no change to the list
  const vector< char > &  getCharVecParameter(const string& name) const;

  //! Get a Vector parameter - no change to the list
  const Vector& getVectorParameter(const string& name) const;

  //! Get a Matrix parameter - no change to the list
  const Matrix& getMatrixParameter(const string& name) const;
  //@}
 
  /** @name Checking parameter existence.
   
    Returns true if the specified parameter exists AND is of the
    specified type.
  */
  //@{

  //! Return true if a parameter with this name exists.
  bool isParameter(const string& name) const;

  //! Return true if a bool parameter with this name exists
  bool isParameterBool(const string& name) const;

  //! Return true if an int parameter with this name exists
  bool isParameterInt(const string& name) const;

  //! Return true if a double parameter with this name exists
  bool isParameterDouble(const string& name) const;

  //! Return true if a string parameter with this name exists
  bool isParameterString(const string& name) const;

  //! Return true if a character vector parameter with this name exists
  bool isParameterCharVec(const string& name) const;

  //! Return true if a sublist with this name exists
  bool isParameterSublist(const string& name) const;

  //! Return true if a Value parameter with this name exists
  bool isParameterValue(const string& name) const;

  //! Return true if a Vector parameter with this name exists
  bool isParameterVector(const string& name) const;

  //! Return true if a Matrix parameter with this name exists
  bool isParameterMatrix(const string& name) const;
  //@}

  //@{ \name Deleting Parameters

  //! Delete the parameter if it exists.
  void deleteParameter(const string& name);

  //@}

  //@{ \name Printing

  //! Pretty-print a list. Indents sublists.
  ostream& print(ostream& stream = cout, int indent = 0) const;

  //@}


#if defined(HAVE_MPI)
  /** @name Pack/Unpack for Inter-Process Communication

  These methods work in conjunction with MPI or PVM objects.
  */

  //@{ 
    
  //! Pack for MPI transmission
  void  pack (GenProcComm &) const;

  //! Unpack from MPI transmission
  void  unpack (GenProcComm &);

  //@}
#endif     //-- HAVE_MPI


private:

  //! Parameter container typedef
  typedef map<string, ParameterEntry> Map;

  //! Parameter container const iterator typedef
  typedef Map::const_iterator ConstIterator;

  //! Parameter container iterator typedef
  typedef Map::iterator Iterator;


  //! Access to name (i.e., returns i->first)
  const string& name(ConstIterator i) const;

  //! Access to ParameterEntry (i.e., returns i->second)
  ParameterEntry& entry(Iterator i);

  //! Access to ParameterEntry (i.e., returns i->second)
  const ParameterEntry& entry(ConstIterator i) const;


  //! Parameter list
  Map params;
 
  //! Used to create a string when the getParameter is called with a
  //! char* nominal value. A new string is created for each such
  //! argument. The whole group of strings is destroyed when this object
  //! is destroyed. This is really annoying, but I don't know a better 
  //! way.
  mutable vector<string> tmpstrings;
};

} // namespace HOPSPACK

#endif


