// $Id: HOPSPACK_Print.hpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_Print.hpp $ 

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
  \file HOPSPACK_Print.hpp
  \brief Class declaration of HOPSPACK::Print.
*/

#ifndef HOPSPACK_PRINT_HPP
#define HOPSPACK_PRINT_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK 
{

//! Printing utilities
class Print 
{

public:

  //! Enumerated type that specifies framework information to print.
  /*!
    If display = n, then everything with a value greater than or equal to n
    is printed.  For example, if display = 3, then the final solution,
    input parameters, and evaluated points are displayed.
  */
  enum PrintType 
    {
      //! Final Solution
      FINAL_SOLUTION = 1,

      //! Input parameters
      INPUT_PARAMETERS = 2,

      //! All evaluated trial points
      EVALUATED_POINTS = 3,

      //! All trial points from citizens
      UNEVALUATED_POINTS = 4,

      //! Trace everything
      MOST_VERBOSE = 5
    };
  //! Conveyor queue lists
  static PrintType  QUEUE_LISTS;


  //! Set static print parameter.
 /*!
  *  @param nDisplayLevel  Display level for framework.  Ignored if not a
  *                        defined PrintType value.
  */
  static void  setDisplayParameter (const PrintType  nDisplayLevel);

  //! Set static print parameter.
  static void  setPrecisionParameter (const int  nPrecision);

  //! Return true if nPrintLevel should be printed based on the current display parameter.
  static bool doPrint(const PrintType nPrintLevel);

  //! Return print precision (number of digits after the decimal).
  static int  getPrecision (void);


  private:

  //! The class is a singleton, so default constructor is hidden.
  Print (void);

  //! Display level
  static PrintType _nDisplayLevel;

  //! Precision for output of real numbers 
  static int precision;		

};

}

#endif
