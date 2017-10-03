// $Id: HOPSPACK_float.hpp 166 2010-03-22 19:58:07Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_float.hpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2010 Sandia Corporation
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
  \file HOPSPACK_float.hpp
  \brief No classes - declare functions for handling double precision values.
*/

#ifndef HOPSPACK_FLOAT_HPP
#define HOPSPACK_FLOAT_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK
{
    //! Use x = HOPSPACK::dne() to say that the value of x does not exist.
    /*!
     *  This is done by returning a special constant, typically DBL_MAX.
     */
    double  dne (void);

    //! Returns true if value was previously set by a call to dne().
    bool  exists (const double  value);

    //! Return the architecture constant for machine epsilon.
    double  getMachineEpsilon (void);

    //! Return a random number x such that 0 <= x < 1.
    /*!
     *  Generate using a linear congruential method for consistent, solid
     *  behavior on all platforms.
     */
    double  genRandomNumber (void);

    //! Return true if the double precision number is valid.
    /*!
     *  Invalid numbers (infinity or NaN) are the result of invalid
     *  arithmetic operations.
     */
    bool  isDoubleValid (const double  d);
}

#endif

