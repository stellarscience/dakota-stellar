// $Id: HOPSPACK_CallbackToMediator.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_CallbackToMediator.cpp $

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
  @file HOPSPACK_CallbackToMediator.cpp
  @brief Partially implement abstract class  HOPSPACK::CallbackToMediator.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CallbackToMediator.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Default constructor for Base class does nothing.
//----------------------------------------------------------------------
CallbackToMediator::CallbackToMediator (void)
{
    return;
}


//----------------------------------------------------------------------
//  Base class needs a destructor even though it's declared pure virtual.
//----------------------------------------------------------------------
CallbackToMediator::~CallbackToMediator (void)
{
    return;
}


}     //-- namespace HOPSPACK
