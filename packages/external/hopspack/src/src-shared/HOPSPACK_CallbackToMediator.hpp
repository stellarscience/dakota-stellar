// $Id: HOPSPACK_CallbackToMediator.hpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_CallbackToMediator.hpp $

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
  @file HOPSPACK_CallbackToMediator.hpp
  @brief Interface declaration for HOPSPACK::CallbackToMediator.
*/

#ifndef HOPSPACK_CALLBACKTOMEDIATOR_HPP
#define HOPSPACK_CALLBACKTOMEDIATOR_HPP

#include "HOPSPACK_common.hpp"

namespace HOPSPACK
{

//! Forward declaration.
class Citizen;

//----------------------------------------------------------------------
//! Interface implemented by the Mediator for adding citizens.
/*! A parent citizen that can dynamically spawn a child citizen uses this
    interface to add the child to the Mediator.  The interface declaration
    is in the "shared" directory to avoid circular compile dependencies.
*/
//----------------------------------------------------------------------
class CallbackToMediator
{
  public:

    //! Destructor.
    virtual ~CallbackToMediator (void) = 0;


    //! Reserve and return a unique citizen ID.
    virtual int  reserveUniqueCitizenID (void) = 0;

    //! Add a dynamically created child citizen to the Mediator.
    /*!
     *  The Mediator that implements this interface will define the method.
     *
     *  @param[in] pCitizen   Fully constructed and validated citizen.
     *  @param[in] nParentID  Parent ID number.
     *  @return false         If the name already exists (and delete pCitizen).
     */
    virtual bool  addChildCitizen (      Citizen * const  pCitizen,
                                   const int              nParentID) = 0;

  protected:

    //! Default constructor for abstract class is visible only to subclasses.
    CallbackToMediator (void);

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_CALLBACKTOMEDIATOR_HPP
