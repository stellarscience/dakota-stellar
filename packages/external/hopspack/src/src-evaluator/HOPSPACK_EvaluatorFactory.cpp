// $Id: TBD $
// $URL: svn+ssh://software.sandia.gov/svn/private/hopspack/trunk/src/src-evaluator/HOPSPACK_EvaluatorFactory.cpp $

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2014 Sandia Corporation
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
  @file HOPSPACK_EvaluatorFactory.cpp
  @brief Implement HOPSPACK::EvaluatorFactory.
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_EvaluatorFactory.hpp"
#include "HOPSPACK_SystemCall.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//  Constructor
//----------------------------------------------------------------------
EvaluatorFactory::EvaluatorFactory (void)
{
    return;
}


//----------------------------------------------------------------------
//  Destructor
//----------------------------------------------------------------------
EvaluatorFactory::~EvaluatorFactory (void)
{
    return;
}


//----------------------------------------------------------------------
//  Method newInstance
//----------------------------------------------------------------------
Evaluator *  EvaluatorFactory::newInstance
                 (const ParameterList &  cEvalParams)
{
    //---- USE THE TYPE PARAMETER TO CHOOSE THE CONCRETE IMPLEMENTATION.
    //---- THIS LEVEL OF INDIRECTION MAKES EVALUATORS EXTENDABLE AT RUNTIME.
    //----
    //---- THE PARAMETER PARSER ENSURES cEvalParams CONTAINS AT MOST ONE
    //---- EVALUATOR TYPE.
    string  sType = cEvalParams.getParameter ("Evaluator Type", "System Call");

    Evaluator *  pResult = NULL;

    if (sType == "System Call")
    {
        try
        {
            pResult = (Evaluator *) new SystemCall (cEvalParams);
        }
        catch (const char * const)     //-- THROWS STRING MESSAGE ON ERROR
        {
            pResult = NULL;
        }
    }
    else
    {
        //---- DO NOT RECOGNIZE THE TYPE OF EVALUATOR.
        cerr << "ERROR: The value '" << sType
             << "' in parameter 'Evaluator Type' is not recognized." << endl;
        cerr << "  Please change parameter 'Evaluator Type' in sublist "
             << "'Evaluator'." << endl;
    }

    return( pResult );
}


}     //-- namespace HOPSPACK
