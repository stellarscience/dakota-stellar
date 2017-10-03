// $Id: HOPSPACK_main_serial.cpp 220 2014-01-02 21:24:59Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-main/HOPSPACK_main_serial.cpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2014 Sandia Corporation
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
  @file HOPSPACK_main_serial.cpp
  @brief Main program that executes HOPSPACK as a serial program
         (no MPI, no multithreading).
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_EvaluatorFactory.hpp"
#include "HOPSPACK_ExecutorSerial.hpp"
#include "HOPSPACK_Hopspack.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_utils.hpp"


//----------------------------------------------------------------------
//  Static declarations
//----------------------------------------------------------------------

//---- RETURN CODES FOR main().
static const int  nSUCCESS =  0;
static const int  nERROR   = -1;


//--------------------------------------------------------------------
//! Main routine for serial version of HOPSPACK.
/*!
 *  The serial component that differs fundamentally from an MPI or multithreaded
 *  implementation is the Executor.  The serial Executor calls the Evaluator
 *  for the next queued point point and waits for it to complete.  This result
 *  is exchanged with the conveyor before the next point is evaluated.
 *
 *  The main loop of work is performed by the Mediator in Hopspack::solve().
 *
 *  @param[in] nArgc   Number of command line arguments (typically, 1).
 *  @param[in] saArgv  Command line arguments (typically, parameters file name).
 */
//--------------------------------------------------------------------
int  main (const int           nArgc,
           const char * const  saArgv[])
{
    using HOPSPACK::parseTextInputFile;   //-- FROM HOPSPACK_utils.hpp
    using HOPSPACK::ParameterList;
    using HOPSPACK::Evaluator;
    using HOPSPACK::EvaluatorFactory;
    using HOPSPACK::ExecutorSerial;
    using HOPSPACK::Hopspack;
    using HOPSPACK::Print;


    if (nArgc < 2)
    {
        cerr << "ERROR: Need the input file name." << endl;
        cout << "Usage:  HOPSPACK_main_serial <input file>" << endl;
        cout << "  The input file contains HOPSPACK parameters in text form."
             << endl;
        return( nERROR );
    }
    ParameterList  cParams;
    if (parseTextInputFile (saArgv[1], cParams) == false)
        return( nERROR );

    //---- CONSTRUCT AN EVALUATOR.
    Evaluator *  pEvaluator
        = EvaluatorFactory::newInstance (cParams.sublist ("Evaluator"));
    if (pEvaluator == NULL)
    {
        cerr << "ERROR: Could not construct Evaluator." << endl;
        return( nERROR );
    }

    //---- CONSTRUCT A SERIAL EXECUTOR.
    ExecutorSerial *  pExecutor = new ExecutorSerial (pEvaluator);

    //---- CONSTRUCT THE OPTIMIZER, CONFIGURE IT, AND RUN IT.
    Hopspack  optimizer (pExecutor);
    if (optimizer.setInputParameters (cParams) == true)
    {
        if (Print::doPrint (Print::FINAL_SOLUTION))
        {
            cout << endl << "Begin solving using the serial HOPSPACK executable."
                 << endl;
        }
        optimizer.solve();
    }

    delete pEvaluator;
    delete pExecutor;

    return( nSUCCESS );
}
