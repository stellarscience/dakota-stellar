// $Id: HOPSPACK_main_threaded.cpp 149 2009-11-12 02:40:41Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-main/HOPSPACK_main_threaded.cpp $

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
  @file HOPSPACK_main_threaded.cpp
  @brief Main program that executes HOPSPACK as using multiple threads on
         a single machine.
*/

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <stdlib.h>
#endif

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_ExecutorMultiThreaded.hpp"
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


//----------------------------------------------------------------------
//  Internal Function allocateThreads_
//----------------------------------------------------------------------
/** Return false if user parameters are inconsistent.
 */
static bool  allocateThreads_ (const HOPSPACK::ParameterList &  cParams,
                                     int                     &  nNumThreads,
                                     int                     &  nNumCtznWrkrs)
{
    //---- HOPSPACK NEEDS 1 MAIN THREAD, AT LEAST 1 EVALUATION THREAD,
    //---- AND A (POSSIBLY EMPTY) POOL OF CITIZEN WORKER THREADS.
    //---- THE NUMBER OF THREADS ON THE OPERATING SYSTEM IS NOT CHECKED.

    HOPSPACK::ParameterList  cMedParams = cParams.sublist ("Mediator");

    //---- ACCEPT EITHER PARAMETER NAME.
    const char  sNUM_THREADS[] = "Number Threads";     //-- OFFICIAL NAME
    const char  sTHREADS[] = "Threads";                //-- ALTERNATE NAME

    if (   (cMedParams.isParameterInt (sTHREADS) == false)
        && (cMedParams.isParameterInt (sNUM_THREADS) == false) )
    {
        cerr << "ERROR: Need '" << sNUM_THREADS
             << "' in 'Mediator' sublist" << endl;
        return( false );
    }
    int  nNumReqThreads = -1;
    if (cMedParams.isParameterInt (sTHREADS) == true)
    {
        nNumReqThreads = cMedParams.getParameter (sTHREADS, -1);
    }
    else if (cMedParams.isParameterInt (sNUM_THREADS) == true)
    {
        nNumReqThreads = cMedParams.getParameter (sNUM_THREADS, -1);
    }
    if (nNumReqThreads < 2)
    {
        cerr << "ERROR: Bad '" << sNUM_THREADS << "' value " << nNumReqThreads
             << " in 'Mediator' sublist" << endl;
        cerr << "  Must be >= 2" << endl;
        if (nNumReqThreads == 1)
            cerr << "  For single-threaded operation run HOPSPACK_main_serial"
                 << endl;
        return( false );
    }
    #if defined(_WIN32)
        //---- WINDOWS THREAD LIMIT IS BASED ON STACK SIZE FOR EACH THREAD.
        //---- HENCE, THE LIMIT COULD BE CHANGED.
        if (nNumReqThreads > 2000)
        {
            cerr << "ERROR: Cannot request more than 2000 threads on Windows"
                 << endl;
            cerr << "  Please change '" << sNUM_THREADS
                 << "' in 'Mediator' sublist" << endl;
            return( false );
        }
    #endif

    const char  sRSRV_WRKRS[] = "Reserved Citizen Workers";
    nNumCtznWrkrs = 0;
    if (cMedParams.isParameterInt (sRSRV_WRKRS) == true)
    {
        nNumCtznWrkrs = cMedParams.getParameter (sRSRV_WRKRS, -1);
    }
    if (nNumCtznWrkrs < 0)
    {
        cerr << "ERROR: Bad '" << sRSRV_WRKRS << "' value " << nNumCtznWrkrs
             << " in 'Mediator' sublist" << endl;
        return( false );
    }
    if (nNumCtznWrkrs > (nNumReqThreads - 2))
    {
        cerr << "ERROR: Bad '" << sRSRV_WRKRS << "' value " << nNumCtznWrkrs
             << "' in 'Mediator' sublist" << endl;
        cerr << "  Cannot exceed '" << sNUM_THREADS << "' = " << nNumReqThreads
             << " - 2" << endl;
        return( false );
    }

    nNumThreads = nNumReqThreads;

    return( true );
}


//--------------------------------------------------------------------
//! Main routine for threaded version of HOPSPACK.
/*!
 *  The serial and multithreaded versions of HOPSPACK are very similar,
 *  differing primarily in the Executor component.  The multithreaded Executor
 *  creates a pool of Evaluators, each in a separate thread.
 *
 *  The main loop of work is performed by the Mediator in Hopspack::solve().
 *  It executes as the main thread of the process.
 *  Two pools of additional threads are created: one for Evaluators and one
 *  for asynchronous citizen workers.  Pool sizes are determined from user
 *  parameters, not the operating system.  Ideally, the user should configure
 *  HOPSPACK to set up one thread per processor/core on the machine.
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
    using HOPSPACK::ExecutorMultiThreaded;
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

    //---- ALLOCATE THREADS FOR VARIOUS TYPES OF WORKERS.
    int  nNumThreads = 0;
    int  nNumCtznWorkers = 0;
    if (allocateThreads_ (cParams, nNumThreads, nNumCtznWorkers) == false)
        return( nERROR );
    int  nNumEvalWorkers = nNumThreads - nNumCtznWorkers - 1;

    //---- CONSTRUCT A MULTITHREADED EXECUTOR.
    ExecutorMultiThreaded *  pExecutor = new ExecutorMultiThreaded();
    if (pExecutor->initialize (nNumEvalWorkers,
                               cParams.sublist ("Evaluator")) == false)
    {
        cerr << "ERROR: Could not construct Executor." << endl;
        return( nERROR );
    }

    //---- CONSTRUCT THE OPTIMIZER, CONFIGURE IT, AND RUN IT.
    Hopspack  optimizer (pExecutor);
    if (optimizer.setInputParameters (cParams) == true)
    {
        if (Print::doPrint (Print::FINAL_SOLUTION))
        {
            cout << endl << "Begin solving using the multithreaded"
                 << " HOPSPACK executable." << endl;
            if (Print::doPrint (Print::INPUT_PARAMETERS))
            {
                cout << "  1 thread for the Mediator" << endl;
                if (nNumCtznWorkers > 0)
                    cout << "  " << nNumCtznWorkers
                         << " threads for citizen workers" << endl;
                cout << "  " << nNumEvalWorkers
                     << " threads for evaluation workers" << endl;
            }
        }
        optimizer.solve();
    }

    if (pExecutor->shutdown() == false)
    {
        //---- ONE OR MORE WORKER THREADS ARE STUCK, SO KILL THE PROCESS.
        #if defined(_WIN32)
            cerr << "ERROR: Executor cannot shut down workers"
                 << " -- killing the process." << endl;
            ExitProcess( nERROR );
        #else
            exit (EXIT_FAILURE);
        #endif
    }
    delete pExecutor;

    return( nSUCCESS );
}
