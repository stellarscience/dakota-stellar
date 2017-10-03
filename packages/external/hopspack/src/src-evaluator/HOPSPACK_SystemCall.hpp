// $Id: HOPSPACK_SystemCall.hpp 220 2014-01-02 21:24:59Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-evaluator/HOPSPACK_SystemCall.hpp $

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
  @file HOPSPACK_SystemCall.hpp
  @brief Declaration for HOPSPACK::SystemCall, implements Evaluator.
*/

#ifndef HOPSPACK_SYSTEMCALL_HPP
#define HOPSPACK_SYSTEMCALL_HPP

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_Evaluator.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Vector.hpp"

namespace HOPSPACK
{


//----------------------------------------------------------------------
//! Implements Evaluator to get information by calling an application executable.
/*!
 *  The implementation must be compatible with serial, MPI, and multithreaded
 *  versions of HOPSPACK.
 *
 *  System call finds the name of a user-supplied executable from user parameters.
 *  A single point is evaluated by writing arguments to a uniquely named file,
 *  making a system call to run the executable as a separate process, and
 *  reading results from another uniquely named file.  The two files are
 *  deleted after results are read and stored by HOPSPACK.
 *
 *  File names are of the form "user_specified_prefix.NNN_TT", where NNN
 *  is a unique integer tag, and TT is the type of information requested.
 *  This information is also passed to the evaluation executable in case
 *  it needs to create uniquely named temporary files.
 *
 *  The input file begins with the SS request code on the first line, and then
 *  the point at which to compute the evaluation.  The point, like all vectors,
 *  is formatted with each component on a separate line, as a floating point
 *  number in exponential notation, with the number of significant digits
 *  determined by a configuration parameter.
 *
 *  The output file, if evaluation was successful, begins with the number
 *  of objectives returned and then their values.  This is followed by the
 *  number of nonlinear equality constraints and their values, and then the
 *  same for nonlinear inequalities.  If there are no nonlinear constraints,
 *  then both parts can be skipped.  If there are only one kind of nonlinear
 *  constraint, the missing type should write a zero in the output file.
 *
 *  Linear constraints need not be evaluated and should not be reported
 *  (citizens already know about them from configuration parameters).
 *  If evaluation was unsuccessful, then the output file should contain a
 *  simple text error message that is passed on.
 */
//----------------------------------------------------------------------
class SystemCall : Evaluator
{
  public:

    //! Constructor.
    /*!
     *  @param[in] cEvalParams    Parameters in the "Evaluator" sublist.
     *  @throws string if parameter error is detected.
     */
    SystemCall (const ParameterList &  cEvalParams);

    //! Destructor.
    ~SystemCall (void);

    //! Evaluate the objective function(s) at a point x.
    /*!
     *  @param[in] nTag   Contains a unique tag for the evaluation which can be
     *                    used to name files, etc.
     *  @param[in] cX     The point at which to evaluate the function(s).
     *  @param[out] cFns  On output, contains a vector of objective function
     *                    values computed at X.  Multiple objectives are allowed.
     *                    If an evaluation failed, return an empty vector or set
     *                    individual elements of the vector to HOPSPACK::dne().
     *  @param[out] sMsg  On output, contains a message about the evaluation;
     *                    typically the word "Success" or an error message.
     */
    void  evalF (const int       nTag,
                 const Vector &  cX,
                       Vector &  cFns,
                       string &  sMsg);

    //! Evaluate the objective functions and nonlinear constraints at a point x.
    /*!
     *  @param[in] nTag     Contains a unique tag for the evaluation which can be
     *                      used to name files, etc.
     *  @param[in] cX       The point at which to evaluate the function(s).
     *  @param[out] cFns    On output, contains a vector of objective function
     *                      values computed at X.  Multiple objectives are
     *                      allowed.  If an evaluation failed, return an empty
     *                      vector or set individual function elements of the
     *                      vector to HOPSPACK::dne().
     *  @param[out] cEqs    On output, contains a vector of nonlinear equality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] cIneqs  On output, contains a vector of nonlinear inequality
     *                      constraint function values computed at X.  If an
     *                      evaluation failed, return an empty vector or set
     *                      individual elements of the vector to HOPSPACK::dne().
     *  @param[out] sMsg    On output, contains a message about the evaluation;
     *                      typically the word "Success" or an error message.
     */
    void  evalFC (const int       nTag,
                  const Vector &  cX,
                        Vector &  cFns,
                        Vector &  cEqs,
                        Vector &  cIneqs,
                        string &  sMsg);

    //! Return the string from parameter 'Evaluator Type'.
    string  getEvaluatorType (void) const;

    //! Print debug information about the Evaluator instance.
    void  printDebugInfo (void) const;

  private:

    //! By design, there is no copy constructor
    SystemCall (const SystemCall &);
    //! By design, there is no assignment operator.
    SystemCall & operator= (const SystemCall &);

    //! Generate input and output file names, and the executable system call.
    /*!
     *  @param[in]  nTag             Tag number for file name suffix.
     *  @param[in]  sReqType         Request type code file file name suffix.
     *  @param[out] sInputFileName   Filled with the input file name.
     *  @param[out] sOutputFileName  Filled with the output file name.
     *  @param[out] sSysCall         Filled with the system call.
     */
    void  generateStrings_ (const int       nTag,
                            const string &  sReqType,
                                  string &  sInputFileName,
                                  string &  sOutputFileName,
                                  string &  sSysCall) const;

    //! Write the input file.
    bool  writeInputFile_ (const string &  sInputFileName,
                           const string &  sReqType,
                           const Vector &  cX) const;

    //! Read a vector from the output file.
    /*!
     *  @param[in,out] fptr       File pointer positioned at start of vector.
     *  @param[in]     sFileName  Name of the file being read.
     *  @param[out]    cV         Filled with vector values.
     *  @param[out]    sMsg       Filled if no vector is found.
     *  @return true              If successful; false means cV is undefined
     *                            and sMsg contains any evaluation message.
     */
    bool  readVector_ (      ifstream &  fptr,
                       const string   &  sFileName,
                             Vector   &  cV,
                             string   &  sMsg) const;

    //! Delete an I/O file.
    void  deleteIOFile_ (const string &  sFileName) const;


    string  _sExecutableName;
    string  _sInputPrefix;
    string  _sOutputPrefix;
    int     _nPrecisionDigits;     //-- SIGNIFICANT DIGITS TO WRITE OUT
    bool    _bSaveIOFiles;         //-- TRUE MEANS SAVE ALL TEMP I/O FILES
    int     _bDebug;

};

}          //-- namespace HOPSPACK

#endif     //-- HOPSPACK_SYSTEMCALL_HPP
