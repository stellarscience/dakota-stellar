// $Id: HOPSPACK_CacheManager.cpp 149 2009-11-12 02:40:41Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-framework/HOPSPACK_CacheManager.cpp $ 

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

/*! \file HOPSPACK_CacheManager.cpp
    \brief Implement HOPSPACK::CacheManager
*/

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_CacheManager.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_utils.hpp"

HOPSPACK::CacheManager::CacheManager (const ParameterList  &  params) :
    isFout(false)
{
  treeptr = new CacheSplayTree<CachePoint>;

  precision = params.getParameter("Cache Output Precision", 14);
  if (precision < 0)
  {
      cerr << "WARNING: Illegal 'Cache Output Precision' value"
           << " in 'Mediator' sublist" << endl;
      cerr << "         Changing 'Cache Output Precision' to zero" << endl;
      precision = 0;
  }

  inname = params.getParameter("Cache Input File", inname);
  parseInputFile(inname);

  outname = params.getParameter("Cache Output File", outname);
  openOutputFile(outname);
} 

HOPSPACK::CacheManager::~CacheManager()
{
  delete treeptr;
  closeOutputFile();
}

bool HOPSPACK::CacheManager::insert(const Vector& x,
                                    const Vector& f,
                                    const Vector& cEqs,
                                    const Vector& cIneqs)
{
  CachePoint cp(x,f,cEqs,cIneqs);
  bool isInserted = treeptr->insert(cp);
  if (isInserted)
    writeToOutputFile(x,f,cEqs,cIneqs);
  return isInserted;
}

bool HOPSPACK::CacheManager::isCached(const Vector& x,
                                            Vector& f,
                                            Vector& cEqs,
                                            Vector& cIneqs)
{
  CachePoint cp(x,f,cEqs,cIneqs);

  if ( ! (treeptr->find(cp)) )
  {
    return false;
  }

  f = cp.getF();
  cEqs = cp.getEqs();
  cIneqs = cp.getIneqs();
  return true;
}

void  HOPSPACK::CacheManager::printDebugInfo (void) const
{
    cout << "  HOPSPACK_Cache" << endl;

    cout << "    Cache Input File:       " << inname;
    if (bCanOpenInname == false)
        cout << "  (could not open file)";
    cout << endl;

    cout << "    Cache Output File:      " << outname << endl;
    cout << "    Cache Output Precision: " << precision << endl;
    cout << "    current num points in cache = "
         << treeptr->getNumNodes() << endl;

    return;
}


// PRIVATE
void HOPSPACK::CacheManager::parseInputFile(const string &  inname)
{
  bCanOpenInname = true;

  if (inname.empty())
    return;

  ifstream fin;
  fin.open(inname.c_str());
  
  if (!fin)
  {
    cerr << "WARNING: Cannot open cache input file '" << inname << "'"
         << endl;
    bCanOpenInname = false;
    return;
  }

  string line;
  while (!fin.eof())
  {
    getline(fin, line);
    if (processInputLine(line) == false)
    {
       cerr << "WARNING: Error parsing cache input line, point is ignored"
            << endl;
    }
  }

  fin.close();
}

// PRIVATE
bool HOPSPACK::CacheManager::processInputLine(string& line)
{

  string::size_type line_pos;
  
  string element;
  

  line_pos = 0;

  // Empty line is OK
  if (!getNextString(line, line_pos, element))
    return( true );

  // Stop reading if line does not have the right format
  if (element != "x=[") 
    return( false );

  // Read elements of X
  Vector x;
  if (readVectorFromLine (line, line_pos, x) == false)
    return( false );

  // X cannot be empty.
  if (x.size() == 0)
    return( false );


  // Read next token
  if (!getNextString(line, line_pos, element))
    return( false );

  // Stop reading if line does not have the right format
  if (element != "f=[") 
    return( false );

  // Read elements of F
  Vector f;
  if (readVectorFromLine (line, line_pos, f) == false)
    return( false );


  // Read next token
  if (!getNextString(line, line_pos, element))
    return( false );

  // Stop reading if line does not have the right format
  if (element != "c_e=[") 
    return( false );

  // Read elements of cEqs
  Vector cEqs;
  if (readVectorFromLine (line, line_pos, cEqs) == false)
    return( false );


  // Read next token
  if (!getNextString(line, line_pos, element))
    return( false );

  // Stop reading if line does not have the right format
  if (element != "c_i=[") 
    return( false );

  // Read elements of cEqs
  Vector cIneqs;
  if (readVectorFromLine (line, line_pos, cIneqs) == false)
    return( false );


  // Insert into the queue
  insert (x, f, cEqs, cIneqs);

  return( true );
}

// PRIVATE
bool  HOPSPACK::CacheManager::readVectorFromLine (string            &  line,
                                                  string::size_type &  line_pos,
                                                  Vector            &  result)
{
  result.resize (0);

  string nextElement;
  while (1)
  {
    // Read the next element
    if (!getNextString (line, line_pos, nextElement))
      return( false );

    // Check if the vector is empty
    if (nextElement == "(empty)")
    {
      if (!getNextString (line, line_pos, nextElement))
        return( false );
      if (nextElement != "]")
        return( false );
      return( true );
    }

    // Check for termination
    if (nextElement == "]")
      return( true );

    // Process into a double value
    string::size_type element_pos = 0;
    double  d;
    if (!getNextDouble (nextElement, element_pos, d))
      return( false );
    result.push_back (d);
  }
}

// PRIVATE
void HOPSPACK::CacheManager::openOutputFile(const string &  filename)
{
  if (filename.empty())
    return;

  // Open file for append
  fout.open(filename.c_str(), std::ios::out|std::ios::app);
  
  if (!fout)
  {
    cerr << "WARNING: Cannot open cache output file '" << filename << "'"
         << endl;
    return;
  }

  isFout = true;
}

// PRIVATE
void HOPSPACK::CacheManager::writeToOutputFile(const Vector& x,
                                               const Vector& f,
                                               const Vector& cEqs,
                                               const Vector& cIneqs)
{
  if (!isFout)
    return;

  fout << "x=[ ";
  x.leftshift (fout,precision);
  fout << " ]";

  fout << " f=[ ";
  f.leftshift (fout,precision);
  fout << " ]";

  fout << " c_e=[ ";
  cEqs.leftshift (fout,precision);
  fout << " ]";

  fout << " c_i=[ ";
  cIneqs.leftshift (fout,precision);
  fout << " ]";

  fout << endl;
  fout.flush();
}

// PRIVATE
void HOPSPACK::CacheManager::closeOutputFile()
{
  if (!isFout)
    fout.close();
}
