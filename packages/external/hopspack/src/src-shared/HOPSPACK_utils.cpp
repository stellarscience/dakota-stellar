// $Id: HOPSPACK_utils.cpp 217 2013-11-25 21:59:49Z tplante $ 
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_utils.cpp $ 

//@HEADER
// ************************************************************************
// 
//         HOPSPACK: Hybrid Optimization Parallel Search Package
//                 Copyright 2009-2013 Sandia Corporation
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
  \file HOPSPACK_utils.cpp
  \brief Implement functions declared in HOPSPACK_utils.hpp.
*/

#include <stdlib.h>    //-- FOR system()
#include <sstream>

#include "HOPSPACK_utils.hpp"
#include "HOPSPACK_common.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_Matrix.hpp"

#ifndef SNL_TFLOPS_ENV


/*!  Get the next quoted string on the given line, starting at
  position pos. 

  \param line - Line of text from which to read

  \param pos - On input, the starting position in the line. On output,
  the next position after the quoted string (which may be
  std::string::npos). If there is any sort of error, this is set to
  std::string::npos upon return.

  \param value - On output, filled in with the quoted string (without
  the quotes). This is an empty string if no quoted string is found.

  \retval Returns true if the quoted string is successfully found,
  false otherwise.
  
*/
static bool getNextQuotedString(const string& line,
                                string::size_type& pos,
                                string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  string::size_type linelength = line.length();

  // Find the location of the first quote
  string::size_type pos1 = line.find('"', pos); 
  
  // Check that the operation succeeded and that we're not at the end of the line
  if ((pos1 == string::npos) || (pos1 == linelength - 1))
  {
    pos = string::npos;
    return false;
  }

  // Advance to first character after the first quote
  pos1 = pos1 + 1;

  // Find the location of the second quote
  string::size_type pos2 = line.find('"', pos1); 

  // Check that the operation was successful
  if (pos2 == string::npos)
  {
    pos = string::npos;
    return false;
  }

  // Compute the length of the expression
  string::size_type length = pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool HOPSPACK::getNextString(const string& line, string::size_type& pos, string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  string::size_type linelength = line.length();

  // Find the location of the first non-space character
  string::size_type pos1 = line.find_first_not_of(' ', pos); 
  
  // Check that the operation succeeded 
  if (pos1 == string::npos) 
  {
    pos = string::npos;
    return false;
  }

  // Find the location of the next space, if any
  string::size_type pos2 = line.find(' ', pos1); 

  // Compute the length of the expression
  string::size_type length = (pos2 == string::npos) ? linelength - pos1 : pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool HOPSPACK::getNextInt(const string& line, string::size_type& pos, int& value)
{
  string field;
  if ((!getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  //---- ACCEPT ONLY DIGITS OR +/-.  IN CONTRAST, sscanf ACCEPTS EVERYTHING
  //---- AND CAN PRODUCE A HUGE POSITIVE OR NEGATIVE INTEGER FOR A STRING INPUT.
  const char *  pCstr = field.c_str();
  for (int  i = 0; i < (int) field.size(); i++)
  {
      if (*pCstr == 0)
          break;
      if (   ((*pCstr < '0') || (*pCstr > '9'))
          && (*pCstr != '+')
          && (*pCstr != '-') )
      {
          value = -1;
          return( false );
      }
      pCstr++;
  }

  //---- USE sscanf TO DO THE PARSING.
  #if defined(HAVE_MSVC_SECURE_STRING_FNS)
      return ((sscanf_s(field.c_str(), "%d", &value)) == 1);
  #else
      return ((sscanf(field.c_str(), "%d", &value)) == 1);
  #endif
}

bool HOPSPACK::getNextDouble(const string& line, string::size_type& pos, double& value)
{
  string field;
  if ((!getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  #if defined(HAVE_MSVC_SECURE_STRING_FNS)
      return ((sscanf_s(field.c_str(), "%le", &value)) == 1);
  #else
      return ((sscanf(field.c_str(), "%le", &value)) == 1);
  #endif
}




/*! Splits line into a vector of strings based upon whitespace delimiter.
\par Example
The following code segment
\code
string line("one two three");
vector<string> linesplit;
tokenize(line,linesplit);
for (int i = 0; i < linesplit.size(); i++)
  cout << linesplit[i] << endl;
\endcode
will output
\verbatim
one
two
three
\endverbatim
*/
static void tokenize(const string& line, vector<string>& linesplit)
{
  linesplit.resize(0);
  string buf; // Have a buffer string
  stringstream ss(line); // Insert the string into a stream
  while (ss >> buf)
    linesplit.push_back(buf);
}


static bool  _bIsFirstErrMsg = true;
static string  _sParseFilename("");
static void  infileErrMsg(const string &  sErrMsg)
{
    if (_bIsFirstErrMsg)
    {
        cerr << "WARNINGS while parsing input file '"
             << _sParseFilename << "'" << endl;
        _bIsFirstErrMsg = false;
    }
    cerr << sErrMsg << endl;
    return;
}


bool HOPSPACK::processTextInputFileLine(const string& lineToProcess,
                                        ParameterList& params,
                                        ParameterList*& subPtr,
                                        ifstream &fin)
{
  string::size_type pos;        // current position inside line

  string field;                 // name of parameter

  int tmpint;                   // for reading in an int
  double tmpdouble;             // for reading in a double
  string tmpstring;             // for reading in a string
  Vector tmpvector;             // for reading in a HOPSPACK::Vector

  string type;                  // parameter type

  string line (lineToProcess);
  if (line.size() > 0)
  {
      //---- A PARAMETER FILE CREATED ON WINDOWS BUT READ ON UNIX
      //---- MAY STILL HAVE AN ASCII 13 AT THE END.
      if (line[line.size() - 1] == 13)
          line = line.substr (0, line.size() - 1);

      //---- REMOVE ANY LEADING WHITESPACE.
      size_t  firstChar = line.find_first_not_of (" \t");
      if (firstChar != string::npos)
          line = line.substr (firstChar, line.size());

      //---- REMOVE ANY CHARACTERS AFTER A MID-LINE COMMENT BEGINS.
      size_t  firstComment = line.find_first_of ("#");
      if ((firstComment != string::npos) && (firstComment > 0))
          line = line.substr (0, firstComment - 1);
  }

  if (line.size() == 0)         // empty line - which is okay
  {
    return true;
  }
  else if (line[0] == '#')      // comment - which is okay
  {
    return true;
  }
  else if (line[0] == '@')      // sublist command
  {
    subPtr = &params;           // reset to the top of the list (no nesting!)
    
    pos = 0;
    
    while (pos != string::npos)
    {
      // reset to sublist pointer if this is opening a new sublist
      if ((getNextQuotedString(line, pos, field)) && (field.size() > 0))
        subPtr = &(subPtr->getOrSetSublist(field));
    }

    return true;
  }
  else if (line[0] == '"')        // new style
  {

    // Get the name
    pos = 0;
    if ((!getNextQuotedString(line, pos, field)) || (field.empty()))
      return false;
    
    // Read in the type
    if (!getNextString(line, pos, type))
      return false;

    if (type == "int")
    {
      if (!getNextInt(line, pos, tmpint))
      {
        return false;
      }
      else
      {
        if (subPtr->isParameterInt(field))
        {
          infileErrMsg(" Parameter '" + field + "' already defined");
          return false;
        }
        subPtr->setParameter(field, tmpint);
        return true;
      }
    }
    else if (type == "bool")
    {
      if (!getNextString(line, pos, tmpstring))
      {
        return false;
      }
      else
      {
          if (subPtr->isParameterBool(field))
          {
              infileErrMsg(" Parameter '" + field + "' already defined");
              return false;
          }
          if (   (tmpstring == "true")
              || (tmpstring == "True")
              || (tmpstring == "TRUE")
              || (tmpstring == "T")
              || (tmpstring == "t") )
          {
              subPtr->setParameter (field, true);
              return true;
          }
          if (   (tmpstring == "false")
              || (tmpstring == "False")
              || (tmpstring == "FALSE")
              || (tmpstring == "F")
              || (tmpstring == "f") )
          {
              subPtr->setParameter (field, false);
              return true;
          }
          return false;
      }
    }
    else if (type == "double")
    {
      if (!getNextDouble(line, pos, tmpdouble))
      {
        return false;
      }
      else
      {
        if (subPtr->isParameterDouble(field))
        {
          infileErrMsg(" Parameter '" + field + "' already defined");
          return false;
        }
        subPtr->setParameter(field, tmpdouble);
        return true;
      }
    }
    else if (type == "string")
    {
      if ((!getNextQuotedString(line, pos, tmpstring)) || (tmpstring.empty()))
      {
        return false;
      }
      else
      {
        if (subPtr->isParameterString(field))
        {
          infileErrMsg(" Parameter '" + field + "' already defined");
          return false;
        }
        subPtr->setParameter(field, tmpstring);
        return true;
      }
    }
    else if (type == "charvec")
    {
      // get the size
      if (!getNextInt(line, pos, tmpint))
        return false;
      
      if (tmpint < 0)
        return false;

      vector< char >  tmpCharVec (tmpint);

      // Split remainder of line.
      string vecline;
      vecline = line.substr(pos);
      vector<string> linesplit;
      tokenize(vecline, linesplit);

      // Check that there is sufficient data present.
      if ((int) linesplit.size() != tmpint)
        return false;

      // Add in vector, checking for one character elements
      for (int i = 0; i < tmpint; i ++)
      {
        if ((int) linesplit[i].size() == 1)
          tmpCharVec[i] = linesplit[i][0];
        else
          return false;
      }
    
      if (subPtr->isParameterCharVec(field))
      {
        infileErrMsg(" Parameter '" + field + "' already defined");
        return false;
      }
      subPtr->setParameter(field, tmpCharVec);
      return true;
    }
    else if (type == "vector")
    {
      // get the size
      if (!getNextInt(line, pos, tmpint))
        return false;
      
      if (tmpint < 0)
        return false;
      
      tmpvector.resize(tmpint);

      // Split remainder of line.
      string vecline;
      vecline = line.substr(pos);
      vector<string> linesplit;
      tokenize(vecline, linesplit);

      // Check that there is sufficient data present.
      if ((int) linesplit.size() != tmpint)
        return false;

      // Add in vector, checking for infinity denoted by DNE.
      for (int i = 0; i < tmpint; i ++)
      {
        if (linesplit[i] == "DNE")
          tmpvector[i] = dne();
      #if defined(HAVE_MSVC_SECURE_STRING_FNS)
        else if (sscanf_s(linesplit[i].c_str(), "%le", &tmpvector[i]) != 1)
      #else
        else if (sscanf(linesplit[i].c_str(), "%le", &tmpvector[i]) != 1)
      #endif
          return false;
      }
    
      if (subPtr->isParameterVector(field))
      {
        infileErrMsg(" Parameter '" + field + "' already defined");
        return false;
      }
      subPtr->setParameter(field, tmpvector);        
      return true;
    }
    else if (type == "matrix")
    {
      int nrows;
      int ncols;

      // get number of rows
      if (!getNextInt(line, pos, nrows))
        return false;

      // get number of columns
      if (!getNextInt(line, pos, ncols))
        return false;
      
      if ( (nrows <= 0) || (ncols <= 0) )
        return false;

      Matrix mat;
      tmpvector.resize(ncols);
  
      // Grab data for matrix.
      string matline;
      vector<string> linesplit;
      for (int i = 0; i < nrows; i++)
      {
        // Process row i of matrix.
        if (fin.eof())
        {
          return false;
        }

        getline(fin, matline);

        // Split line.
        tokenize(matline, linesplit);          
        if ((int) linesplit.size() != ncols)
        {
          return false;
        }

        // Copy into matrix.
        for (int j = 0; j < ncols; j++)
        {
          if (linesplit[j] == "DNE")
            tmpvector[j] = dne();
        #if defined(HAVE_MSVC_SECURE_STRING_FNS)
          else if ((sscanf_s(linesplit[j].c_str(), "%le", &tmpvector[j])) != 1)
        #else
          else if ((sscanf(linesplit[j].c_str(), "%le", &tmpvector[j])) != 1)
        #endif
            return false;
        }
        mat.addRow(tmpvector);
      }      

      if (subPtr->isParameterMatrix(field))
      {
        infileErrMsg(" Parameter '" + field + "' already defined");
        return false;
      }
      subPtr->setParameter(field, mat);        
      return true;
    }
    else
    {
      return false;
    }

  } 
  else
  {
    return false;
  } 

  return false;
}


bool HOPSPACK::parseTextInputFile(const string filename, ParameterList& params)
{
  // Open the input file
  ifstream fin;
  fin.open(filename.c_str());

  if (!fin) 
  {
    cerr << "ERROR: Cannot find input file '" << filename << "'" << endl;
    cerr << "       Current working directory is " << system ("pwd") << endl;
    return false;
  }

  // Store the filename for any error messages.
  if (_sParseFilename.size() == 0)
  {
    _sParseFilename.assign(filename);
  }

  string line;                  // one line from input file
  ParameterList* subPtr;        // sublist pointer

  subPtr = &params;

  while (!fin.eof())
  {
    getline(fin, line);
    if (!processTextInputFileLine(line, params, subPtr,fin))
    {
      infileErrMsg(" Ignoring line: " + line);
    }
  }

  fin.close();

  return true;
}

#endif     //-- #ifndef SNL_TFLOPS_ENV
