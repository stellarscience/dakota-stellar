// $Id: HOPSPACK_LinConstr.cpp 166 2010-03-22 19:58:07Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-shared/HOPSPACK_LinConstr.cpp $

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
  @file HOPSPACK_LinConstr.cpp
  @brief Implementation of HOPSPACK::LinConstr.
*/

#include <iomanip>
#include <math.h>     //-- FOR fabs

#include "HOPSPACK_common.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_float.hpp"
#include "HOPSPACK_Matrix.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_ProblemDef.hpp"
#include "HOPSPACK_SolveLinConstrProj.hpp"
#include "HOPSPACK_Vector.hpp"


namespace HOPSPACK
{

LinConstr::LinConstr (const ProblemDef &  probDef)
    :
    probDef (probDef),
    scaling (probDef.getVarScaling())
{
    //---- SET DEFAULTS IN CASE initialize IS NOT CALLED.

    //---- THIS VALUE IS JUST SUFFICIENT FOR example2 ON ALL TEST MACHINES:
    //----    _dActiveTol = 4.0 * HOPSPACK::getMachineEpsilon();
    //---- HOWEVER, IT IS TOO TIGHT IF VARIABLES ARE SCALED TO THEIR LIMITS.
    _dActiveTol = 1.0e-12;

    _nDisplayFlag = 0;
    return;
}


LinConstr::LinConstr (const LinConstr &  cOther,
                      const bool         bDropEqs)
    :
    probDef (cOther.probDef),
    _dActiveTol (cOther._dActiveTol),
    _nDisplayFlag (cOther._nDisplayFlag),
    scaling (cOther.probDef.getVarScaling()),
    aIneq (cOther.aIneq),
    bIneqLower (cOther.bIneqLower),
    bIneqUpper (cOther.bIneqUpper)
{
    if (bDropEqs == false)
    {
        aEq = cOther.aEq;
        bEq = cOther.bEq;
    }

    if (setupScaledSystem() == false)
    {
        throwError ("constructor", "cannot set up scaled system");
    }

    return;
}


LinConstr::~LinConstr()
{
}

bool  LinConstr::initialize (const ParameterList &  cLinConstrParams)
{
  _dActiveTol = cLinConstrParams.getParameter ("Active Tolerance", _dActiveTol);

  _nDisplayFlag = cLinConstrParams.getParameter ("Display", _nDisplayFlag);
  if (_nDisplayFlag < 0)
    _nDisplayFlag = 0;
  if (_nDisplayFlag > 2)
    _nDisplayFlag = 2;

  if (setupMatrix (cLinConstrParams) == false)
      return( false );
  if (setupRhs (cLinConstrParams) == false)
      return( false );
  if (setupScaledSystem() == false)
      return( false );

  return( true );
}


// PRIVATE
bool LinConstr::setupMatrix(const ParameterList& params)
{
  if (params.isParameterMatrix("Inequality Matrix"))
  {
    #if !defined(HAVE_LAPACK)
      cerr << "ERROR: Cannot use linear constraints, no LAPACK in build" << endl;
      throw INTERNAL_ERROR;
    #endif

    aIneq = params.getMatrixParameter("Inequality Matrix");
    if ((aIneq.empty() == false) && (aIneq.getNcols() != scaling.size()))
    {
      cerr << "ERROR: Number of columns in 'Inequality Matrix' = "
           << aIneq.getNcols() << " does not match number variables = "
           << scaling.size() << endl;
      return( false );
    }

    for (int  i = 0; i < aIneq.getNrows(); i++)
    {
      const Vector  v = aIneq.getRow (i);
      for (int  j = 0; j < v.size(); j++)
      {
        if (exists (v[j]) == false)
        {
          cerr << "ERROR: DNE value is not allowed in 'Inequality Matrix'"
               << endl;
          return( false );
        }
      }
    }
  }

  if (params.isParameterMatrix("Equality Matrix"))
  {
    #if !defined(HAVE_LAPACK)
      cerr << "ERROR: Cannot use linear constraints, no LAPACK in build" << endl;
      throw INTERNAL_ERROR;
    #endif

    aEq = params.getMatrixParameter("Equality Matrix");
    if ((aEq.empty() == false) && (aEq.getNcols() != scaling.size()))
    {
      cerr << "ERROR: Number of columns in 'Equality Matrix' = "
           << aEq.getNcols() << " does not match number variables = "
           << scaling.size() << endl;
      return( false );
    }

    for (int  i = 0; i < aEq.getNrows(); i++)
    {
      const Vector  v = aEq.getRow (i);
      for (int  j = 0; j < v.size(); j++)
      {
        if (exists (v[j]) == false)
        {
          cerr << "ERROR: DNE value is not allowed in 'Equality Matrix'"
               << endl;
          return( false );
        }
      }
    }
  }

  return( true );
}

// PRIVATE
bool LinConstr::setupRhs(const ParameterList& params)
{  
  if (params.isParameterVector("Inequality Lower"))
    bIneqLower = params.getVectorParameter("Inequality Lower");
  else
    bIneqLower.assign( aIneq.getNrows(), dne());
  if (bIneqLower.size() != aIneq.getNrows())
  {
    cerr << "ERROR: Length of 'Inequality Lower' = " << bIneqLower.size()
         << " does not match 'Inequality Matrix' = " << aIneq.getNrows()
         << endl;
    return( false );
  }

  if (params.isParameterVector("Inequality Upper"))
    bIneqUpper = params.getVectorParameter("Inequality Upper");
  else
    bIneqUpper.assign(aIneq.getNrows(), dne());
  if (bIneqUpper.size() != aIneq.getNrows())
  {
    cerr << "ERROR: Length of 'Inequality Upper' = " << bIneqUpper.size()
         << " does not match 'Inequality Matrix' = " << aIneq.getNrows()
         << endl;
    return( false );
  }

  for (int i = 0; i < aIneq.getNrows(); i++)
  {
    if ((exists (bIneqLower[i]) == false) && (exists (bIneqUpper[i]) == false))
    {
      cerr << "ERROR: No bounds defined for inequality [" << (i + 1)
           << "] in sublist 'Linear Constraints'" << endl;
      return( false );
    }
    if (   (exists (bIneqLower[i])) && (exists (bIneqUpper[i]))
        && (bIneqLower[i] > bIneqUpper[i]) )
    {
      cerr << "ERROR: Bounds are inconsistent for inequality [" << (i + 1)
           << "] in sublist 'Linear Constraints'" << endl;
      return( false );
    }
  }

  if (params.isParameterVector("Equality Bounds"))
  {
    bEq = params.getVectorParameter("Equality Bounds");
    if (bEq.size() != aEq.getNrows())
    {
      cerr << "ERROR: Length of 'Equality Bounds' = " << bEq.size()
           << " does not match 'Equality Matrix' = " << aEq.getNrows()
           << endl;
      return( false );
    }
    for (int  i = 0; i < bEq.size(); i++)
    {
      if (exists (bEq[i]) == false)
      {
        cerr << "ERROR: No bound defined for equality [" << (i + 1)
             << "] in sublist 'Linear Constraints'" << endl;
        return( false );
      }
    }
  }
  else if (aEq.empty() == false)
  {
    cerr << "ERROR: Need 'Equality Bounds' to go with 'Equality Matrix'" << endl;
    return( false );
  }

  return( true );
}

// PRIVATE
bool LinConstr::setupScaledSystem()
{  
  // *** Form lHat ***
  const Vector &  bLower = probDef.getLowerBnds();
  for (int j = 0; j < scaling.size(); j++)
  {
    if (exists(bLower[j]))
      lHat.push_back(bLower[j]);
    else
      lHat.push_back(0);
  }

  // *** Form bHatLower and bHatUpper ***

  // Insert components corresponding to lower bounds
  for (int j = 0; j < scaling.size() ; j++)
  {
    if (exists(bLower[j]))
      bHatLower.push_back((bLower[j] - lHat[j]) / scaling[j]);
    else
      bHatLower.push_back(dne());
  }

  // Insert components corresponding to upper bounds
  const Vector &  bUpper = probDef.getUpperBnds();
  for (int j = 0; j < scaling.size(); j++)
  {
    if (exists(bUpper[j]))
      bHatUpper.push_back((bUpper[j] - lHat[j]) / scaling[j]);
    else
      bHatUpper.push_back(dne());
  }

  // Insert component corresponding to general linear constraints
  if (!aIneq.empty())
  {
    // Compute ail = aIneq * lhat
    Vector ail(aIneq.getNrows());
    aIneq.multVec(lHat, ail);

    
    for (int j = 0; j < aIneq.getNrows(); j++)
    {
      // Form scaled left-hand side.
      if (exists(bIneqLower[j]))
        bHatLower.push_back(bIneqLower[j] - ail[j]);
      else
        bHatLower.push_back(dne());

      // Form scaled right-hand side.
      if (exists(bIneqUpper[j]))
        bHatUpper.push_back(bIneqUpper[j] - ail[j]);
      else
        bHatUpper.push_back(dne());
    }
  }
  

  // *** Form aHat ***

  // First add identity
  aHat.setToIdentity(scaling.size());
  // Now add in aTildeIneq
  aHat.addMatrix(aIneq, scaling);
  
  // *** Form aTildeEq and bTildeEq ***
  // Nullspace matrix, Z, of aTildeEq is needed to compute distance to
  // inequality constraints within nullspace.
  Matrix ZT;
  if (!aEq.empty())
  {
    // ael = aEq * lhat
    Vector ael(aEq.getNrows());
    aEq.multVec(lHat,ael);

    for (int i = 0; i < aEq.getNrows(); i++)
      bTildeEq.push_back(bEq[i] - ael[i]);

    aTildeEq.scale(aEq, scaling);
    aTildeEq.nullSpace(ZT, _dActiveTol);
  }

  // Error checks.
  if (bHatLower.size() != aIneq.getNrows() + scaling.size())
  {
    cerr << "ERROR: Incorrect length for bHatLower  <LinConstr.initialize()>"
         << endl;
    return( false );
  }
  if (bHatUpper.size() != aIneq.getNrows() + scaling.size())
  {
    cerr << "ERROR: Incorrect length for bHatUpper  <LinConstr.initialize()>"
         << endl;
    return( false );
  }
  if (   (aHat.getNrows() != aIneq.getNrows() + scaling.size())
      || (aHat.getNcols() != scaling.size()) )
  {
    cerr << "ERROR: Incorrect length for aHat  <LinConstr.initialize()>"
         << endl;
    return( false );
  }

  // Store norms of rows for determining distances.
  aHatZNorm.resize(aHat.getNrows());
  // Project matrix into nullspace of scaled equality constraints.
  Matrix AZ = aHat;
  if (!ZT.empty())
    AZ.multMat(ZT, Matrix::TRANSPOSE);
  for (int i = 0; i < AZ.getNrows(); i++)
    aHatZNorm[i] = AZ.getRow(i).norm();

  return( true );
}


void LinConstr::scale(Vector& x) const
{
  if (x.size() != scaling.size())
    throwError("scale", "x vector has incorrect length");

  for (int i = 0; i < scaling.size(); i++)
    x[i] = (x[i] - lHat[i]) / scaling[i];
}

void LinConstr::unscale(Vector& w) const
{ 
  if (w.size() != scaling.size())
    throwError("scale", "w vector has incorrect length");

  for (int i = 0; i < scaling.size(); i++)
    w[i] = (scaling[i] * w[i]) + lHat[i];
}

double LinConstr::getActiveTol() const
{
  return _dActiveTol;
}

bool LinConstr::hasLinearConstraints() const
{
  if ((aIneq.getNrows() > 0) || (aEq.getNrows() > 0))
    return( true );
  return( false );
}

const Matrix & LinConstr::getAhat() const
{
  return aHat;
}

const Vector & LinConstr::getBhatLower() const
{
  return bHatLower;
}

const Vector & LinConstr::getBhatUpper() const
{
  return bHatUpper;
}

const Matrix & LinConstr::getAtildeEq() const
{
  return aTildeEq;
}

const Vector & LinConstr::getBtildeEq() const
{
  return bTildeEq;
}


bool LinConstr::isFeasible(const Vector& x,
                           const bool bPrintViolationInfo) const
{
  int n = scaling.size();
  // Error check
  if (x.size() != n)
    throwError("isFeasible", "x vector has incorrect length");

  // Create scaled x
  Vector xTilde(x);
  scale(xTilde);

  if (!isInequalityFeasible(xTilde, bPrintViolationInfo))
    return false;

  if (!isEqualityFeasible(xTilde, bPrintViolationInfo))
    return false;
  
  return true;
}

double  LinConstr::getL2Norm (const Vector &  x) const
{
    //---- COMPUTING WITH dnrm2 MAY BE MORE ACCURATE, BUT THEN
    //---- THE METHOD REQUIRES LAPACK.

    double  dResult = 0.0;

    //---- LOOP THRU VARIABLE BOUNDS.
    const Vector &  bLower = probDef.getLowerBnds();
    const Vector &  bUpper = probDef.getUpperBnds();
    for (int  i = 0; i < x.size(); i++)
    {
        if (exists (bLower[i]))
        {
            double  dViolation = bLower[i] - x[i];
            if (dViolation > dResult)
                dResult += dViolation * dViolation;
        }
        if (exists (bUpper[i]))
        {
            double  dViolation = x[i] - bUpper[i];
            if (dViolation > dResult)
                dResult += dViolation * dViolation;
        }
    }

    //---- LOOP THRU EQUALITY CONSTRAINTS.
    for (int  i = 0; i < aEq.getNrows(); i++)
    {
        const Vector &  nextCon = aEq.getRow (i);
        double  dConX = x.dot (nextCon);
        double  dViolation = fabs (dConX - bEq[i]);
        dResult += dViolation * dViolation;
    }

    //---- LOOP THRU INEQUALITY CONSTRAINTS.
    for (int  i = 0; i < aIneq.getNrows(); i++)
    {
        const Vector &  nextCon = aIneq.getRow(i);
        double  dConX = x.dot (nextCon);
        double  dViolation = 0.0;
        if (exists (bIneqLower[i]))
        {
            if (dConX < bIneqLower[i])
                dViolation = bIneqLower[i] - dConX;
        }
        if (exists (bIneqUpper[i]))
        {
            if (dConX > bIneqUpper[i])
                dViolation = dConX - bIneqUpper[i];
        }
        dResult += dViolation * dViolation;
    }

    return( sqrt (dResult) );
}

double  LinConstr::getLInfNorm (const Vector &  x) const
{
    double  dResult = 0.0;

    //---- LOOP THRU VARIABLE BOUNDS.
    const Vector &  bLower = probDef.getLowerBnds();
    const Vector &  bUpper = probDef.getUpperBnds();
    for (int  i = 0; i < x.size(); i++)
    {
        if (exists (bLower[i]))
        {
            double  dViolation = bLower[i] - x[i];
            if (dViolation > dResult)
                dResult = dViolation;
        }
        if (exists (bUpper[i]))
        {
            double  dViolation = x[i] - bUpper[i];
            if (dViolation > dResult)
                dResult = dViolation;
        }
    }

    //---- LOOP THRU EQUALITY CONSTRAINTS.
    for (int  i = 0; i < aEq.getNrows(); i++)
    {
        const Vector &  nextCon = aEq.getRow (i);
        double  dConX = x.dot (nextCon);
        double  dViolation = fabs (dConX - bEq[i]);
        if (dViolation > dResult)
            dResult = dViolation;
    }

    //---- LOOP THRU INEQUALITY CONSTRAINTS.
    for (int  i = 0; i < aIneq.getNrows(); i++)
    {
        const Vector &  nextCon = aIneq.getRow(i);
        double  dConX = x.dot (nextCon);
        double  dViolation = 0.0;
        if (exists (bIneqLower[i]))
        {
            dViolation = bIneqLower[i] - dConX;
            if (dViolation > dResult)
                dResult = dViolation;
        }
        if (exists (bIneqUpper[i]))
        {
            dViolation = dConX - bIneqUpper[i];
            if (dViolation > dResult)
                dResult = dViolation;
        }
    }

    return( dResult );
}


double  LinConstr::maxStep(const Vector& x,
                           const Vector& d,
                           double maxLength) const
{
  double maxStep = maxLength;
  int n = scaling.size();

  Vector xTilde(x);
  scale(xTilde);

  const Vector &  bLower = probDef.getLowerBnds();
  const Vector &  bUpper = probDef.getUpperBnds();

  // Determine max step for bounds component

  for (int i = 0; i < scaling.size(); i++)
  {
    if (d[i] < -scaling[i] * _dActiveTol) // Points to lower bound
      switch (getIneqState(i, LOWER_BOUND, xTilde))
      {
      case ACTIVE:      // Can't move if the constraint is active
        return 0;
      case VIOLATED:    // Or violated
        return 0;
      case INACTIVE:
        maxStep = min(maxStep, (bLower[i] - x[i]) / d[i]);
        break;
      case DNE:         // This means there is no lower bound
        break;
      }
    else if (d[i] > scaling[i] * _dActiveTol) // Points to upper bound
      switch (getIneqState(i, UPPER_BOUND, xTilde))
      {
      case ACTIVE:      // Can't move if the constraint is active
        return 0;
      case VIOLATED:    // Or violated
        return 0;
      case INACTIVE:
        maxStep = min(maxStep, (bUpper[i] - x[i]) / d[i]);
        break;
      case DNE:         // This means there is no upper bound
        break;
      }
  }

  // Determine max step for inequality component.

  if (!aIneq.empty())
  {
    int p = aIneq.getNrows();

    // aix = aIneq * x
    Vector aix(p);
    aIneq.multVec(x, aix);

    // aid = aIneq * d
    Vector aid(p);
    aIneq.multVec(d, aid);

    for (int j = 0; j < p; j++)
    {
      if (aid[j] < -_dActiveTol) // Points to lower bound
        switch (getIneqState(n + j, LOWER_BOUND, xTilde))
        {
        case ACTIVE:     // Can't move if the constraint is active
          return 0;
        case VIOLATED:   // Or violated
          return 0;
        case INACTIVE:
          maxStep = min(maxStep, (bIneqLower[j] - aix[j]) / aid[j]);
          break;
        case DNE:        // This means there is no lower bound
          break;
        }
      else if (aid[j] > _dActiveTol)  // Points to upper bound
        switch (getIneqState(n + j, UPPER_BOUND, xTilde))
        {
        case ACTIVE:     // Can't move if the constraint is active
          return 0;
        case VIOLATED:   // Or violated
          return 0;
        case INACTIVE:
          maxStep = min(maxStep, (bIneqUpper[j] - aix[j]) / aid[j]);
          break;
        case DNE:        // This means there is no upper bound
          break;
        }
    }
  }

  // Process equality constraints.

  if (!aEq.empty())
  {
    int p = aEq.getNrows();

    // z = aEq * d
    Vector z(p);
    aEq.multVec(d,z);

    // Check that direction stays on hyperplane defined by the
    // equality constraints
    for (int i = 0; i < p; i++)
      if (fabs(z[i]) > _dActiveTol)
        return 0;
  }

  return maxStep;
}

void LinConstr::getActiveIneqIndices(const Vector& xdist,
                                     double epsilon,
                                     vector<ActiveType>& index) const

{
  int nrows = aHat.getNrows();
  // Size the index
  index.resize(nrows);

  // Determine which bounds are near.
  for (int i = 0; i < nrows; i++)
  {
    // Check if upper or lower bounds are near.
    bool nearlow = (xdist[i] < epsilon);
    bool nearupp = (xdist[i+nrows] < epsilon);
    
    // Update index.
    if (nearlow && nearupp)
      index[i] = BOTH_ACTIVE;
    else if (nearlow)     // and not nearupp
      index[i] = LOWER_ACTIVE;
    else if (nearupp)     // and not nearlow
      index[i] = UPPER_ACTIVE;
    else
      index[i] = NEITHER_ACTIVE;
  }
}

void LinConstr::formDistanceVector(const Vector& x, Vector& xdist) const
{
  // Transform x to scaled space..
  Vector xTilde = x;
  scale(xTilde);
  
  int nrows = aHat.getNrows();

  // z = aHat*xTilde.
  Vector z(nrows);
  aHat.multVec(xTilde, z);
  
  // xdist must be twice the size of z to hold upper and lower bound info.
  xdist.resize(2*nrows);
  for (int i = 0; i < nrows; i++)
  {
    // Process lower bound.
    if (exists(bHatLower[i]))
    {
      if (aHatZNorm[i] > _dActiveTol)
        xdist[i] = fabs( z[i]  - bHatLower[i] ) / aHatZNorm[i];
      else if (fabs( z[i] - bHatLower[i] ) < _dActiveTol)
        xdist[i] = 0;
      else
        xdist[i] = dne();
    }
    else
      xdist[i] = dne();
    
    // Process upper bound.
    if (exists(bHatUpper[i]))
    {
      if (aHatZNorm[i] > _dActiveTol)
        xdist[i+nrows] = fabs( z[i]  - bHatUpper[i] ) / aHatZNorm[i];
      else if (fabs( z[i] - bHatUpper[i] ) < _dActiveTol)
        xdist[i+nrows] = 0;
      else
        xdist[i+nrows] = dne();
    }
    else
      xdist[i+nrows] = dne();
  }
}

void LinConstr::snapToBoundary(Vector& x,
                               double esnap) const
{
  // Form scaled x.
  Vector xTilde = x;
  scale(xTilde);

  // Form snap-to system.
  Matrix Asnap;
  Vector bsnap;
  formSnapSystem(xTilde, esnap, Asnap, bsnap);

  // Find closest point satisfying snap constraints.
  if (Asnap.specialConstrainedLSQR (xTilde, bsnap) == false)
      return;

  // Return the result.
  unscale(xTilde);
  x = xTilde;
  return;
}

// PRIVATE
void LinConstr::formSnapSystem(const Vector& xTilde, double esnap,
                               Matrix& Asnap,
                               Vector& bsnap) const
{
  Asnap.clear();
  bsnap.resize(0);
  
  if (aHat.empty())
  {
    if (!isEqualityFeasible(xTilde))
    {
      Asnap = aTildeEq;
      bsnap = bTildeEq;
    }
    return;
  }

  int nrow = aHat.getNrows();
  int nvar = aHat.getNcols();
  // Form z = Ahat*xTilde;
  Vector z(nrow);
  aHat.multVec(xTilde, z);
  
  // Measure and sort distance in terms of scaled space.
  multimap<double,int> dmap;
  multimap<double,int>::iterator iter;
  for (int i = 0; i < nrow; i++)
  { 
    if (exists(bHatLower[i]))
    { 
      double distlow = fabs( z[i] - bHatLower[i] ) / aHatZNorm[i];
      dmap.insert(std::pair<const double, int>(distlow, i));
    }
    if (exists(bHatUpper[i]))
    {
      // Indices for upper bounds in dmap denoted by i+nrow.
      // Symbolically stacking matrix.
      double distupp = fabs( z[i] - bHatUpper[i] ) / aHatZNorm[i];
      dmap.insert(std::pair<const double, int>(distupp, i+nrow));
    }
  }

  // Erase constraints which are outside the esnap ball.
  dmap.erase(dmap.upper_bound(esnap),dmap.end());

  if (   (dmap.lower_bound(_dActiveTol) == dmap.end())
      && isEqualityFeasible(xTilde) )
    return; // There are no constraints within distance esnap which are inactive.

  // Equality constraints are always added.
  Asnap = aTildeEq;
  bsnap = bTildeEq;

  // Now add inequality constraints within distance esnap.
  // Recall that we have symbolically stacked aHat, so we must mod out by nrow.
  for (iter=dmap.begin(); iter != dmap.end(); iter++)
  {
    if (bsnap.size() >= nvar)
      break;
    int irow = iter->second;
    if (irow < nrow)
    { // Snap to lower bound.
      Asnap.addRow(aHat.getRow(irow));
      bsnap.push_back(bHatLower[irow]);
    }
    else
    { // Snap to upper bound.
      Asnap.addRow(aHat.getRow(irow % nrow));
      bsnap.push_back(bHatUpper[irow % nrow]);
    }
  }

  // Now ensure that rank(Asnap) = Asnap.getNrows().
  Asnap.pruneDependentRows(bsnap,_dActiveTol);
}


bool  LinConstr::projectToFeasibility (Vector &  x) const
{
    SolveLinConstrProj  solver;
    Vector  xSolution;
    if (solver.solve (probDef, *this, x, xSolution) == false)
        return( false );

    x = xSolution;
    return( true );
}


void LinConstr::printDefinition (const bool  bDisplayFull) const
{
    if (_nDisplayFlag <= 0)
        return;

    if ((_nDisplayFlag < 2) || (bDisplayFull == false))
    {
        cout << "Linear Constraints" << endl;
    }
    else
    {
        cout << "Linear Constraints (full display)" << endl;
        cout << "  (Variable bounds are displayed in the Problem Definition)"
             << endl;
    }

    printCounts_();
    cout << "  Tolerance for feasibility, active constraints = "
         << setw(14) << setprecision(6)
         << setiosflags (ios::scientific) << getActiveTol() << endl;

    if ((_nDisplayFlag == 2) && (bDisplayFull == true))
    {
        //---- PRINT A LIST WITH EVERY INEQUALITY.
        if (aIneq.empty() == false)
        {
            cout << "  Inequality constraints:" << endl;
            for (int  i = 0; i < aIneq.getNrows(); i++)
            {
                cout << "    ";

                if (exists (bIneqLower[i]))
                {
                    cout << setw(14) << setprecision(6)
                         << setiosflags (ios::scientific) << bIneqLower[i];
                    cout << " <= ";
                }
                else
                    cout << "              " << "    ";

                printIneqName_ (i);

                if (exists (bIneqUpper[i]))
                {
                    cout << " <= ";
                    cout << setw(14) << setprecision(6)
                         << setiosflags (ios::scientific) << bIneqUpper[i];
                }
                else
                    cout << "    "  << "              ";

                cout << endl;
            }
            cout << "    A_ineq = ";
            aIneq.formattedPrint ("    ", cout);
            cout << endl;
        }

        //---- PRINT A LIST WITH EVERY EQUALITY.
        if (bEq.empty() == false)
        {
            cout << "  Equality constraints:" << endl;
            for (int  i = 0; i < bEq.size(); i++)
            {
                cout << "    ";
                printEqName_ (i);
                cout << " = " << setw(14) << setprecision(6)
                     << setiosflags (ios::scientific) << bEq[i] << endl;
            }
            cout << "    A_eq = ";
            aEq.formattedPrint ("    ", cout);
            cout << endl;
        }

        cout << "End of Linear Constraints (full display)" << endl;
    }

    cout << endl;
    return;
}

// PRIVATE
void LinConstr::printCounts_ (void) const
{
    int  nNumIneqLower = 0;
    int  nNumIneqUpper = 0;
    for (int i = 0; i < bIneqLower.size(); i++)
    {
        if (exists (bIneqLower[i]))
            nNumIneqLower++;
        if (exists(bIneqUpper[i]))
            nNumIneqUpper++;
    }

    cout << "  Constraint count summary:" << endl;
    cout << "  " << setw (5) << scaling.size() << " variables" << endl;
    cout << "  " << setw (5) << (nNumIneqLower + nNumIneqUpper)
         << " inequality constraints" << endl;
    cout << "  " << setw (5) << bEq.size() << " equality constraints" << endl;
    return;
}

// PRIVATE
void LinConstr::printIneqName_ (const int  nIneqNum) const
{
    cout << "c_ineq[" << setw (3) << nIneqNum << "]";
    return;
}

// PRIVATE
void LinConstr::printEqName_ (const int  nEqNum) const
{
    cout << "c_eq[" << setw (3) << nEqNum << "]";
    return;
}

// PRIVATE
void LinConstr::throwError(const string& fname, const string& msg) const
{
  cerr << "ERROR: " << msg << "  <" << fname << ">" << endl;
  throw INTERNAL_ERROR;
}

// PRIVATE
LinConstr::StateType  LinConstr::getIneqState
                          (const int i,
                           const BoundType bType,
                           const Vector& xTilde,
                           const bool bPrintViolationInfo) const
{
  // a = constraint row
  const Vector& a = aHat.getRow(i);
  double anorm = aHatZNorm[i];
  // b = lhs
  double b = (bType == LOWER_BOUND) ? bHatLower[i] : bHatUpper[i];

  // Check finiteness of the bound
  if (!exists(b))
    return DNE;

  // z = a' * xTilde
  double z;
  z = xTilde.dot(a);
  double  xnorm = xTilde.norm();

  // Check if the constraint is epsilon-active
  if ( fabs(z - b) < (_dActiveTol * max (anorm, xnorm)) )
    return ACTIVE;

  // Check if the constraint is otherwise satisfied
  if (((bType == LOWER_BOUND) && ( b <= z )) ||
      ((bType == UPPER_BOUND) && ( z <= b )))
    return INACTIVE;

  // Otherwise, it must be violated.
  if (bPrintViolationInfo)
  {
      cout << "     Inequality[" << i << "] violated by " << fabs(z - b)
           << " (tolerance = " << (_dActiveTol * max (anorm, xnorm)) << ")"
           << endl;
  }
  return VIOLATED;
}

// PRIVATE
LinConstr::StateType  LinConstr::getEqState
                          (const int i,
                           const Vector& xTilde,
                           const bool bPrintViolationInfo) const
{
  // a = constraint row
  const Vector& a = aTildeEq.getRow(i);
  double anorm = a.norm();

  // b = rhs
  double b = bTildeEq[i];

  // z = a' * xTilde
  double z;
  z = xTilde.dot(a);
  double  xnorm = xTilde.norm();
  // Check if the constraint is epsilon-active
  if ( fabs(z - b) < (_dActiveTol * max (anorm, xnorm)) )
    return ACTIVE;

  // Otherwise, it must be violated.
  if (bPrintViolationInfo)
  {
      cout << "     Equality[" << i << "] violated by " << fabs(z - b)
           << " (tolerance = " << (_dActiveTol * max (anorm, xnorm)) << ")"
           << endl;
  }
  return VIOLATED;
}

// PRIVATE
bool LinConstr::isEqualityFeasible(const Vector& xTilde,
                                   const bool    bPrintViolationInfo) const
{  
  // Check equality constraints
  for (int i = 0; i < aTildeEq.getNrows(); i ++)
  {
    if (VIOLATED == getEqState(i, xTilde, bPrintViolationInfo))
      return false;
  }
  
  return true;
}

// PRIVATE
bool LinConstr::isInequalityFeasible(const Vector& xTilde,
                                     const bool    bPrintViolationInfo) const
{  
  // Check inequality constraints
  for (int i = 0; i < aHat.getNrows(); i ++)
  {
    if (VIOLATED == getIneqState(i, UPPER_BOUND, xTilde, bPrintViolationInfo))
      return false;
    if (VIOLATED == getIneqState(i, LOWER_BOUND, xTilde, bPrintViolationInfo))
      return false;
  }
  
  return true;
}

}     //-- namespace HOPSPACK
