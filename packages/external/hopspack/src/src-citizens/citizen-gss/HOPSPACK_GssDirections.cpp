// $Id: HOPSPACK_GssDirections.cpp 216 2013-11-13 23:34:51Z tplante $
// $URL: https://software.sandia.gov/svn/hopspack/trunk/src/src-citizens/citizen-gss/HOPSPACK_GssDirections.cpp $

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
  \file HOPSPACK_GssDirections.cpp
  \brief Implemtation of HOPSPACK::GssDirections.
*/

#include <stdlib.h>    //-- FOR free()
#include <iomanip>

#include "HOPSPACK_common.hpp"
#ifdef HAVE_CDDLIB
  #include "HOPSPACK_CddLibWrapper.h"
#endif
#include "HOPSPACK_GssDirections.hpp"
#include "HOPSPACK_LinConstr.hpp"
#include "HOPSPACK_ParameterList.hpp"
#include "HOPSPACK_Print.hpp"
#include "HOPSPACK_ProblemDef.hpp"

namespace HOPSPACK
{

GssDirections::GssDirections (const ProblemDef    &  probDef,
                              const LinConstr     &  linConstr,
                                    ParameterList &  params) :
  probDef(probDef),
  linConstr(linConstr),
  nDimensions(probDef.getVarScaling().size()),
  zero(nDimensions, static_cast<double>(0)),
  nDirections(0),
  nCached(0),
  nLapack(0),
  nCddLib(0),
  nMaxDirections(0),
  nAppend(0)
{
  stepTolerance = params.getOrSetParameter("Step Tolerance", 0.01);
  minStep = params.getOrSetParameter("Minimum Step", 2 * stepTolerance);
  theta = params.getOrSetParameter("Contraction Factor", 0.5);
  epsilonMax = params.getOrSetParameter("Epsilon Max", 2 * stepTolerance);
  withNormals = params.getOrSetParameter("Add Projected Normals", true);
  withCompass = params.getOrSetParameter("Add Projected Compass", false);

  // Check parameters
  if (stepTolerance <= 0)
  {
    cerr << "ERROR: 'Step Tolerance' must be positive  <GssDirections>." << endl;
    throw "GSS Error";
  }

  if (minStep <= stepTolerance)
  {
    cerr << "ERROR: 'Minimum Step' must be greater than 'Step Tolerance'"
         << "  <GssDirections>." << endl;
    throw "GSS Error";
  }

  if ((theta <= 0) || (theta >= 1))
  {
    cerr << "ERROR: 'Contraction Factor' must be strictly between zero and one"
         << "  <GssDirections>." << endl;
    throw "GSS Error";
  }

  epsilonMin = epsilonMax;
}

GssDirections::~GssDirections()
{
}

const Vector& GssDirections::getDirection(int i) const
{
  return directionsMatrix.getRow(i);
}

double GssDirections::getStep(int i) const
{
  return step[i];
}

void GssDirections::getDirectionIndices(vector<int>& idvector) const
{
  idvector.resize(0);
  for (int i = 0; i < nDirections; i ++)
    if ((step[i] >= stepTolerance) && (tag[i] == -1))
      idvector.push_back(i);
}

void GssDirections::computeNewDirections(const GssPoint& newPoint)
{
  // Grab new point info.
  const Vector& x = newPoint.getX();
  
  // Update distance vector for new point.
  linConstr.formDistanceVector( x, xDistance );

  double newStep = max(newPoint.getStepLength(), minStep);

  do 
  {
    // Update active state vector, member constraintState.
    if (updateConstraintState(newStep))
    {    
      // Generate search directions.
      generateForLinear(directionsMatrix);
    }

    if (!directionsMatrix.empty())
    {
      // Update step, trueStep, epsilonMin, and nDirections.
      updateDirectionInfo(newStep);
      break;
    }

    // If directions are empty, we will want a smaller value of newStep.
    newStep = theta * newStep;
  } while (newStep >= stepTolerance);
  
  if (directionsMatrix.empty())
  {
    cerr << "ERROR: Cannot compute generators for epsilon-tangent cone" << endl
         << "       <GssDirections::computeNewDirections()>." << endl
         << "       Most likely the problem is one of the following:" << endl
         << "       (1) Parameter 'Step Tolerance' is too large" << endl
         << "       (2) No feasible search directions exist at the current point."
         << endl;
    throw "GSS Error";
  }
}

void GssDirections::appendNewDirections()
{
  double newEps = getSmallestStep();
  // If newEps >= epsilonMin active set remains the same.
  if (newEps >= epsilonMin)
    return;

  // If constraintState remains the same, do nothing.
  if (!updateConstraintState(newEps))
  {
    epsilonMin = newEps;
    return;
  }

  // Compute new directions.  Note that distance vector has not been updated
  // as the current point is assumed the same.
  Matrix newDirectionsMatrix;
  generateForLinear(newDirectionsMatrix);

  // Now add new directions to the current list.
  directionsMatrix.addUniqueRows(newDirectionsMatrix, 1.0e-12);
  directionsMatrix.addUniqueRows(newDirectionsMatrix, linConstr.getActiveTol());
  
  // Update step, trueStep, espilonk, and nDirections.
  updateDirectionInfo(newEps, true);
}

void GssDirections::setStepConverged(int i)
{
  step[i] = stepTolerance/2;
}

void GssDirections::setTrueStepAndTag(int i, double trueStep_in, int tag_in)
{
  trueStep[i] = trueStep_in;
  tag[i] = tag_in;
}


void GssDirections::print(const string label) const
{
  if (!label.empty())
    cout << label << ":" << endl;

  int  nPrec = Print::getPrecision();

  for (int i = 0; i < nDirections; i ++)
  {
    cout << setw(4) << i << " : ";

    cout << "d =[";
    directionsMatrix.getRow(i).leftshift (cout);
    cout << "] ";

    cout.setf(ios::scientific);
    cout.precision(nPrec);
    cout << "step = " << setw(nPrec + 7) << step[i] << " ";
    cout.unsetf(ios::scientific);

    if (tag[i] != -1) 
    {
      cout << "tag = " << setw(4) << tag[i] << " ";

      cout.setf(ios::scientific);
      cout.precision(nPrec);
      cout << "trueStep = " << setw(nPrec + 7) << trueStep[i];
      cout.unsetf(ios::scientific);
    }

    cout << endl;
  }
  cout << " Number of times directions calculated by..." << endl
       << "  LAPACK: "    << nLapack << endl
       << "  CDDLIB: "    << nCddLib << endl
       << "  Cached: " << nCached  << endl
       << " Max directions in single iteration : " << nMaxDirections << endl
       << " Number of times directions appended: " << nAppend << endl;
}

bool GssDirections::isStepConverged() const
{
  for (int i = 0; i < nDirections; i ++)
  {
    if (step[i] >= stepTolerance)
      return false;
  }

  return true;
}

bool GssDirections::empty() const
{
  return (directionsMatrix.empty() == true);
}

void GssDirections::reduceStep(int i)
{
  double tmpStep = theta * step[i];

  step[i] = tmpStep;
  trueStep[i] = -1;
  tag[i] = -1;
}

//PRIVATE
void GssDirections::buildNormalCone(Matrix& VpT,
                                    Matrix& VlT) const
{
  // Always insert equality constraints in VlT.
  VlT.addMatrix(linConstr.getAtildeEq());
  
  // Insert elements into VpT and VlT from aHat.
  const Matrix& aHat = linConstr.getAhat();
  for (int i = 0; i < (int) constraintState.size(); i++)
  {
    if (constraintState[i] == LinConstr::BOTH_ACTIVE)
      VlT.addRow(aHat.getRow(i));
    else if (constraintState[i] == LinConstr::LOWER_ACTIVE)
      VpT.addRow(aHat.getRow(i), -1.0);
    else if (constraintState[i] == LinConstr::UPPER_ACTIVE)
      VpT.addRow(aHat.getRow(i));
  }
}

//PRIVATE
void GssDirections::buildTangentCone(const Matrix& VpT,
                                     const Matrix& VlT,
                                     Matrix& T)
{
#if !defined(HAVE_LAPACK)
  buildWithNothing(T);
  return;
#endif
  
  if ((VpT.empty()) && (VlT.empty()))
  {
    // Problem looks locally unconstrained. Generate all compass directions.
    generateUnconstrained(T);
    return;
  }
  
  if (buildWithLapack(VpT, VlT, T))
  {
    nLapack++;
    return;
  }
  
#ifdef HAVE_CDDLIB
  // Either constraints degenerate or LAPACK unavailable.  Try CDDLIB.
  if (buildWithCddLib(VpT, VlT, T))
  {
    nCddLib++;
    return;
  }
#endif  

  cerr << "ERROR: Cannot compute generators for epsilon-tangent cone" << endl
       << "       <GssDirections::buildTangentCone()>." << endl
       << "       Most likely the problem has degenerate constraints, and CDDLIB is" << endl
       << "       (1) not configured with HOPSPACK, or" << endl
       << "       (2) failed to compute generators (please send a bug report!)"
       << endl;
  throw "GSS Error";
}

//PRIVATE
void GssDirections::buildWithNothing(Matrix& D)
{
  D.clear();
  const Vector& scaling = probDef.getVarScaling();
  for (int i = 0; i < (int) constraintState.size(); i++)
  {
    tmpVector = zero;
    if (constraintState[i] == LinConstr::NEITHER_ACTIVE)
    {
      // Add +e_i
      tmpVector[i] = scaling[i];
      D.addRow(tmpVector);
      // Add -e_i
      tmpVector[i] = -1 * scaling[i];
      D.addRow(tmpVector);
    }
    else if (constraintState[i] == LinConstr::LOWER_ACTIVE)
    {
      // Add +e_i
      tmpVector[i] = scaling[i];
      D.addRow(tmpVector);
    }
    else if (constraintState[i] == LinConstr::UPPER_ACTIVE)
    {
      // Add -e_i
      tmpVector[i] = -1 * scaling[i];
      D.addRow(tmpVector);
    }
  }
}

#ifdef HAVE_CDDLIB
//PRIVATE
bool GssDirections::buildWithCddLib( const Matrix& VpT,
                                     const Matrix& VlT,
                                     Matrix& Tcone)
{
  Tcone.clear();  
  const Vector& scaling = probDef.getVarScaling();
  int nvar = scaling.size();

  // Get row pointers needed by call to compute_cone_generators.
  vector< double *> VpTptr;
  const_cast<Matrix &>(VpT).getModifiableRowPointers(VpTptr);

  vector< double *> VlTptr;
  const_cast<Matrix &>(VlT).getModifiableRowPointers(VlTptr);
  
  int num_pointy = 0;
  int num_lineality = 0;
  double **P;
  double **L;

  // If vectors are empty, cannot access zeroth element.
  double **Eq;
  double **Iq;
  if (VlTptr.size() == 0)
    Eq = NULL;
  else
    Eq = &VlTptr[0];

  if (VpTptr.size() == 0)
    Iq = NULL;
  else
    Iq = &VpTptr[0];

  if (compute_cone_generators (&num_pointy, &P, &num_lineality, &L, nvar,
                               VlT.getNrows(), Eq, VpT.getNrows(), Iq, 0) != 0)
  {
    return false;
  }

  // Add in lineality directions.
  Matrix Lmat(L, num_lineality, nvar);
  Lmat.normalize();
  Lmat.scale(scaling);
  Tcone.addMatrix(Lmat);
  Tcone.addMatrix(Lmat, -1.0);

  // Add in pointy directions.
  Matrix Pmat(P, num_pointy, nvar);
  Pmat.normalize();
  Pmat.scale(scaling);
  Tcone.addMatrix(Pmat);
  
  // Free memory that was allocated by compute_cone_generators().
  for (int i = 0; i < num_pointy; i++)
    free(P[i]);
  free(P);
  for (int i = 0; i < num_lineality; i++)
    free(L[i]);
  free(L);
  
  return true;
}
#endif


//PRIVATE
bool GssDirections::buildWithLapack(const Matrix& VpT, 
                                    const Matrix& VlT,
                                    Matrix& Tcone)
{ 
#if !defined(HAVE_LAPACK)
  return false;
#endif
  Tcone.clear();
  // Get orthonormal basis matrix Z for nullspace of VlT if Vlt is non-empty.
  // Otherwise Z becomes the idenity.
  Matrix ZT;
  if (!VlT.empty())
  {
    VlT.nullSpace(ZT, linConstr.getActiveTol());
    if (ZT.empty())
      return true;
  }
  else
  {
    // The are no equality constraints.  null(VlT) = I (trivially).
    ZT.setToIdentity(VpT.getNcols());
  }

  // Only equality constraints are present, so +/- Z generates space.
  if (VpT.empty()) 
  {
    // ZTscaled = Z^T * S
    Matrix ZTscaled(ZT, probDef.getVarScaling());
    Tcone.addMatrix(ZTscaled);
    Tcone.addMatrix(ZTscaled, -1.0);
    return true;
  }
  
  // Compute product VpTZ = VpT * Z.
  Matrix VpTZ;
  VpT.multMat(ZT, VpTZ, Matrix::TRANSPOSE);
  
  // Compute matrices R and N such that VpTZ*R = I and VptZ*N = 0;
  Matrix RT, NT;
  if (!VpTZ.getRightInvAndNullBasis(RT, NT, linConstr.getActiveTol()))
    return false;

  // Now form matrix Z*R.  Since we are working with transposes, we really form RT*ZT.
  Matrix ZRT(RT);
  ZRT.multMat(ZT);
  // Columns of Z*R are not necessarily unit length.
  ZRT.normalize();
  // Now scale columns and add to direction.
  ZRT.scale(probDef.getVarScaling());
  Tcone.addMatrix(ZRT, -1.0);

  // Now form Z*N.  Since we are working with transposes, we really form NT*ZT.
  // N may be empty if Z'*Vp is a square system, implying null(Z'*Vp) = {0}.
  if (!NT.empty())
  {
    Matrix ZNT(NT);
    ZNT.multMat(ZT); 
    // Normalizing unnecessary, N and Z both have orthonormal columns.
    // Now scale columns and add to direction.
    ZNT.scale(probDef.getVarScaling());
    Tcone.addMatrix(ZNT);
    Tcone.addMatrix(ZNT, -1.0);
  }
  
  return true;
}

//PRIVATE
void GssDirections::generateUnconstrained(Matrix& D)
{
  D.clear();
  const Vector& scaling = probDef.getVarScaling();

  for (int i = 0; i < nDimensions; i ++)
  {
    tmpVector = zero;    
    // Add +e_i
    tmpVector[i] = scaling[i];
    D.addRow(tmpVector);
    // Add -e_i
    tmpVector[i] = -1 * scaling[i];
    D.addRow(tmpVector);
  }
}

//PRIVATE
void GssDirections::generateForLinear(Matrix& D)
{
  D.clear();

  // Check if directions are cached.
  cacheIter = directionCache.find(constraintState);
  if (cacheIter != directionCache.end())
  {
    nCached++;
    D = cacheIter->second;
    return;
  }
  
  // Directions are not cached.  Form the normal cone.
  Matrix VpT;
  Matrix VlT;
  buildNormalCone(VpT, VlT);
  buildTangentCone(VpT, VlT, D);
  
  // If there are no tangent direction, Deltak must be reduced.
  if (D.empty())
  {
    // We still want to cache direction here--the work to figure out D is 
    // empty may be nontrivial.
    directionCache[constraintState] = D;
    return;
  }
  
  // Add normals.
  if (withNormals)
    addNormalDirections(VpT, VlT, D);

  if (withCompass)
    addCompassDirections(VlT, D);
  
  // Cache newly computed directions.
  directionCache[constraintState] = D;
}

//PRIVATE
void GssDirections::addNormalDirections(const Matrix& VpT,
                                        const Matrix& VlT,
                                              Matrix& D)
{
  // Add in normal directions, projecting with Z*ZT if necessary.
  if (VpT.empty())
    return;
  
  Matrix ND(VpT);
  if (!VlT.empty())
  {
#if !defined(HAVE_LAPACK)
    {
      cerr << "ERROR: Cannot add projected normals without LAPACK"
           << "       <GssDirections::addNormalDirections()>." << endl;
      throw "GSS Error";
    }
#endif
    Matrix ZT;
    VlT.nullSpace(ZT, linConstr.getActiveTol());
    if (ZT.empty()) // null(ZT) is empty, projecting to empty set.
      return;
    ND.multMat(ZT, Matrix::TRANSPOSE);
    ND.multMat(ZT);
  }
  
  ND.normalize();
  ND.scale(probDef.getVarScaling());
  D.addMatrix(ND);
}

//PRIVATE
void GssDirections::addCompassDirections(const Matrix& VlT,
                                               Matrix& D)
{
  //First form compass search directions.
  Matrix PCD;
  generateUnconstrained(PCD);

  //No project if necessary
  if (!VlT.empty())
  {    
#if !defined(HAVE_LAPACK)
    {
      cerr << "ERROR: Cannot add projected normals without LAPACK"
           << "       <GssDirections::addNormalDirections()>." << endl;
      throw "GSS Error";
    }
#endif
    Matrix ZT;
    VlT.nullSpace(ZT, linConstr.getActiveTol());
    if (ZT.empty()) // null(ZT) is empty, projecting to empty set.
      return;
    PCD.multMat(ZT, Matrix::TRANSPOSE);
    PCD.multMat(ZT);
  }
  
  PCD.normalize();
  PCD.scale(probDef.getVarScaling());
  D.addUniqueRows(PCD, linConstr.getActiveTol());
}

//PRIVATE
bool GssDirections::updateConstraintState(double newStep)
{
  // Ensure epsilon < epsilonMax.
  double epsilon = min(epsilonMax, newStep);
  // Determine indices of near active constraints wrt epsilon.
  vector< LinConstr::ActiveType > newState;
  linConstr.getActiveIneqIndices(xDistance, epsilon, newState);

  // If constraintState remains the same and is not updated return false.
  if (newState == constraintState)
    return false;
  
  // Otherwise update constraintState to newState and return true.
  constraintState = newState;
  return true;
}

//PRIVATE
void GssDirections::updateDirectionInfo(double newStep, bool isAppend)
{
  //Note: newStep \ge minStep if 1) isAppend = false, and 2) the corresponding
  //tangent cone is empty.  newStep is bumped up to minStep initially inside
  //computeNewDirections.  An empty tangent cone corresponds to a failed iteration
  //and hence we do not want to bump up to minStep in this case.
  if (isAppend)
  {
    // Only update new directions.
    int numnew = directionsMatrix.getNrows() - nDirections;
    if (numnew > 0)
      nAppend++;
    nDirections = directionsMatrix.getNrows();
    step.append(numnew, newStep);
    trueStep.append(numnew, -1);
    tag.insert(tag.end(), numnew, -1);
  }
  else
  {
    // Update all directions.
    nDirections = directionsMatrix.getNrows();
    step.assign(nDirections, newStep);
    trueStep.assign(nDirections, -1);
    tag.assign(nDirections, -1);
  }
  nMaxDirections=max(nMaxDirections, nDirections);
  // epsilonMin records smallest value of epsilon used so far.
  epsilonMin = min(epsilonMax, getSmallestStep());
}

double GssDirections::getSmallestStep() const
{
  double sstep = step.maxElement();
  for (int i=0; i < step.size(); i++)
  {
    if (step[i] >= stepTolerance)
      sstep = min(sstep, step[i]);
  }
  return sstep;
}

}     //-- namespace HOPSPACK
