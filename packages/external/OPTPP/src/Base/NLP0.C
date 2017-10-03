//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cfloat>
#include <cstring>
#else
#include <float.h>
#include <string.h>
#endif

#ifdef OPTPP_HAVE_MPI
#include "mpi.h"
#endif

#include "NLP0.h"
#include "TOLS.h"
#include "cblas.h"
#include "ioformat.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using namespace std;

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//----------------------------------------------------------------------------
// Evaluate the Hessian using finite differences
// No analytical gradients available so use function values
//----------------------------------------------------------------------------

SerialSymDenseMatrix<int,double> NLP0::FD2Hessian(SerialDenseVector<int,double> & sx) 
{
  double mcheps = DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  double hieps, eta;
  int i;
  int nr = getDim();

  double xtmpi, xtmpj;
  double fii, fij, fx;
  
  SerialDenseVector<int,double> fhi(nr), step(nr);
  SerialSymDenseMatrix<int,double> H(nr);

  // do we need this??? Dougm xc = getXc();
  fx = getF();

  for (i=0; i<nr; i++) {
    hieps = max(mcheps,fcn_accrcy(i));
    eta   = pow(hieps,0.333333);
    step(i) = eta*max(fabs(mem_xc(i)),sx(i));
    step(i) = copysign(step(i),mem_xc(i));
    xtmpi = mem_xc(i);
    mem_xc(i) = xtmpi + step(i);
    fhi(i) = evalF(mem_xc);
    mem_xc(i) = xtmpi;
  }
  
  for (i=0; i<nr; i++) {
    xtmpi = mem_xc(i);
    mem_xc(i) = mem_xc(i) + step(i)*2.0;
    fii = evalF(mem_xc); 
    H(i,i) = ((fx - fhi(i)) + (fii - fhi(i))) / (step(i)*step(i));
    mem_xc(i) = xtmpi + step(i);
    for (int j=i+1; j<nr; ++j) {
      xtmpj = mem_xc(j);
      mem_xc(j) = mem_xc(j) + step(j);
      fij = evalF(mem_xc);
      H(i,j) = ((fx - fhi(i)) + (fij - fhi(j))) / (step(i)*step(j));
      mem_xc(j) = xtmpj;
    }
    mem_xc(i) = xtmpi;
  } 

  return H;
}

// Compute gradient using backward finite differences

SerialDenseVector<int,double> NLP0::BDGrad(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& x,
			  double& fx, SerialDenseVector<int,double>& grad)
{
  int i, j, gradStart, gradEnd, nBcasts;
  double fminus, hi;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  const int tmpSize = (int) ceil((double) ndim/nprocs);
  double *tmpGradMinus = new double[tmpSize];
  SerialDenseVector<int,double> xcurrent(x.length());
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SpecOption SpecPass = getSpecOption();
  CompoundConstraint* constraints = getConstraints();
  bool scaleStep = false;

#ifdef OPTPP_HAVE_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::BDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(x);
      setSpecOption(SpecPass);
#ifdef OPTPP_HAVE_MPI
      if (nprocs > 1)
	MPI_Bcast(&fx, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);
#endif
    }
    gradStart = 1;
    gradEnd = min(ndim, nprocs-1);
    nBcasts = min(ndim, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    gradStart = nprocs;
    gradEnd = ndim;
    nBcasts = min(gradEnd-gradStart+1, nprocs);
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    nBcasts = min(ndim, nprocs);
    if (SpecPass != NoSpec) {
    cerr << "NLP0::BDGrad: Invalid speculative gradient option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute my piece of the gradient.

  for (i=me+gradStart-1; i<=gradEnd-1; i+=nprocs) {

    xcurrent = perturbX(i, x, sx(i), *constraints, 
			fcn_accrcy(i), hi, scaleStep, BackwardDiff);
    setSpecOption(NoSpec);
    fminus = evalF(xcurrent);
    setSpecOption(SpecPass);
#ifdef OPTPP_HAVE_MPI
    if (SpecPass == Spec1)
      MPI_Bcast(&fx, 1, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
#endif
    grad(i) = (fx - fminus) / hi;
  }

  // Share my piece of the gradient with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+gradStart; j<=gradEnd; j+=nprocs)
	tmpGradMinus[(j-me-gradStart)/nprocs] = grad(j-1);

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpGradMinus, tmpSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=nprocs)
	grad(j-1) = tmpGradMinus[(j-i-gradStart)/nprocs];
    }
  }

  if (tmpGradMinus != NULL)
    delete[] tmpGradMinus;

  return grad;
}

// Compute gradient using forward finite differences

SerialDenseVector<int,double> NLP0::FDGrad(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& x,
			  double& fx, SerialDenseVector<int,double>& grad) 
{
  int i, j, gradStart, gradEnd, nBcasts;
  double fplus, hi;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  const int tmpSize = (int) ceil((double) ndim/nprocs);
  double *tmpGradPlus = new double[tmpSize];
  SerialDenseVector<int,double> xcurrent(x.length());
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SpecOption SpecPass = getSpecOption();
  CompoundConstraint* constraints = getConstraints();
  bool scaleStep = false;

#ifdef OPTPP_HAVE_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::FDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif
  
  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(x);
      setSpecOption(SpecPass);
#ifdef OPTPP_HAVE_MPI
      if (nprocs > 1)
	MPI_Bcast(&fx, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);
#endif
    }
    gradStart = 1;
    gradEnd = min(ndim, nprocs-1);
    nBcasts = min(ndim, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    gradStart = nprocs;
    gradEnd = ndim;
    nBcasts = min(gradEnd-gradStart+1, nprocs);
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    nBcasts = min(ndim, nprocs);
    if (SpecPass != NoSpec) {
    cerr << "NLP0::FDGrad: Invalid speculative gradient option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the gradient.

  for (i=me+gradStart-1; i<=gradEnd-1; i+=nprocs) {
    xcurrent = perturbX(i, x, sx(i), *constraints, 
			fcn_accrcy(i), hi, scaleStep, ForwardDiff);
    setSpecOption(NoSpec);
    fplus = evalF(xcurrent);
    setSpecOption(SpecPass);
#ifdef OPTPP_HAVE_MPI
    if (SpecPass == Spec1)
      MPI_Bcast(&fx, 1, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
#endif
    grad(i) = (fplus - fx) / hi;
  }

  // Share my piece of the gradient with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+gradStart; j<=gradEnd; j+=nprocs)
	tmpGradPlus[(j-me-gradStart)/nprocs] = grad(j-1);

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpGradPlus, tmpSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=nprocs)
	grad(j-1) = tmpGradPlus[(j-i-gradStart)/nprocs];

    }
  }

  if (tmpGradPlus != NULL)
    delete[] tmpGradPlus;

  return grad;
}

// Compute gradient using central differences

SerialDenseVector<int,double> NLP0::CDGrad(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& x,
			  double& fx, SerialDenseVector<int,double>& grad) 
{
  int i, gradStart, gradEnd, myStart, inc, nBcasts;
  double fplus, fminus, xtmp, h1, h2;
  int j, tmpSize;

  int me = 0;
  int nprocs = 1;
  int ndim = getDim();
  SerialDenseVector<int,double> xcurrent(x.length());
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SpecOption SpecPass = getSpecOption();
  CompoundConstraint* constraints = getConstraints();
  bool scaleStep = false;

#ifdef OPTPP_HAVE_MPI

  char buffer[MPI_MAX_ERROR_STRING];
  int error, resultlen, flag;

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "NLP0::CDGrad: MPI Error - " << buffer << endl;
    }
  }

#endif
  
  // Set loop endpoints, f, and x according to which pass of
  // speculative gradient evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      setSpecOption(NoSpec);
      fx = evalF(x);
      setSpecOption(SpecPass);
    }
    gradStart = 1;
    gradEnd = min(ndim, (int) floor((double) (nprocs-1)/2));
    if (nprocs > 1)
      inc = (int) floor((double) (nprocs-1)/2);
    else
      inc = 1;
    nBcasts = min(ndim, (int) floor((double) (nprocs-1)/2));
  }
  else if (SpecPass == Spec2) {
    gradStart = (int) ceil((double) nprocs/2);
    gradEnd = ndim;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(gradEnd-gradStart+1, (int) floor((double) nprocs/2));
  }
  else {
    gradStart = 1;
    gradEnd = ndim;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(ndim, (int) floor((double) nprocs/2));
    if (SpecPass != NoSpec) {
    cerr << "NLP0::FDGrad: Invalid speculative gradient option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the gradient.

  myStart = (int) floor((double) me/2) + gradStart;

  for (i=myStart-1; i<=gradEnd-1; i+=inc) {

#ifdef OPTPP_HAVE_MPI
    if (nprocs > 1) {

      // For multiple processors, even processors look forward, and
      // odd look backward.

      // ***This still does not respect bounds!  Fix later.***

      if (me%2 == 0)
	xcurrent(i) = xtmp + hi;
      else
	xcurrent(i) = xtmp - hi;

      setSpecOption(NoSpec);
      grad(i) = evalF(xcurrent)/(2*hi);
      setSpecOption(SpecPass);
    }
    else {
      // Otherwise, do the same as in the serial case.
#endif
      xcurrent = perturbX(i, x, sx(i), *constraints, 
			fcn_accrcy(i), h1, scaleStep, CentralDiff1);
      setSpecOption(NoSpec);
      fplus = evalF(xcurrent);
      setSpecOption(SpecPass);

      h2 = h1;
      xcurrent = perturbX(i, x, sx(i), *constraints, 
			fcn_accrcy(i), h2, scaleStep, CentralDiff2);
      setSpecOption(NoSpec);
      fminus = evalF(xcurrent);
      setSpecOption(SpecPass);

      grad(i)= (fplus - fminus) / (h1+h2);
#ifdef OPTPP_HAVE_MPI
    }
    xcurrent(i) = xtmp;
#endif
  }

  if (nprocs > 1) {

    if (nprocs%2 == 0)
      tmpSize = (int) ceil((double) (2*ndim)/nprocs);
    else
      // If there are an odd number of processors, the last one doesn't
      // count.
      tmpSize = (int) ceil((double) (2*ndim)/(nprocs-1));

    double *tmpGradPlus = new double[tmpSize];
    double *tmpGradMinus = new double[tmpSize];

    for (i=0; i<nBcasts; i++) {

      for (j=myStart; j<=gradEnd; j+=inc) {
	if (me%2 == 0)
	  tmpGradPlus[(j-myStart)/inc] = grad(j-1);
	else
	  tmpGradMinus[(j-myStart)/inc] = grad(j-1);
      }

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpGradPlus, tmpSize, MPI_DOUBLE, 2*i, MPI_COMM_WORLD);
      MPI_Bcast(tmpGradMinus, tmpSize, MPI_DOUBLE, (2*i)+1, MPI_COMM_WORLD);
#endif

      for (j=i+gradStart; j<=gradEnd; j+=inc)
	grad(j-1) = tmpGradPlus[(j-i-gradStart)/inc] - 
	                     tmpGradMinus[(j-i-gradStart)/inc];
    }
    if (tmpGradPlus != NULL)
      delete[] tmpGradPlus;
    if (tmpGradMinus != NULL)
      delete[] tmpGradMinus;
  }
  return grad;
}

// Compute gradient of nonlinear constraints using backward finite differences
SerialDenseMatrix<int,double> NLP0::CONBDGrad(const SerialDenseVector<int,double>& sx) 
{
  double mcheps = DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps;
  SerialDenseVector<int,double> fminus(ncnln), fx(ncnln);
  
  n = dim;
  SerialDenseMatrix<int,double> grad(n,ncnln), gtmp(ncnln,n);
  fx = evalCF(mem_xc);
  //fx = getConstraintValue();

  for (i=0; i<n; i++) {
    hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
    hi = hieps*max(fabs(mem_xc(i)),sx(i));
    hi = copysign(hi,mem_xc(i));
    xtmp = mem_xc(i);
    mem_xc(i) = xtmp - hi;
    fminus = evalCF(mem_xc);
    for(int j = 0; j<ncnln; j++)
      { gtmp(j,i) = fx(j);
	gtmp(j,i) -= fminus(j);
	gtmp(j,i) *= 1/hi;
      }
    mem_xc(i) = xtmp;
  }
  // grad = gtmp.trans();
  for(i=0; i<ncnln; i++)
    for(int j=0;j<n;j++)
      grad(i,j) = gtmp(i,j);
  return grad;
}

// Compute gradient of nonlinear constraints using forward finite differences
SerialDenseMatrix<int,double> NLP0::CONFDGrad(const SerialDenseVector<int,double>& sx) 
{
  double mcheps = DBL_EPSILON;
  SerialDenseVector<int,double>  fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps;
  SerialDenseVector<int,double> fx(ncnln), fplus(ncnln);
  
  n = dim;
  SerialDenseVector<int,double> xcurrent(n);
  SerialDenseMatrix<int,double> grad(n,ncnln), gtmp(ncnln,n);
  xcurrent = getXc();
  fx = evalCF(xcurrent);
  //fx = getConstraintValue();

  for (i=0; i<n; i++) {
    hieps = sqrt(max(mcheps,fcn_accrcy(i) ));
    hi = hieps*max(fabs(xcurrent(i)),sx(i));
    hi = copysign(hi,xcurrent(i));
    xtmp = xcurrent(i);
    xcurrent(i) = xtmp + hi;
    fplus = evalCF(xcurrent);
    for(int j = 0; j<ncnln; j++)
      { gtmp(j,i) = fplus(j);
	gtmp(j,i) -= fx(j);
	gtmp(j,i) *= 1/hi;
      }
    xcurrent(i) = xtmp;
  }
  // grad = gtmp.t();
  for(i=0; i<ncnln; i++)
    for(int j=0;j<n;j++)
      grad(i,j) = gtmp(j,i);
  return grad;
}

// Compute gradient of nonlinear constraints using central differences
SerialDenseMatrix<int,double> NLP0::CONCDGrad(const SerialDenseVector<int,double>& sx) 
{
  double mcheps = DBL_EPSILON;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  int i, n;
  double xtmp, hi, hieps; 
  SerialDenseVector<int,double> fplus(ncnln), fminus(ncnln);
  
  n = dim;
  SerialDenseMatrix<int,double> grad(n, ncnln), gtmp(ncnln,n);

  for (i=0; i<n; i++) {

    hieps = max(mcheps,fcn_accrcy(i) );
    hieps = pow(hieps,0.333333);

    hi = hieps*max(fabs(mem_xc(i)),sx(i));
    hi = copysign(hi,mem_xc(i));

    xtmp   = mem_xc(i);
    mem_xc(i)  = xtmp + hi;
    fplus  = evalCF(mem_xc);

    mem_xc(i)  = xtmp - hi;
    fminus = evalCF(mem_xc);

        for(int j = 0; j<ncnln; j++)
	  { gtmp(j,i) = fplus(j);
	    gtmp(j,i) -= fminus(j);
	    gtmp(j,i) *= 1/(2*hi);
      }
    mem_xc(i) = xtmp;
  }
  // grad = gtmp.t();
  for(i=0; i<ncnln; i++)
    for(int j=0;j<n;j++)
      grad(i,j) = gtmp(j,i);
  return grad;
}

SerialDenseVector<int,double> NLP0::perturbX(int& index, const SerialDenseVector<int,double>& xcurrent, const double& sx, CompoundConstraint& constraints, double& fcn_accrcy, double& hi, bool& scaleStep, const DerivOption typeDiff)
{
  SerialDenseVector<int,double> xPert(xcurrent);
  SerialDenseVector<int,double> dist_to_lower(xcurrent.length());
  SerialDenseVector<int,double> dist_to_upper(xcurrent.length());

  double mcheps = DBL_EPSILON;
  double hieps;

  if (typeDiff == ForwardDiff || typeDiff == BackwardDiff ||
      typeDiff == CentralDiff1) {
    hieps = sqrt(max(mcheps, fcn_accrcy));
    if (typeDiff == CentralDiff1)
      hieps = pow(hieps, 0.333333);
    hi = hieps * max(fabs(xcurrent(index)), sx);
    hi = copysign(hi, xcurrent(index));
  }

  if (typeDiff == ForwardDiff || typeDiff == CentralDiff1) {
    xPert(index) = xcurrent(index) + hi;
    constraints.computeDistanceToBounds(xPert, dist_to_lower, dist_to_upper);
    if ( (xcurrent(index) < 0.0 && dist_to_lower(index) < 0.0) ||
	 (xcurrent(index) >= 0.0 && dist_to_upper(index) < 0.0) ) {
      xPert(index) = xcurrent(index) - hi;
      constraints.computeDistanceToBounds(xPert, dist_to_lower, dist_to_upper);
      if ( (xcurrent(index) < 0.0) && (dist_to_upper(index) >= 0.0) ||
	   (xcurrent(index) >= 0.0) && (dist_to_lower(index) >= 0.0) )
	hi = -hi;
      else
	scaleStep = true;
    }
  }

  if (typeDiff == BackwardDiff || typeDiff == CentralDiff2) {
    xPert(index) = xcurrent(index) - hi;
    constraints.computeDistanceToBounds(xPert, dist_to_lower, dist_to_upper);
    if ( (xcurrent(index) < 0.0 && dist_to_upper(index) < 0.0) ||
	 (xcurrent(index) >= 0.0 && dist_to_lower(index) < 0.0) ) {
      xPert(index) = xcurrent(index) + hi;
      constraints.computeDistanceToBounds(xPert, dist_to_lower, dist_to_upper);
      if ( (xcurrent(index) < 0.0) && (dist_to_lower(index) >= 0.0) ||
	   (xcurrent(index) >= 0.0) && (dist_to_upper(index) >= 0.0) )
	hi = -hi;
      else
	scaleStep = true;
    }
  }

  if ( (typeDiff == ForwardDiff || typeDiff == BackwardDiff ||
	typeDiff == CentralDiff1) && scaleStep) {
    xPert(index) = xcurrent(index);
    constraints.computeDistanceToBounds(xPert, dist_to_lower, dist_to_upper);
    if (dist_to_lower(index) < dist_to_upper(index))
      hi = dist_to_upper(index);
    else
      hi = -dist_to_lower(index);
    xPert(index) = xcurrent(index) + hi;
  }

  return xPert;
}

//-------------------------------------------------------------------------
// Output Routines
//-------------------------------------------------------------------------

void NLP0::printState(const char * s) 
{ // Print out current state: x current, gradient and Function value
  cout << "\n\n=========  " << s << "  ===========\n\n";
  cout << "\n    i\t   x  \t      grad   \t\t fcn_accrcy \n\n";
  for (int i=0; i<dim; i++) 
    cout << d(i,5) << "\t" << e(mem_xc(i),12,4)<< "\t\t"
         << e(mem_fcn_accrcy(i),12,4) << "\n";
  cout <<"Function Value     = " << e(fvalue,12,4) << "\n";
  //cout <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  cout <<"\n\n===================================================\n\n";
}

void NLP0::fPrintState(ostream *nlpout, const char * s) 
{ // Print out current state: x current, gradient and Function value
  (*nlpout) << "\n\n=========  " << s << "  ===========\n\n";
  (*nlpout) << "\n    i\t   x  \t      grad   \t\t fcn_accrcy \n\n";
  for (int i=0; i<dim; i++) 
    (*nlpout) << d(i,5) << "\t" << e(mem_xc(i),12,4) << "\t\t"
              << e(mem_fcn_accrcy(i),12,4) << "\n";
  (*nlpout) <<"Function Value     = " << e(fvalue,12,4) << "\n";
 // (*nlpout) <<"Function Accuracy  = " << e(mem_fcn_accrcy,12,4) << "\n";
  (*nlpout) <<"\n\n===================================================\n\n";
}

void NLP0::saveState() 
{ // Save current state: x current, gradient and Function value
  cout << dim << "\n";
  for (int i=0; i<dim; i++) 
     cout << e(mem_xc(i),24,16) << "\t" << e(mem_fcn_accrcy(i),24,16) << "\n";
  cout << e(fvalue,24,16) << "\n" 
	<< nlp_name << "\n"
	<< nfevals << "\n"
	<< is_expensive << "\n"
	<< debug_ << "\n"
	<< e(function_time,24,16) << "\n";
}

bool NLP0::hasConstraints()
{
  bool nonempty = false;
  if( constraint_)
     if (constraint_->getNumOfSets())
     nonempty = true;
  return nonempty;
}

void NLP0::printConstraints()
{
   
  constraint_->printConstraints();
  cout <<"\n\n===================================================\n\n";

}

} // namespace OPTPP
