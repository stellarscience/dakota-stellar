
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

#include "LSQNLF.h"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

using namespace std;
using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;


namespace OPTPP {

//-------------------------------------------------------------------------
// LSQNLF method routines.
// Derived from the NLP2 class.
//-------------------------------------------------------------------------
void LSQNLF::reset() // Reset parameter values  
{
  init_flag = false;    
  nfevals   = ngevals = 0;
#ifdef OPTPP_HAVE_MPI
  SpecFlag  = Spec1;
#else
  SpecFlag  = NoSpec;
#endif
  application.reset();
}

void LSQNLF::initFcn() // Initialize Function
{
  if (init_flag == false)  {
    init_fcn(dim, mem_xc);
    init_flag = true;
  }
  else  {
    cout << "LSQNLF:initFcn: Warning - initialization called twice\n";
    init_fcn(dim, mem_xc);
  }
}

double LSQNLF::evalF() // Evaluate Function
{
  int result   = 0;
  double time0 = get_wall_clock_time();

  if(fcn0_v != NULL){
     if(SpecFlag == NoSpec) {
        if (!application.getLSQF(mem_xc,fvector)) {
           fcn0_v(dim, mem_xc, fvector, result, vptr);
	   application.lsq_update(NLPFunction,dim,lsqterms_,mem_xc,fvector);
           nfevals++;
           Jacobian_current = false;
        }
     } 
     else {
       SpecFlag = Spec1;
       (void) evalG();
       SpecFlag = Spec2;
     }
  }
  else if(fcn1_v != NULL){
     SerialDenseMatrix<int,double> jac(lsqterms_,dim); 
     if (!application.getLSQF(mem_xc,fvector)) {
        fcn1_v(NLPFunction, dim, mem_xc, fvector, jac, result, vptr);
        application.lsq_update(result,dim,lsqterms_,mem_xc,fvector,jac);
        nfevals++;
        Jacobian_current = false;
     }
  } 
  else{
    cerr << "Error: A function has not been declared. \n";
    exit(1);
  }

  fvalue = fvector.dot(fvector);  

  setFcnResidual(fvector);

  function_time = get_wall_clock_time() - time0;
  if(debug_)
  cout << "LSQNLF::evalF( )\n" 
       << "nfevals       = " << nfevals << "\n"
       << "fvalue        = " << fvalue << "\n"
       << "function time = " << function_time << "\n";
  return fvalue;

}

double LSQNLF::evalF(const SerialDenseVector<int,double>& x) // Evaluate Function at x
{
  int result = 0;
  SerialDenseVector<int,double> fx(lsqterms_);
  double ftmp, time0 = get_wall_clock_time();

  if (fcn0_v != NULL){
     if (SpecFlag == NoSpec) {
         if (!application.getLSQF(x,fx)) {
            fcn0_v(dim, x, fx, result, vptr);
	    application.lsq_update(NLPFunction,dim,lsqterms_,x,fx);
            nfevals++;
            Jacobian_current = false;
         }
     } 
     else {
       SpecFlag = Spec1;
       (void) evalG(x);
       fx = specLSQF;
       SpecFlag = Spec2;
     }
  } 
  else if(fcn1_v != NULL){
     SerialDenseMatrix<int,double> gx(lsqterms_,dim);
     if (!application.getLSQF(x,fx)) {
        fcn1_v(NLPFunction, dim, x, fx, gx, result, vptr);
        application.lsq_update(result,dim,lsqterms_,x,fx,gx);
        nfevals++;
        Jacobian_current = false;
     }
  } 
  else{
    cerr << "Error: A function has not been declared. \n";
    exit(1);
  }

  ftmp = fx.dot(fx);  
  /*
   * In the search strategy routines, an evalF(xplus) is
   * followed with a setF(xplus) and evalG().  The setF
   * call updates the function value.  A similiar routine
   * is needed for the function residuals otherwise we would have
   * to recompute the function residuals.  My kluge is
   * to update setFcnResidual each time the function is evaluated.
   */
  setFcnResidual(fx);

  function_time = get_wall_clock_time() - time0;

  if(debug_)
  cout << "LSQNLF::evalF(x)\n" 
       << "nfevals       = " << nfevals << "\n"
       << "fvalue        = " << ftmp << "\n"
       << "function time = " << function_time << "\n";
  return ftmp;

}

SerialDenseVector<int,double> LSQNLF::evalG() 
{

  int result;

  if(fcn0_v != NULL){
    SerialDenseVector<int,double> sx(dim);
    sx = 1.0;

    if (!application.getLSQF(mem_xc,fvector)) {
        fcn0_v(dim, mem_xc, fvector, result, vptr);
	application.lsq_update(NLPFunction,dim,lsqterms_,mem_xc,fvector);
	nfevals++;
    }
    else
	fvector = getFcnResidual();

    if(finitediff == ForwardDiff)
       Jacobian_ = LSQFDJac(sx, mem_xc, fvector, partial_jac);
    else if(finitediff == BackwardDiff)
       Jacobian_ = LSQBDJac(sx, mem_xc, fvector, partial_jac);
    else if(finitediff == CentralDiff)
       Jacobian_ = LSQCDJac(sx, mem_xc, fvector, partial_jac);
    else{
       cout << "LSQNLF::evalG: Unrecognized difference option\n";
       cout << "LSQNLF::evalG: Using forward difference option\n";
       Jacobian_ = LSQFDJac(sx, mem_xc, fvector, partial_jac);
    }
    // mem_grad = 2*Jacobian_.t()*fvector;
    mem_grad.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, fvector,0.0);

  }
  else if(fcn1_v != NULL){
    if (!application.getLSQF(mem_xc,fvector) || !application.getLSQJac(mem_xc,Jacobian_)){
       int mode = NLPGradient;
       if (!application.getLSQF(mem_xc,fvector)) {
	 mode = NLPFunction | NLPGradient;
	 nfevals++;
       }
       fcn1_v(mode, dim, mem_xc, fvector, Jacobian_, result, vptr);
       application.lsq_update(result,dim,lsqterms_,mem_xc,fvector,Jacobian_);
       //mem_grad = 2*Jacobian_.t()*fvector;  
 mem_grad.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, fvector,0.0);       
ngevals++;
    }
    else {
       /*
        * If an external package overrode the mode
        * and computed the function and gradient simultaneously
        * the if block wouldn't be activated and hence the
        * gradient would never change.  We grab the stored
        * residual values.
        */ 
      // mem_grad = 2*Jacobian_.t()*getFcnResidual();
      mem_grad.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, getFcnResidual(),0.0);
    }
  }

  Jacobian_current = true;
  return mem_grad;

}

SerialDenseVector<int,double> LSQNLF::evalG(const SerialDenseVector<int,double>& x) 
{
  int result = 0;
  SerialDenseVector<int,double> fx(lsqterms_), gtmp(dim);
  SerialDenseMatrix<int,double> gx(lsqterms_,dim),HessTmp(Hessian.numRows(),Hessian.numCols());

  if(fcn0_v != NULL){
    SerialDenseVector<int,double> sx(dim);
    sx = 1.0;

    if(SpecFlag == NoSpec){
       if(!application.getLSQF(x,specLSQF)){
         fcn0_v(dim, x, specLSQF, result, vptr);
         nfevals++;
       }
    }

    if(finitediff == ForwardDiff)
       gx = LSQFDJac(sx, x, specLSQF, partial_jac);
    else if(finitediff == BackwardDiff)
       gx = LSQBDJac(sx, x, specLSQF, partial_jac);
    else if(finitediff == CentralDiff)
       gx = LSQCDJac(sx, x, specLSQF, partial_jac);
    else{
       cout << "LSQNLF::evalG: Unrecognized difference option\n";
       cout << "LSQNLF::evalG: Using forward difference option\n";
       gx =  LSQFDJac(sx, x, specLSQF, partial_jac);
    }
    // gtmp    = 2*gx.t()*specLSQF;
    gtmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, specLSQF, 0.0);
    // Hessian <<  ( gx.t()*gx ) * 2.0;
    // Hessian.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, gx, 0.0);
    HessTmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, gx, 0.0);
    for(int i=0; i<Hessian.numRows();i++)
      for(int j=0; j<=i;j++)
	{Hessian(i,j) = HessTmp(i,j);}
    // Hessian = HessTmp;

  }
  else if(fcn1_v != NULL){
    if (!application.getLSQF(x,specLSQF) || !application.getLSQJac(x,gx)){
       int mode = NLPGradient;
       if (!application.getLSQF(x,specLSQF)) {
	 mode = NLPFunction | NLPGradient;
	 nfevals++;
       }
       fcn1_v(mode, dim, x, fx, gx, result, vptr);
       application.lsq_update(result,dim,lsqterms_,x,fx,gx);
       // gtmp        = 2*gx.t()*fx;
       gtmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, fx, 0.0);
       ngevals++;
    }
    else {
       /*
        * If an external package overrode the mode
        * and computed the function and gradient simultaneously
        * the if block wouldn't be activated and hence the
        * gradient would never change.  We grab the stored
        * residual values.
        */
      // gtmp = 2*gx.t()*getFcnResidual();
      gtmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, getFcnResidual(), 0.0);
    }
    // Hessian <<  ( gx.t()*gx ) * 2.0;
    // Hessian.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, gx, 0.0);
    HessTmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, gx, gx, 0.0);
    for(int i=0; i<=Hessian.numRows();i++)
      for(int j=0; j<=i;j++)
	{Hessian(i,j) = HessTmp(i,j);}
  }

  Jacobian_current = true;
  return gtmp;
}

SerialSymDenseMatrix<int,double> LSQNLF::evalH() 
{
  SerialDenseMatrix<int,double> HessTmp(Hessian.numRows(),Hessian.numCols());
  if (!application.getLSQJac(mem_xc,Jacobian_))
    (void) evalG();
  // Hessian << (Jacobian_.t()*Jacobian_)*2.0;
  // Hessian.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, Jacobian_, 0.0);
    HessTmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, Jacobian_, 0.0);
    for(int i=0; i<Hessian.numRows();i++)
      for(int j=0; j<=i;j++)
	{Hessian(i,j) = HessTmp(i,j);}

  return Hessian;
}

SerialSymDenseMatrix<int,double> LSQNLF::evalH(SerialDenseVector<int,double>& x) 
{
  SerialDenseMatrix<int,double> gx(lsqterms_,dim);

  if (!application.getLSQJac(x,gx))
    (void) evalG(x);
  return Hessian;
}

void LSQNLF::eval()
{
  (void) evalG();
  SerialDenseMatrix<int,double> HessTmp(Hessian.numRows(),Hessian.numCols());
  fvalue = fvector.dot(fvector);  
  setFcnResidual(fvector);

  // Hessian << (Jacobian_.t()*Jacobian_)*2.0;
  //Hessian.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, Jacobian_, 0.0);
 HessTmp.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 2.0, Jacobian_, Jacobian_, 0.0);
    for(int i=0; i<Hessian.numRows();i++)
      for(int j=0; j<=i;j++)
	{Hessian(i,j) = HessTmp(i,j);}

}

double LSQNLF::evalLagrangian(const SerialDenseVector<int,double>& xc , 
                          SerialDenseVector<int,double>& multiplier,
                          const SerialDenseVector<int,double>& type)
{
   double result = evalF(xc);
   if( hasConstraints()){
     SerialDenseVector<int,double> resid(constraint_->evalResidual(xc));
     //     SerialDenseVector<int,double> resid(constraint_->getNumOfCons());
     //	resid = constraint_->evalResidual(xc);
      result  -=  resid.dot(multiplier);
   }
   return result;
}

SerialDenseVector<int,double> LSQNLF::evalLagrangianGradient(const SerialDenseVector<int,double>& xc, 
                                          const SerialDenseVector<int,double>& multiplier,
					  const SerialDenseVector<int,double>& type) 
{
   SerialDenseVector<int,double> grad = evalG(xc);
   SerialDenseVector<int,double> tmult2(grad.length());
   if(hasConstraints()){
     SerialDenseVector<int,double> tmult(multiplier.length());
     tmult = multiplier;
      for (int i = 0; i < getNumOfCons(); i++){
         if(type(i) == NLineq || type(i) == Lineq)
            tmult(i)*= -1;
      }
      tmult2.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,constraint_->evalGradient(xc),tmult,0.0);      
      grad += tmult2;

   }
   return grad;
}

SerialDenseVector<int,double> LSQNLF::evalCF(const SerialDenseVector<int,double>& x) // Evaluate Function at x
{
 
  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
       << "for nonlinear constraints.  Please select a different   \n" 
       << "NLF object, say an FDNLF.  " 
       << endl;
  exit(1);
  return(x);

}

SerialDenseMatrix<int,double> LSQNLF::evalCG(const SerialDenseVector<int,double>& x) // Evaluate the gradient at x
{

  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
       << "for nonlinear constraints.  Please select a different   \n" 
       << "NLF object, say an FDNLF.  " 
       << endl;
  exit(1);
  return(x);

}

SerialSymDenseMatrix<int,double> LSQNLF::evalCH(SerialDenseVector<int,double>& x) // Evaluate the Hessian at x
{
  
    cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
         << "for nonlinear constraints.  Please select a different   \n" 
         << "NLF object, say an FDNLF.  " 
         << endl;
  
    exit(1);
    SerialSymDenseMatrix<int,double> H(dim);
    H = 0;
    return(H);

}

OptppArray<SerialSymDenseMatrix<int,double> > LSQNLF::evalCH(SerialDenseVector<int,double>& x, int darg) // Evaluate the Hessian at x
{
    cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
       << "for nonlinear constraints.  Please select a different   \n" 
       << "NLF object, say an FDNLF.  " 
       << endl;
    exit(1);
    SerialSymDenseMatrix<int,double> H(dim);
    H = 0;
    OptppArray<SerialSymDenseMatrix<int,double> > HH;
    HH.append(H);
    return(HH);
  
}

void LSQNLF::evalC(const SerialDenseVector<int,double>& x)
{
  cerr << "Error: OPT++ does not support the Gauss Newton operator \n"
       << "for nonlinear constraints.  Please select a different   \n" 
       << "NLF object, say an FDNLF.  " 
       << endl;
  exit(1);
}

// Compute Jacobian of function vector using backward finite differences
SerialDenseMatrix<int,double> LSQNLF::LSQBDJac(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& xc,
	                SerialDenseVector<int,double>& fx, SerialDenseMatrix<int,double>& jac)
{
  int i, j, k, jacStart, jacEnd, nBcasts;
  double xtmp, hi, hieps;
  SerialDenseVector<int,double> fminus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  const int tmpSize = (int) ceil((double) n/nprocs);
  double *tmpJacMinus = new double[tmpSize*lsqterms_];
  double *tmpF = new double[lsqterms_];
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SerialDenseVector<int,double> xcurrent(xc.length());
  xcurrent  = xc;
  double mcheps = DBL_EPSILON;
  SpecOption SpecPass = getSpecOption();

#ifdef OPTPP_HAVE_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQBDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
#ifdef OPTPP_HAVE_MPI
      if (nprocs > 1) {
	for (i=0; i<lsqterms_; i++)
	  tmpF[i] = fx(i);
	MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, me, MPI_COMM_WORLD);
      }
#endif
    }
    jacStart = 1;
    jacEnd = min(n, nprocs-1);
    nBcasts = min(n, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    jacStart = nprocs;
    jacEnd = n;
    nBcasts = min(jacEnd-jacStart+1, nprocs);
  }
  else {
    jacStart = 1;
    jacEnd = n;
    nBcasts = min(n, nprocs);
    if (SpecPass != NoSpec) {
    cerr << "LSQNLF::LSQBDJac: Invalid speculative Jacobian option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  for (i=me+jacStart-1; i<=jacEnd-1; i+=nprocs) {

    hieps = sqrt(max(mcheps,fcn_accrcy(i)));
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp          = xcurrent(i);
    xcurrent(i)   = xtmp - hi;

    fcn0_v(n, xcurrent, fminus, result, vptr);
#ifdef OPTPP_HAVE_MPI
    if (SpecPass == Spec1) {
      MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
      for (j=0; j<lsqterms_; j++)
	fx(j) = tmpF[j];
    }
#endif
    // jac.Column(i) = (fx - fminus) / hi;
      for(int j=0; j<jac.numRows(); j++)
      {jac(j,i)=(fx(j)-fminus(j))/hi;}
    xcurrent(i)   = xtmp;
  }

  // Share my piece of the Jacobian with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  tmpJacMinus[k+lsqterms_*((j-me-jacStart)/nprocs)] = jac(k+1,j);
      }

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpJacMinus, tmpSize*lsqterms_, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  jac(k,j-1) = tmpJacMinus[k+lsqterms_*((j-i-jacStart)/nprocs)];
      }
    }
  }

  if (tmpJacMinus != NULL)
    delete[] tmpJacMinus;
  if (tmpF != NULL)
    delete[] tmpF;

  return jac;
}

// Compute Jacobian of function vector using forward finite differences
SerialDenseMatrix<int,double> LSQNLF::LSQFDJac(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& xc,
	                SerialDenseVector<int,double>& fx, SerialDenseMatrix<int,double>& jac)
{
  int i, j, k, jacStart, jacEnd, nBcasts;
  double xtmp, hi, hieps;
  SerialDenseVector<int,double> fplus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  const int tmpSize = (int) ceil((double) n/nprocs);
  double *tmpJacPlus = new double[tmpSize*lsqterms_];
  double *tmpF = new double[lsqterms_];
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SerialDenseVector<int,double> xcurrent(xc.length());
  xcurrent   = xc;
  double mcheps = DBL_EPSILON;
  SpecOption SpecPass = getSpecOption();

#ifdef OPTPP_HAVE_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQFDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
#ifdef OPTPP_HAVE_MPI
      if (nprocs > 1) {
	for (i=0; i<lsqterms_; i++)
	  tmpF[i] = fx(i);
	MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, me, MPI_COMM_WORLD);
      }
#endif
    }
    jacStart = 1;
    jacEnd = min(n, nprocs-1);
    nBcasts = min(n, nprocs-1);
  }
  else if (SpecPass == Spec2) {
    jacStart = nprocs;
    jacEnd = n;
    nBcasts = min(jacEnd-jacStart+1, nprocs);
  }
  else {
    jacStart = 1;
    jacEnd = n;
    nBcasts = min(n, nprocs);
    if (SpecPass != NoSpec) {
    cerr << "LSQNLF::LSQFDJac: Invalid speculative Jacobian option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  for (i=me+jacStart-1; i<=jacEnd-1; i+=nprocs) {

    hieps = sqrt(max(mcheps,fcn_accrcy(i)));
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp          = xcurrent(i);
    xcurrent(i)   = xtmp + hi;

    fcn0_v(n, xcurrent, fplus, result, vptr);
#ifdef OPTPP_HAVE_MPI
    if (SpecPass == Spec1) {
      MPI_Bcast(tmpF, lsqterms_, MPI_DOUBLE, nprocs-1, MPI_COMM_WORLD);
      for (j=0; j<lsqterms_; j++)
	fx(j) = tmpF[j];
    }
#endif
    // jac.Column(i) = (fplus - fx) / hi;
    
    for(int j=0; j<jac.numRows(); j++)
      {jac(j,i)=(fplus(j)-fx(j))/hi;}
    xcurrent(i)   = xtmp;
  }

  // Share my piece of the Jacobian with everyone else, and
  // incorporate their pieces.

  if (nprocs > 1) {

    for (i=0; i<nBcasts; i++) {

      for (j=me+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  tmpJacPlus[k+lsqterms_*((j-me-jacStart)/nprocs)] = jac(k,j-1);
      }

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpJacPlus, tmpSize*lsqterms_, MPI_DOUBLE, i, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=nprocs) {
	for (k=0; k<lsqterms_; k++)
	  jac(k,j-1) = tmpJacPlus[k+lsqterms_*((j-i-jacStart)/nprocs)];
      }
    }
  }

  if (tmpJacPlus != NULL)
    delete[] tmpJacPlus;
  if (tmpF != NULL)
    delete[] tmpF;

  return jac;
}

// Compute Jacobian of function vector using central differences
SerialDenseMatrix<int,double> LSQNLF::LSQCDJac(const SerialDenseVector<int,double>& sx, const SerialDenseVector<int,double>& xc,
	                SerialDenseVector<int,double>& fx, SerialDenseMatrix<int,double>& jac)
{
  int i, j, k, tmpSize, jacStart, jacEnd, myStart, inc, nBcasts;
  double xtmp, hi, hieps;
  SerialDenseVector<int,double> fplus(lsqterms_), fminus(lsqterms_); 

  int me = 0;
  int nprocs = 1;
  int n = getDim(), result = 0;
  SerialDenseVector<int,double> fcn_accrcy(getFcnAccrcy().length());
  fcn_accrcy = getFcnAccrcy();
  SerialDenseVector<int,double> xcurrent(xc.length());
  xcurrent   = xc;
  double mcheps = DBL_EPSILON;
  SpecOption SpecPass = getSpecOption();
#ifdef OPTPP_HAVE_MPI

  int error, resultlen, flag;
  char buffer[MPI_MAX_ERROR_STRING];

  // Check to see if MPI has been initialized.

  error = MPI_Initialized(&flag);
  if (error != MPI_SUCCESS) {
    MPI_Error_string(error, buffer, &resultlen);
    cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
  }

  // If it has, reset me and nprocs accordingly.

  if (flag == 1) {
    error = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (error != MPI_SUCCESS) {
      MPI_Error_string(error, buffer, &resultlen);
      cerr << "LSQNLF::LSQCDJac: MPI Error - " << buffer << endl;
    }
  }

#endif

  // Set loop endpoints, f, and x according to which pass of
  // speculative Jacobian evaluation this is.

  if (SpecPass == Spec1) {
    if (me == nprocs-1) {
      fcn0_v(n, xcurrent, fx, result, vptr);
    }
    jacStart = 1;
    jacEnd = min(n, (int) floor((double) (nprocs-1)/2));
    if (nprocs > 1)
      inc = (int) floor((double) (nprocs-1)/2);
    else
      inc = 1;
    nBcasts = min(n, (int) floor((double) (nprocs-1)/2));
  }
  else if (SpecPass == Spec2) {
    jacStart = (int) ceil((double) nprocs/2);
    jacEnd = n;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(jacEnd-jacStart+1, (int) floor((double) nprocs/2));
  }
  else {
    jacStart = 1;
    jacEnd = n;
    if (nprocs > 1)
      inc = (int) floor((double) nprocs/2);
    else
      inc = 1;
    nBcasts = min(n, (int) floor((double) nprocs/2));
    if (SpecPass != NoSpec) {
    cerr << "LSQNLF::LSQCDJac: Invalid speculative Jacobian option - "
	 << "SpecFlag = " << SpecPass << "\n"
	 << "Assuming NoSpec..." << endl;
    }
  }

  // Compute only my piece of the Jacobian.

  myStart = (int) floor((double) me/2) + jacStart;

  for (i=myStart-1; i<=jacEnd-1; i+=inc) {

    hieps = max(mcheps,fcn_accrcy(i) );
    hieps = pow(hieps,0.333333);
    hi    = hieps*max(fabs(xcurrent(i)),sx(i));
    hi    = copysign(hi,xcurrent(i));

    xtmp  = xcurrent(i);

#ifdef OPTPP_HAVE_MPI
    if (nprocs > 1) {

      // For multiple processors, even processors look forward, and
      // odd look backward.

      if (me%2 == 0)
	xcurrent(i) = xtmp + hi;
      else
	xcurrent(i) = xtmp - hi;

      fcn0_v(n, xcurrent, fplus, result, vptr);
      // jac.Column(i) = fplus/(2*hi);
    
    for(int j=0; j<jac.numRows(); j++)
      {jac(j,i)=fplus(j)/(2*hi);}
    }
    else {
      // Otherwise, do the same as in the serial case.
#endif

      xcurrent(i) = xtmp + hi;
      fcn0_v(n, xcurrent, fplus, result, vptr);

      xcurrent(i) = xtmp - hi;
      fcn0_v(n, xcurrent, fminus, result, vptr);

      // jac.Column(i) = (fplus - fminus) / (2*hi);
   
    for(int j=0; j<jac.numRows(); j++)
      {jac(j,i)=(fplus(j)-fminus(j))/(2*hi);}

#ifdef OPTPP_HAVE_MPI
    }
#endif
    xcurrent(i) = xtmp;
  }

  if (nprocs > 1) {

    if (nprocs%2 == 0)
      tmpSize = (int) ceil((double) (2*n)/nprocs);
    else
      // If there are an odd number of processors, the last one doesn't
      // count.
      tmpSize = (int) ceil((double) (2*n)/(nprocs-1));

    double *tmpJacPlus = new double[tmpSize*lsqterms_];
    double *tmpJacMinus = new double[tmpSize*lsqterms_];

    for (i=0; i<nBcasts; i++) {

      for (j=myStart; j<=jacEnd; j+=inc) {
	for (k=0; k<lsqterms_; k++) {
	  if (me%2 == 0)
	    tmpJacPlus[k+lsqterms_*((j-myStart)/inc)] = jac(k,j-1);
	  else
	    tmpJacMinus[k+lsqterms_*((j-myStart)/inc)] = jac(k,j-1);
	}
      }

#ifdef OPTPP_HAVE_MPI
      MPI_Bcast(tmpJacPlus, tmpSize*lsqterms_, MPI_DOUBLE, 2*i, MPI_COMM_WORLD);
      MPI_Bcast(tmpJacMinus, tmpSize*lsqterms_, MPI_DOUBLE, (2*i)+1, MPI_COMM_WORLD);
#endif

      for (j=i+jacStart; j<=jacEnd; j+=inc) {
	for (k=0; k<lsqterms_; k++)
	  jac(k,j-1) = tmpJacPlus[k+lsqterms_*((j-i-jacStart)/inc)] - 
	               tmpJacMinus[k+lsqterms_*((j-i-jacStart)/inc)];
      }
    }

    if (tmpJacPlus != NULL)
      delete[] tmpJacPlus;
    if (tmpJacMinus != NULL)
      delete[] tmpJacMinus;
  }

  return jac;
}

} // namespace OPTPP
