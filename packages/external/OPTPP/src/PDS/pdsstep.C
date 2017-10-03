//JWG

/*********************************************************************
 *
 * File:  pdsstep.c
 *
 * Purpose:
 *    Computes a trust region step using PDS
 *
 *********************************************************************/


#include "OptPDS.h"
#include "NLF.h"
#include "ioformat.h"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

using std::streambuf;
using std::endl;
using std::flush;

extern "C" {
int pdsdgn(int ndim, double *s, double *a, double *work1, 
	    double *qraux, int *jpvt, double *rcond);
}

namespace OPTPP {

int pdsstep(NLP1* nlp, std::ostream *fout, SerialSymDenseMatrix<int,double>& Hessian,
	    SerialDenseVector<int,double>& grad, SerialDenseVector<int,double>& pds_step,
	    SerialDenseVector<int,double>& sx, double& step_norm, double& TR_size,
	    double step_max, double& radius, bool init_value,
	    int sss)
{
  /*******************************************************************
   *
   * This routine solves the trust region subproblem by using PDS
   *
   * Parameters
   *
   *    nlp          -->  pointer to nonlinear problem object 
   *
   *    Hessian      -->  Hessian Matrix
   *
   *    grad         -->  Vector of length n which specifies the
   *                      gradient on input
   *
   *    sx           -->  Scaling vector for x gradient on input
   *
   *    pds_step    <-->  Vector of length n which specifies the 
   *                      newton direction on input. On output it
   *                      will contain the step 
   *
   *    step_norm   <--   Norm of the step taken
   *
   *    TR_size      -->  Trust region size
   *
   *    step_max     -->  Maximum step size allowed
   *
   *    debug        -->  debug flag
   *
   *    outfile      -->  File pointer for output messages
   *
   *******************************************************************/

// Local variables

  int   i;
  double  gBg, lambda;
  double  current_val, cauchy_val, pds_val, newton_val;
  double  c_decrease, n_decrease, check_norm, check_norm1;
  double  norm_cauchy, norm_newton, grad_norm;
  double  a1, b, c, tau;
  int ndim = nlp->getDim();
  SerialDenseVector<int,double> cauchy_point(ndim), cauchy_step(ndim);
  SerialDenseVector<int,double> newton_point(ndim), newton_step(ndim);
  SerialDenseVector<int,double> scaled_grad(ndim), diff(ndim), diff1(ndim);
  SerialDenseVector<int,double> xc(ndim);
  SerialDenseVector<int,double> new_x(ndim);

  double CP_length;
  bool dogleg;
  bool debug = nlp->getDebug();

  static bool first_time = true;

  // Use a TR 2 times bigger than the original trust region in order
  // to allow for twice the Newton direction.

  double PDS_TR_size = 2.0*TR_size;

  //
  //  set up the problem
  //

  xc           = nlp->getXc();
  current_val  = nlp->getF();
  newton_step  = pds_step;
  newton_point = xc;
  newton_point += newton_step;
  // norm_newton  = Norm2(newton_step);
  norm_newton = sqrt(newton_step.dot(newton_step));
  //
  // Compute Scaled gradient
  //
  int n = sx.length();
   for (i=0; i<ndim; i++) scaled_grad(i) = sx(i)*grad(i);
  
  //
  // Compute the Cauchy point in cauchy_point.  If it lies outside
  // the PDS trust region, then take a step to the boundary
int  beta = Hessian.numRows();
  int chi = Hessian.numCols();
  SerialDenseMatrix<int,double> Hessian2(beta, chi);
  SerialDenseVector<int,double> sg(ndim);

for (i=0; i<beta; i++)
    for (int j=0; j<chi; j++)
      Hessian2(i,j) = Hessian(i,j)*sx(i);
  for (i=0; i<ndim; i++)
    sg(i) = sx(i)*scaled_grad(i);
  cauchy_step.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,1.0, Hessian2, sg, 0.0);
  
  //cauchy_step = Dx*Hessian*Dx*scaled_grad;
  gBg         = scaled_grad.dot(cauchy_step);
  //grad_norm   = Norm2(scaled_grad);
  grad_norm = sqrt(scaled_grad.dot(scaled_grad));
  CP_length   = (grad_norm * grad_norm * grad_norm) / gBg;

  if (debug) 
    *fout << "Compute Cauchy Point\n" 
          << "gBg       = " << e(gBg,12,4) << "\n"
          << "grad_norm = " << e(grad_norm,12,4) << "\n"
          << "CP_length = " << e(CP_length,12,4) << "\n";

  // If PDS_TR_size has not been set
  // set it equal to the minimum of CP length and max step allowed
  
  if (PDS_TR_size == 0.0) PDS_TR_size = min(2*CP_length,step_max);

  // Check to make sure that the CP is inside the trust region
  // If it's outside, then just take the scaled gradient to the TR
  // boundary.  In this case, can't use the dogleg point if needed.

  if (CP_length >= TR_size) {
    dogleg = false;
    lambda    = TR_size / grad_norm;
    if (debug) 
      *fout <<"Cauchy point outside trust region\n"
            <<"lambda    = " << e(lambda,12,4) << "\n";
  } 

  //
  // O.K., the CP is inside the trust region - can use dogleg.
  //

  else {
    dogleg = true;
    lambda       = (grad_norm * grad_norm / gBg);
  }

		   cauchy_step  = scaled_grad;
		   cauchy_step *= -lambda;
		   cauchy_point = xc;
		   cauchy_point += cauchy_step;
		   // norm_cauchy  = Norm2(cauchy_step);
		   norm_cauchy = sqrt(cauchy_step.dot(cauchy_step));

  if (norm_newton > TR_size) {

    if (debug) 
      *fout <<"Newton point outside trust region\n";

    // If the Newton point is outside the trust region, compute the
    // dogleg point if allowable.  Otherwise, take the scaled Newton
    // direction to the TR boundary.

    if (dogleg) {
      (*fout) << "using DOGLEG\n";
      newton_step -= cauchy_step;
      a1 = newton_step.dot(newton_step);
      b = 2.0 * (newton_step.dot(cauchy_step));
      c = (TR_size * TR_size) - (norm_cauchy * norm_cauchy);
      tau = (-b + sqrt(b*b + 4.0*a1*c)) / (2.0*a1);
      newton_step *= tau;
      newton_step += cauchy_step;
    }
    else {
      (*fout) << "PROJECT Newton onto TR\n";
      tau = (TR_size/norm_newton);
      newton_step *= tau;
    }

    // norm_newton = Norm2(newton_step);
    norm_newton = sqrt(newton_step.dot(newton_step));
		       newton_point = xc;
		       newton_point += newton_step;
  }

  (*fout) << "|| S_CP || =" << e(norm_cauchy,12,4) << endl
	  << "|| S_N  || =" << e(norm_newton,12,4) << endl
	  << " TR SIZE   =" << e(TR_size,12,4) << endl;

  //
  // set up the PDS problem
  // a must be dimensioned ndim*ndim
  // work1, qraux, jpvt must be dimensioned ndim

  int j;
  int *jpvt = new int[ndim];
  double *a = new double[ndim*ndim];
  double *work1 = new double[ndim];
  double *qraux = new double[ndim];
  double rcond, rcond_tol, ratio;
  double cauchy_angle, newton_angle;

  SerialDenseVector<int,double> vscale(ndim);
  SerialDenseMatrix<int,double> init_simplex(ndim,ndim+1);
  const char *schemefilename = {"myscheme"};

  double *stemparray = new double[(ndim+1)*ndim];


  NLP1 *nlp_subproblem;
  nlp_subproblem = nlp;

  OptPDS subproblem(nlp_subproblem);      
  *(fout) << flush;
  streambuf* save_buffer = fout->rdbuf();
  subproblem.setOutputFile(*fout);
  fout->rdbuf(save_buffer);

   //Initialize max_proc and search_scheme_size
   //in case there is no opt.input file
  

   int max_proc = 64;
  int search_scheme_size = 4*ndim;

  if (search_scheme_size > max_proc)
  search_scheme_size = max_proc;

  subproblem.setSSS(sss);

  // set default values 
  //

  subproblem.setFcnTol(0.9);
  subproblem.setMaxIter(10);

  for (i=0; i < ndim; i++) {
    diff(i) = newton_point(i) - cauchy_point(i);
  }

//   If this is the first nonlinear iteration, Cauchy and Newton
//   directions will be the same (assuming initial Hessian is scaled
//   identity).  Use the Cauchy point and the current point as
//   vertices in the initial simplex.  get remaining vertices from
//   right-angle simplex constructed around the Cauchy point.  Edges
//   are set to be the original trust region size.

  if (init_value) {

    for (i=0; i < ndim; i++) {
      init_simplex(i,0) = cauchy_point(i);
      init_simplex(i,1) = xc(i);
    }
  
    for (j=2; j < ndim+1; j++) {

      for (i=0; i < ndim; i++) {
	init_simplex(i,j) = init_simplex(i,0);
      }
  
      init_simplex(j-1,j) = init_simplex(j-1,j) + TR_size;
  
    }
  }

  // Otherwise, use the Cauchy, Newton, and current points as
  // vertices.  Obtain the remaining points from a right-angle simplex
  // constructed around the Newton point.  Edges are set to be the
  // original trust region size.

  else {

    for (i=0; i < ndim; i++) {
      init_simplex(i,0) = cauchy_point(i);
      init_simplex(i,1) = newton_point(i);
      init_simplex(i,2) = xc(i);
    }

    for (j=3; j < ndim+1; j++) {

      for (i=0; i < ndim; i++) {
	init_simplex(i,j) = init_simplex(i,1);
      }

      init_simplex(j-1,j) = init_simplex(j-1,j) + TR_size;
    }
  }

  // Check to see if the initial simplex is degenerate.  The choice of
  // tolerance is somewhat ad hoc - it is a problem-dependent
  // quantity that demonstrated a correspondence to the acceptable
  // condition number of the simplex, and it worked well in
  // experiments.
  int count = 0;
  rcond_tol = min(10*gBg, 10/gBg);
  for(i=0;i<init_simplex.numCols();i++)
    for(j=0;j<init_simplex.numRows();j++)
      {stemparray[count] = init_simplex(j,i);
	count = count+1;}



  // pdsdgn(ndim, stemp.Store(), a, work1, qraux, jpvt, &rcond);
pdsdgn(ndim, stemparray, a, work1, qraux, jpvt, &rcond);

  if (debug) (*fout) << "\npdsstep: initial rcond = " << rcond
		     << "\n";

  // If the original simplex is not degenerate, use it.
 
  if (rcond > rcond_tol) {
    (*fout) << "\npdsstep: Using Newton Simplex\n";
    subproblem.setSimplexType(4);
  }

  // Otherwise, scale the Cauchy edge and the right-angle edges to be
  // the same length as the Newton edge.

  else {
    ratio = norm_cauchy/norm_newton;
    (*fout) << "ratio =" << e(ratio,12,4) << endl;

    for (i=0; i < ndim; i++)
      cauchy_step(i) = (1/ratio)*cauchy_step(i);

    cauchy_point = xc;
    cauchy_point += cauchy_step;
    //norm_cauchy = Norm2(cauchy_step);
    norm_cauchy = sqrt(cauchy_step.dot(cauchy_step));
    (*fout) << "\npdsstep: Newton Simplex degenerate. Using default "
	    << "simplex\n";

    if (init_value) {

      for (i=0; i < ndim; i++) {
	init_simplex(i,0) = cauchy_point(i);
	init_simplex(i,1) = xc(i);
      }

      for (j=2; j < ndim+1; j++) {

	for (i=0; i < ndim; i++) {
	  init_simplex(i,j) = init_simplex(i,0);
	}

	init_simplex(j-1,j) = init_simplex(j-1,j) + norm_newton;
      }
    }
    else {

      for (i=0; i < ndim; i++) {
	init_simplex(i,0) = cauchy_point(i);
	init_simplex(i,1) = newton_point(i);
	init_simplex(i,2) = xc(i);
      }

      for (j=3; j < ndim+1; j++) {

	for (i=0; i < ndim; i++) {
	  init_simplex(i,j) = init_simplex(i,1);
	}

	init_simplex(j-1,j) = init_simplex(j-1,j) + norm_newton;
      }
    }

    // Check to see if the new simplex has a better condition number.
    // Currently, this information is not used other than being
    // printed out.

    //    stemp = init_simplex.t();
//     int eta = init_simplex.numRows();
//     int gamma = init_simplex.numCols();
//     for(int i=0; i<gamma; i++)
//       for(int j=0; j<eta; j++)
// 	{stemp(i,j) = init_simplex(j,i);}

    //double stemparray[(ndim+1)*ndim];  
  
    // pdsdgn(ndim, stemp.Store(), a, work1, qraux, jpvt, &rcond);
   
  count = 0;
  
  for(i=0;i<init_simplex.numCols();i++)
    for(j=0;j<init_simplex.numRows();j++)
      {stemparray[count] = init_simplex(j,i);
	count = count+1;}

 pdsdgn(ndim, stemparray, a, work1, qraux, jpvt, &rcond);

    if (debug) (*fout) << "\npdsstep:  new rcond = " << rcond << "\n";
    subproblem.setSimplexType(4);
  }

  // Print out simplex vertices.

  if (debug) {

    if (init_value) 
      (*fout) << "\npdsstep: Simplex order: CP, XC, ...\n";
    else
      (*fout) << "pdsstep: Simplex order: CP, XN, XC, ...\n";

    (*fout) << "pdsstep: Initial Simplex\n";
    FPrint(fout, init_simplex);
  }

  (*fout) << flush;

  // Initialize parameters for PDS.

  vscale = 1.0;


  subproblem.setScale(vscale);
  subproblem.setSimplex(init_simplex);

  //
  //  Create scheme file on first time through
  //

  if (first_time) {
    subproblem.setCreateFlag();
    first_time = false;
  }
  else {
    subproblem.setCreateFlag(false);
  }

  subproblem.setSchemeFileName((char*)schemefilename);
  subproblem.setTRSize(PDS_TR_size);
  subproblem.setNonIter(init_value);
  subproblem.setTRPDS(true);

  (*fout) << "\n*************************";
  (*fout) << "  Begin Inner Iteration  ";
  (*fout) << "*************************\n";
  (*fout) << flush;

  subproblem.optimize();

  std::string status_string("Solution from PDS subproblem");
  if (debug) subproblem.printStatus( (char*)status_string.c_str() );
  subproblem.cleanup();

  (*fout) << "\n*************************";
  (*fout) << "   End Inner Iteration   ";
  (*fout) << "*************************\n";
  (*fout) << flush;

  //
  // What is the step we took
  //

  new_x   = nlp_subproblem->getXc();
  pds_val = nlp->getF();
  radius  = subproblem.getSimplexSize();





  pds_step = new_x;
  pds_step -=  xc;
 //  step_norm = Norm2(pds_step);
  step_norm = sqrt(pds_step.dot(pds_step));
  // If debug is set, print out the function values at each of the
  // primary points in the simplex, the fraction of Cauchy/ Newton
  // decrease obtained, and the angle between the step and the
  // Cauchy/Newton directions.  Notice that this requires extra
  // function evaluations because the values for the Cauchy and Newton
  // points are not readily available.

  if (debug) {
    (*fout) << "\npdsstep: Current point \tNew point \tPDS step"
	    << "\tNewton step\n";

    for (i=0; i<ndim; i++) {
      (*fout) << d(i,4) << e(xc(i),18,6) << e(new_x(i),14,6) 
	      << e(pds_step(i),14,6) << e(newton_step(i),14,6)
	      << endl;
    }

    SpecOption SpecTmp = nlp->getSpecOption();
    nlp->setSpecOption(NoSpec);
    cauchy_val  = nlp->evalF(cauchy_point);
    newton_val  = nlp->evalF(newton_point);
    nlp->setSpecOption(SpecTmp);

    for (i=0; i<ndim; i++)
      { diff(i) = pds_step(i);
	diff(i) -= newton_step(i);}

    // check_norm   = Norm2(diff)/Norm2(newton_step);
    check_norm = sqrt(diff.dot(diff))/sqrt(newton_step.dot(newton_step));
    c_decrease   = (current_val - pds_val)/(current_val - cauchy_val);
    n_decrease   = (current_val - pds_val)/(current_val - newton_val);
    cauchy_angle = (cauchy_step.dot(pds_step))/(norm_cauchy*step_norm);
    newton_angle = (newton_step.dot(pds_step))/(norm_newton*step_norm);

    (*fout) << "\nf(X_PDS)                 = " << e(pds_val,19,6) << endl
            << "f(X_CP)                  = " << e(cauchy_val,19,6) << endl
            << "fraction Cauchy decrease = " << e(c_decrease,19,6) << endl
	    << "angle from Cauchy dir.   = " << e(cauchy_angle,19,6) << endl
            << "f(X_N)                   = " << e(newton_val,19,6) << endl
            << "fraction Newton decrease = " << e(n_decrease,19,6) << endl
	    << "angle from Newton dir.   = " << e(newton_angle,19,6) << endl
	    << "||S_PDS - S_N||/||S_N||  = " << e(check_norm,19,6) << endl;
  }

  // Determine if the step is the unit Newton step, the unit Cauchy
  // step, or a different step.

  for (i=0; i<ndim; i++) {
    diff(i)  = new_x(i) - cauchy_point(i);
    diff1(i) = new_x(i) - newton_point(i);
  }

    // check_norm  = Norm2(diff);
    // check_norm1 = Norm2(diff1);
    check_norm = sqrt(diff.dot(diff));
    check_norm1 = sqrt(diff1.dot(diff1));

  if (jpvt != NULL)
    delete[] jpvt;
  if (a != NULL)
    delete[] a;
  if (work1 != NULL)
    delete[] work1;
  if (qraux != NULL)
    delete[] qraux;
  if (stemparray != NULL)
    delete[] stemparray;

   if (check_norm < 1.e-8)
     return(0);
   else if (check_norm1 < 1.e-8)
     return(2);
   else
    return(1);
}

} // namespace OPTPP


