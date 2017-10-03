//JWG

//Figure out how to deal with Norm2 command.

/*********************************************************************
*   File:  dogleg.c
*
*   The source code in this file computes a dogleg step using
*   Powell's method.
*
*   Code development:
*      12 Oct 94 - Originated by T. Plantenga, Sandia National Labs.
*       4 Nov 94 - Converted to C++ by J.C. Meza, Sandia National Labs.
**********************************************************************/
#include "Opt.h"
#include "ioformat.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialSymDenseMatrix;

namespace OPTPP {

int dogleg(NLP1* nlp, std::ostream *fout, 
	   SerialSymDenseMatrix<int,double>& Hessian, SerialDenseVector<int,double>& grad, 
	   SerialDenseVector<int,double>& dogleg_step, SerialDenseVector<int,double>& sx, 
	   double& dnorm, double& TR_size, double step_max)       
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   This routine solves the trust region subproblem by Powell's
 *   dogleg method.  It constructs the Cauchy and Newton steps
 *      d_cp = -alpha grad
 *      d_nw = -H grad
 *   using L-BFGS Hessian approximations.  The dogleg path runs from
 *   dogleg_step = 0 to d_cp to d_nw, and its intersection with the
 *   spherical trust region constraint is the dogleg solution point.
 *   The Newton step should always be further away than the Cauchy
 *   point, but it may lie inside the trust region.
 *
 *   Parameters
 *     nlp          -->  pointer to nonlinear problem object 
 *
 *     Hessian      -->  Hessian Matrix
 *
 *     grad         -->  Vector of length n which specifies the 
 *                       gradient on input
 *     sx           -->  Scaling vector for x
 *                       gradient on input
 *     dogleg_step <-->  Vector of length n which specifies the 
 *                       newton direction on input. On output it will
 *                       contain the dogleg step 
 *     dnorm       <--   Norm of the dogleg step
 *     TR_size      -->  Trust region size
 *     step_max     -->  Maximum step size allowed
 *     debug        -->  Debug flag
 *     outfile      -->  File pointer for output messages
 *
*********************************************************************/
{
// Local variables
  int   i;
  double  gBg, alpha;
  double  norm_cp, norm_nw, gnorm;
  double  a, b, c, tau;
  int n = nlp->getDim();
  SerialDenseVector<int,double> tmp_vec(n);
  SerialDenseVector<int,double> scaled_grad(n);
  double CP_length;
  int trace = 0;
  int  beta = Hessian.numRows();
  int chi = Hessian.numCols();
  SerialDenseMatrix<int,double> Hessian2(beta, chi);
  SerialDenseVector<int,double> sg(n);
// Return the Newton step if it's inside the trust region. 
  // norm_nw = Norm2(dogleg_step);
  norm_nw = sqrt(dogleg_step.dot(dogleg_step));
  
  if (norm_nw <= TR_size) {
    dnorm = norm_nw;
    return(Newton_Step);
  }

  if (trace) *fout << "\nNewton step outside of trust region\n";
//
// Compute Scaled gradient
//
  for (i=0; i<n; i++) scaled_grad(i) = sx(i)*grad(i);


// Compute the Cauchy point in tmp_vec.  If it lies outside
// the trust region, then use take a step to the boundary


  for (i=0; i<beta; i++)
    for (int j=0; j<chi; j++)
      Hessian2(i,j) = Hessian(i,j)*sx(i);
  for (i=0; i<n; i++)
    sg(i) = sx(i)*scaled_grad(i);
  tmp_vec.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,1.0, Hessian2, sg, 0.0);
  gBg = scaled_grad.dot(tmp_vec);
  //gnorm = Norm2(scaled_grad);
  gnorm = sqrt(scaled_grad.dot(scaled_grad));
  CP_length = (gnorm * gnorm * gnorm) / gBg;

// If TR_size has not been set
// set it equal to the minimum of CP length and max step allowed

  if (TR_size == 0.0) TR_size = min(CP_length,step_max);

// Check to make sure that the CP is inside the trust region
// If it's outside, then just take the scaled gradient to the TR boundary

  if (trace) 
    *fout << "Compute Cauchy Point\n" 
         << "gBg       = " << e(gBg,12,4) << "\n"
         << "gnorm     = " << e(gnorm,12,4) << "\n"
         << "CP_length = " << e(CP_length,12,4) << "\n";

  if (CP_length >= TR_size) {
    alpha = -TR_size / gnorm;
    dogleg_step = scaled_grad;
    dogleg_step *= alpha;
    dnorm = sqrt(dogleg_step.dot(dogleg_step));
    //dnorm = Norm2(dogleg_step);
    if (trace) 
      *fout  <<"Cauchy point outside trust region\n"
            <<"alpha = " << e(alpha,12,4) << "\n"
            <<"dnorm = " << e(dnorm,12,4) << "\n";
    return(Cauchy_Step);
  } 

// O.K.,  the CP is inside the trust region
//
  alpha = -(gnorm * gnorm / gBg);

  tmp_vec = scaled_grad;
  tmp_vec *= alpha;

  norm_cp = sqrt(tmp_vec.dot(tmp_vec));

  // norm_cp = Norm2(tmp_vec);

// The final case is to return the intersection of the segment
// connecting the Cauchy point to the Newton point with the
// trust region boundary; i.e., find tau such that
//     || d_cp + tau*(d_nw - d_cp) ||_2 = TR_size.
//  This requires solving a single quadratic equation.

  dogleg_step -= tmp_vec;

  a = dogleg_step.dot(dogleg_step);
  b = 2.0 *(dogleg_step.dot(tmp_vec));
  c = (TR_size * TR_size) - (norm_cp * norm_cp);
  tau = (-b + sqrt(b*b + 4.0*a*c)) / (2.0*a);
  dogleg_step *= tau;
  dogleg_step += tmp_vec;
  dnorm = sqrt(dogleg_step.dot(dogleg_step));

  // dnorm = Norm2(dogleg_step);

  if (trace) *fout << "Taking a dogleg step\n"
                  << "dnorm = " << e(dnorm,12,4) << "\n";
  return(Dogleg_Step);
}

} // namespace OPTPP

