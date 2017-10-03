#ifndef optbcnewton2_h
#define optbcnewton2_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "OptBCNewtonLike.h"


namespace OPTPP {

/**
 * OptBaNewton implements a bound constrained Newton method with a 
 * logarithmic barrier term. Barrier methods transform a constrained into an
 * unconstrained problem via a barrier term.  In this implementation, the 
 * unconstrained subproblem's objective function consists of the original 
 * objective function and a logarithmic barrier term that prevents iterates 
 * from becoming infeasible.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 11/2005
 */

class OptBaNewton: public OptBCNewton2Deriv {

  /// Barrier parameter
  double          mu;
  /// Value of barrier function
  double          fvalue_barrier, fprev_barrier, fprev_outer;
  /// Gradient of barrier function
  Teuchos::SerialDenseVector<int,double>    grad_barrier, gprev_barrier; 
  /// Hessian of barrier function 
  Teuchos::SerialSymDenseMatrix<int,double> Hess_barrier; 

 public:

// Analytic Hessian case

 /**
   * @param p a pointer to an NLP2 object.
   * @see OptNIPSLike(NLP2* p, UPDATEFCN u)
   * @see OptNIPSLike(NLP2* p, TOLS t)
   */
  OptBaNewton(NLP2* p): OptBCNewton2Deriv(p)
    {strcpy(method,"Bound constrained Newton with barrier");}

 /**
   * @param p a pointer to an NLP2 object.
   * @param u a function pointer.
   * @see OptNIPSLike(NLP2* p)
   * @see OptNIPSLike(NLP2* p, TOLS t)
   */
  OptBaNewton(NLP2* p, UPDATEFCN u): OptBCNewton2Deriv(p, u)
    {strcpy(method,"Bound constrained Newton with barrier");}

 /**
   * @param p a pointer to an NLP2 object.
   * @param t tolerance class reference.
   * @see OptNIPSLike(NLP2* p)
   * @see OptNIPSLike(NLP2* p, UPDATEFCN u)
   */
  OptBaNewton(NLP2* p, TOLS t): OptBCNewton2Deriv(p, t)
    {strcpy(method,"Bound constrained Newton with barrier");}

 /**
  * Destructor 
  */
  virtual ~OptBaNewton(){ }

//--------------------------------------
// These are defined elsewhere
//--------------------------------------
  void            initOpt();
  void            initHessian();
  int             checkInnerConvg(int);
  int             checkConvg();
  void            acceptStep(int,int);
  void            optimize();
  void            updateBarrierMultiplier();
  Teuchos::SerialSymDenseMatrix<int,double> updateH(Teuchos::SerialSymDenseMatrix<int,double>& H, int k);
  double          compute_Barrier_Fvalue(double,Teuchos::SerialDenseVector<int,double>&);
  Teuchos::SerialDenseVector<int,double>    compute_Barrier_Gradient(Teuchos::SerialDenseVector<int,double>&,
    Teuchos::SerialDenseVector<int,double>&);
  Teuchos::SerialSymDenseMatrix<int,double> compute_Barrier_Hessian(Teuchos::SerialSymDenseMatrix<int,double>&,
    Teuchos::SerialDenseVector<int,double>&);
  Teuchos::SerialDenseVector<int,double>    computeSearch2(Teuchos::SerialSymDenseMatrix<int,double>&, 
    Teuchos::SerialDenseVector<int,double>&);
  int             computeStep(Teuchos::SerialDenseVector<int,double> );
  double          computeMaxStep(Teuchos::SerialDenseVector<int,double> &);
  double          scalarNewton(double,double,double,double,double);
  void            setAsideCurrentVariables();

};
} // namespace OPTPP
#endif
