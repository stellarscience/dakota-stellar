#ifndef OptNIPSLike_h
#define OptNIPSLike_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 P.J. Williams, Sandia National Laboratories, pwillia@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifndef OptConstrNewtonLike_h
#include "OptConstrNewtonLike.h"
#endif

namespace OPTPP {

/**
 * @class OptNIPSLike 
 * OptNIPSLike is a derived class of OptConstrNewtonLike.
 * OptNIPSLike provides common data and functionality for OptFDNIPS,
 * OptQNIPS, and OptNIPS.
 *
 * The OptNIPS algorithm is a C++ implementation of NIPSM, a nonlinear
 * interior-point code developed under MATLAB by Amr El-Bakry at Rice
 * University and NIPSF, a Fortran implementation of the same code written
 * by Frederik Saaf.  Additional features include the merit functions 
 * proposed by Miguel Argaez and Richard Tapia in "Global Convergence of a
 * Primal-Dual Newton Interior-Point Method for Nonlinear Programming
 * Using a Modified Augmented Lagrange Function" as well as 
 * Robert Vanderbei and David Shanno in "An Interior-Point Algorithm For
 * Nonconvex Nonlinear Programming".
 *
 * @author P.J. Williams, Sandia National Laboratories, pwillia@sandia.gov
 */

class OptNIPSLike: public OptConstrNewtonLike {
 protected:
  virtual NLP1* nlprob() const= 0; ///< pointer to NLP1
  double          beta_;  ///< penalty parameter for merit function 3
  double          dirder_; ///< directional derivative of a merit function
  double          mu_; 	///< pertubation parameter 
  double          penalty_; ///< penalty parameter for merit function 2
  double		sigmin_; ///< centering parameter 
  double		taumin_; ///< percentage of steplength to boundary
  const double	rho_;    ///< constant set to .5 
  const double	sw_;	///<  constant

 public:
 /**
  * Default Constructor
  * @see OptNIPSLike(int n)
  * @see OptNIPSLike(int n, UPDATEFCN u)
  * @see OptNIPSLike(int n, TOLS t)
  */
  OptNIPSLike(): OptConstrNewtonLike(), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(0.0e0),
    sigmin_(0.0e0), taumin_(0.0e0), rho_(0.0e0), sw_(0.0e0)
    {strcpy(method,"Nonlinear Interior-Point Method");}
 /**
  * @param n an integer argument.
  */
  OptNIPSLike(int n): OptConstrNewtonLike(n), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param n an integer argument.
  * @param u a function pointer.
  */
  OptNIPSLike(int n, UPDATEFCN u): OptConstrNewtonLike(n,u), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }
 /**
  * @param n an integer argument.
  * @param t tolerance class reference.
  */
  OptNIPSLike(int n, TOLS t): OptConstrNewtonLike(n,t), beta_(0.0e0), 
    dirder_(0.0e0), mu_(0.0e0), penalty_(1.0e2),
    sigmin_(0.1e0), taumin_(0.95e0), rho_(5.0e-1), sw_(1.0e2)
    { strcpy(method,"Nonlinear Interior-Point Method"); }

 /**
  * Destructor
  */
  virtual ~OptNIPSLike(){}

//-------------------------------------------------------------------
// Virtual Functions 
//-------------------------------------------------------------------

  virtual void setMeritFcn(MeritFcn option);
  virtual Teuchos::SerialDenseMatrix<int,double> setupMatrix(const Teuchos::SerialDenseVector<int,double>& xc);
  virtual double merit(int flag, const Teuchos::SerialDenseVector<int,double>& xc, const Teuchos::SerialDenseVector<int,double>& yc,
                     Teuchos::SerialDenseVector<int,double>& zc, Teuchos::SerialDenseVector<int,double>& sc);
  virtual Teuchos::SerialDenseVector<int,double> setupRHS(const Teuchos::SerialDenseVector<int,double>& xc, double mu);
  virtual Teuchos::SerialDenseVector<int,double> setupRHS(const Teuchos::SerialDenseVector<int,double>& xplus,
                                const Teuchos::SerialDenseVector<int,double>& yplus,
				const Teuchos::SerialDenseVector<int,double>& zplus,
				const Teuchos::SerialDenseVector<int,double>& splus, double mu);

  virtual Teuchos::SerialSymDenseMatrix<int,double> updateH(Teuchos::SerialSymDenseMatrix<int,double>& H, int k) = 0;

  virtual int checkConvg();
  virtual int checkDeriv();
  virtual int computeStep(Teuchos::SerialDenseVector<int,double> step);

  virtual void initOpt();		///< Initialize algorithmic parameters
  virtual void initHessian();		///< Initialize Hessian of Lagrangian 
  virtual void optimize();		///< Call the interior-point method
  virtual void printStatus(char *s);	///< Print status of opt. method
  virtual void readOptInput();		///< Read user-specified input options

//-------------------------------------------------------------------
// Accessor Methods
//-------------------------------------------------------------------

/**
 * @return The value of mu_, the pertubation parameter
 */
  double getMu() const { return mu_;}
/**
 * Set the value of the perturbation parameter
 */
  void setMu(double newMu) { mu_ = newMu;}

/**
 * Sets the value of the centering parameter.
 */
  void setCenteringParameter(double newSigma) { sigmin_ = newSigma;}

/**
 * Sets the percentage of step taken towards the boundary 
 */
  void setStepLengthToBdry(double newTau) { taumin_ = newTau;}

//-------------------------------------------------------------------
// These are used by the derived classes 
//-------------------------------------------------------------------

 /**
  * Takes zero arguments with void return.
  * Resets parameter values.
  */
  virtual void reset();

  void recoverFeasibility(Teuchos::SerialDenseVector<int,double> xinit, CompoundConstraint* constraints, 
                          double ftol);
  Teuchos::SerialDenseVector<int,double> computeSearch2(Teuchos::SerialDenseMatrix<int,double>& Jacobian, Teuchos::SerialDenseVector<int,double>& rhs);
  /**
   * Takes two arguments and returns a Teuchos::SerialDenseVector<int,double>.
   * @param df a Teuchos::SerialDenseVector<int,double> - gradient of obj. function
   * @param dcon a Matrix  - gradient of constraints
   * @return The initial value of the Lagrange multiplier z 
   */
  Teuchos::SerialDenseVector<int,double> initMultipliers(const Teuchos::SerialDenseVector<int,double>& df, Teuchos::SerialDenseMatrix<int,double>& dcon);

  /**
   * Takes one arguments and updates the perturbation parameter
   * @param k an integer - iteration counter 
   */
  void updateMu(int k);

  /**
   * Takes five arguments and returns a double value.
   * @param flag an integer argument
   * @param xc a Teuchos::SerialDenseVector<int,double> 
   * @param yc a Teuchos::SerialDenseVector<int,double> of Lagrange multipliers 
   * @param zc a Teuchos::SerialDenseVector<int,double> of Lagrange multipliers 
   * @param sc a Teuchos::SerialDenseVector<int,double> of slack variables
   * @see merit(flag,xc,yc,zc,sc)
   * @return The value of the Argaez-Tapia merit function.
   */
  double merit2(int flag, const Teuchos::SerialDenseVector<int,double>& xc, const Teuchos::SerialDenseVector<int,double>& yc,
              Teuchos::SerialDenseVector<int,double>& zc, Teuchos::SerialDenseVector<int,double>& sc);
  /**
   * Takes four arguments and returns a double value.
   * @param flag an integer argument
   * @param xc a Teuchos::SerialDenseVector<int,double> 
   * @param zc a Teuchos::SerialDenseVector<int,double> of Lagrange multipliers 
   * @param sc a Teuchos::SerialDenseVector<int,double> of slack variables
   * @see merit(flag,xc,yc,zc,sc)
   * @return The value of the Vanderbei et al merit function.
   */
  double merit3(int flag, const Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& zc,
              Teuchos::SerialDenseVector<int,double>& sc);

  /**
   * Takes three arguments and void return.
   * @param sk a Teuchos::SerialDenseVector<int,double> with contains the search direction 
   * @param xc a Teuchos::SerialDenseVector<int,double> of current point 
   * @param derivative a Teuchos::SerialDenseVector<int,double> of derivative of cost function 
   */
  void computeDirDeriv(Teuchos::SerialDenseVector<int,double>& sk, const Teuchos::SerialDenseVector<int,double>& xc,
                       Teuchos::SerialDenseVector<int,double>& derivative);
  /**
   * Takes one arguments and returns a real value.
   * @param step a Teuchos::SerialDenseVector<int,double> with contains the search direction 
   */
  double dampenStep(Teuchos::SerialDenseVector<int,double>& step);

};

} // namespace OPTPP

#endif
