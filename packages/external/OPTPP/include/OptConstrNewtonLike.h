#ifndef OptConstrNewtonLike_h
#define OptConstrNewtonLike_h

/*--------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.
  J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 -----------------------------------------------------------------------*/

#include "Opt.h"

namespace OPTPP {

/**
 * Constrained Newton abstract data classes
 * OptConstrNewtonLike
 * OptConstrNewton1Deriv
 * OptConstrNewton2Deriv
 *
 * OptConstrNewtonLike provides common data and functionality for
 * OptConstrFDNewton, OptConstrQNewton, OptConstrNewton, OptFDNIPS,
 * OptQNIPS, and OptNIPS methods.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams
 * @date 03/27/2007 
 */


class OptConstrNewtonLike: public OptimizeClass {
protected:
  virtual NLP1* nlprob() const = 0; ///< returns an NLP1 pointer
  int me;			///< Number of equality constraints
  int mi;			///< Number of inequality constraints
  int grad_evals;		///< Number of gradient evaluations 
  Teuchos::SerialDenseVector<int,double> gprev;		///< Gradient at previous iteration
  Teuchos::SerialDenseVector<int,double> z;		///< Lagrange mult. corr. to ineq constraints 
  Teuchos::SerialDenseVector<int,double> y;		///< Lagrange mult. corr. to eq constraints
  Teuchos::SerialDenseVector<int,double> s;		///< Slack variables corr. to ineq. constraints 
  Teuchos::SerialDenseVector<int,double> constrType;	///< Vector which contains the type of constraint
  Teuchos::SerialDenseVector<int,double> constraintResidual; ///< Constraint residual at xc
  Teuchos::SerialDenseVector<int,double> gradl;		///< Current gradient of the lagrangian
  Teuchos::SerialDenseVector<int,double> gradlprev;	///< Previous gradient of the lagrangian
  Teuchos::SerialDenseMatrix<int,double> constraintGradient;	///< Current constraint gradient 
  Teuchos::SerialDenseMatrix<int,double> constraintGradientPrev; ///< Previous constraint gradient
  Teuchos::SerialSymDenseMatrix<int,double> Hessian;	///< Current Hessian
  Teuchos::SerialSymDenseMatrix<int,double> hessl;	///< Current Hessian of the lagrangian
  SearchStrategy strategy;	///< User-specified globalization strategy
  DerivOption finitediff;	///< User-specified derivative option
  MeritFcn    mfcn;		///< Merit function option
  double TR_size;			///< Size of the trust region radius 
  double gradMult;		///< Gradient multiplier to compute TR_size
  int searchSize;               ///< Search pattern size for TRPDS
  double cost;			///< Value of the merit function
  void defaultAcceptStep(int, int);
  Teuchos::SerialDenseVector<int,double> defaultComputeSearch(Teuchos::SerialSymDenseMatrix<int,double>& );
  bool WarmStart;
  bool feas_flag;		///< Switch to turn on the feasibility recovery method
  int max_feas_iter;            ///< Maximize number of iterations to achieve feasibility 

public:
 /**
  * Default Constructor
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */

  OptConstrNewtonLike() {}
 /**
  * @param n an integer argument
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */
  OptConstrNewtonLike(int n): 
    OptimizeClass(n), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {z = 0; y = 0; s = 0;}
 /**
  * @param n an integer argument
  * @param u a function pointer.
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, TOLS t)
  */
  OptConstrNewtonLike(int n, UPDATEFCN u): 
    OptimizeClass(n), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {update_fcn = u; z = 0; y = 0; s = 0;}
 /**
  * @param n an integer argument
  * @param t tolerance class reference.
  * @see OptConstrNewtonLike(int n)
  * @see OptConstrNewtonLike(int n, UPDATEFCN u)
  */
  OptConstrNewtonLike(int n, TOLS t): 
    OptimizeClass(n,t), me(0), mi(0), grad_evals(0),   gprev(n), 
    z(n), y(n), s(n), constrType(n), constraintResidual(n), 
    gradl(n), gradlprev(n),constraintGradient(n,n), constraintGradientPrev(n,n),
    Hessian(n), hessl(n), strategy(TrustRegion), finitediff(ForwardDiff), 
    mfcn(ArgaezTapia), TR_size(0.0), 
    gradMult(0.1), searchSize(64), cost(0.0), WarmStart(false),
    feas_flag(false), max_feas_iter(3)
    {z = 0; y = 0; s = 0;}
  
 /**
  * Destructor
  */
  virtual ~OptConstrNewtonLike(){}

 //---------------------------------------
 // Virtual functions 
 //---------------------------------------

// These have default values
  /// Accept this step and update the nonlinear model 
  virtual void acceptStep(int k, int step_type)
    {defaultAcceptStep(k, step_type);}

  /// Solve for the Newton direction
  virtual Teuchos::SerialDenseVector<int,double> computeSearch(Teuchos::SerialSymDenseMatrix<int,double>& H )
    {return defaultComputeSearch(H);}

  virtual void updateModel(int k, int ndim, Teuchos::SerialDenseVector<int,double> x)
    {OptimizeClass::defaultUpdateModel(k, ndim, x);}

  virtual void reset();

 // These have to be defined by other classes

  /// Check to see if algorithm satisfies the convergence criterion
  virtual int  checkConvg();

  /// Compare the analytic gradient with the finite-difference approximation 
  virtual int  checkDeriv();

  /// Compute the steplength along the search direction.  
  /// If an acceptable step not found, returns an error code = -1. 
  virtual int  computeStep(Teuchos::SerialDenseVector<int,double> sk);

  virtual double computeMaxStep(Teuchos::SerialDenseVector<int,double> sk);

  /// Initialize algorithmic parameters 
  virtual void initOpt();

  /// Compute the Hessian or its approximation at the initial point
  virtual void initHessian();

  /// Initialize the radius of the trust-region 
  virtual double initTrustRegionSize() const;

  /// Invoke a constrained Newton's method 
  virtual void optimize();

  /// Read user-specified input options from a file
  virtual void readOptInput();

// These are used by all derived classes

  /// Check finite difference approximation with analytic derivatives
  int checkAnalyticFDGrad();

  /**
   * @return The number of function evaluations
   */
  int getFevals() const {return fcn_evals;}

  /**
   * @return The number of gradient evaluations
   */
  int getGevals() const {return grad_evals;}

  /**
   * @return Trust region radius 
   */
  double getTRSize() const {return TR_size;}
  /// Set radius of trust-region
  void setTRSize(double delta) {TR_size = delta;}

  /**
   * @return Gradient multiplier to compute radius of trust-region 
   */
  double getGradMult() const {return gradMult;}
  /// Set gradient multiplier which is used to compute radius of trust-region 
  void setGradMult(double tau) {gradMult = tau;}

  /**
   * @return Number of points in search scheme 
   * (only relevant when trustpds search strategy is selected)
   */
  int getSearchSize() const {return searchSize;}
  /// Set number of points in search scheme for trustpds strategy
  void setSearchSize(int sss) {searchSize = sss;}

  bool getWarmStart() const {return WarmStart;}
  void UseWarmStart(Teuchos::SerialSymDenseMatrix<int,double>& H) {Hessian = H; WarmStart = true;}

  /**
   * @return Globalization strategy for optimization algorithm 
   */
  void setSearchStrategy(SearchStrategy search) {strategy = search;}
  /// Set globalization strategy for optimization algorithm
  SearchStrategy getSearchStrategy() const {return strategy;}

  /**
   * @return Type of finite difference approximation
   */
  void setDerivOption(DerivOption d) {finitediff = d;}
  /// Set type of finite difference approximation (forward, backward, or central)
  DerivOption getDerivOption() const {return finitediff;}

  /**
   * @return Lagrangian multipliers corresponding to equality constraints 
   */
  Teuchos::SerialDenseVector<int,double> getEqualityMultiplier() const { return y;}
  /// Store the current Lagrangian multipliers corresponding to equality constraints
  void setEqualityMultiplier(const Teuchos::SerialDenseVector<int,double>& ymult) { y  = ymult;}

  /**
   * @return Lagrangian multipliers corresponding to inequality constraints 
   */
  Teuchos::SerialDenseVector<int,double> getInequalityMultiplier() const { return z;}
  /// Store the current Lagrangian multipliers corresponding to inequality constraints
  void setInequalityMultiplier(const Teuchos::SerialDenseVector<int,double>& zmult){ z = zmult;}

  /**
   * @return Slack variables associated with inequality constraints 
   */
  Teuchos::SerialDenseVector<int,double> getSlacks() const { return s;}
  /// Store the current slack vector 
  void setSlacks(const Teuchos::SerialDenseVector<int,double>& slackVar) { s = slackVar;}

  /**
   * @return Merit function 
   */
  MeritFcn getMeritFcn() const  { return mfcn;}
  /// Specify the merit function to used in step acceptance test
  virtual void setMeritFcn(MeritFcn option) { mfcn = option;}

  /**
   * @return Gradient of the Lagrangian at the current iteration 
   */
  Teuchos::SerialDenseVector<int,double> getGradL() const  { return gradl;}
  /// Store the gradient of the Lagrangian at the current iteration
  virtual void setGradL(Teuchos::SerialDenseVector<int,double> gradl_value) { gradl = gradl_value;
}

  /**
   * @return Gradient of the Lagrangian at the previous iteration 
   */
  Teuchos::SerialDenseVector<int,double> getGradLPrev() const  { return gradlprev;}
  /// Store the gradient of the Lagrangian at the previous iteration
  virtual void setGradLPrev(Teuchos::SerialDenseVector<int,double> gradl_value) { gradlprev = gradl_value;}

  /**
   * @return Residuals of the constraints at the current iteration 
   */
  Teuchos::SerialDenseVector<int,double> getConstraintResidual() const  { return constraintResidual;}
  /// Store the residuals of the constraints at the current iteration
  virtual void setConstraintResidual(const Teuchos::SerialDenseVector<int,double>& constraint_value) 
                       { constraintResidual = constraint_value;}

  /**
   * @return Gradient of the constraints at the current iteration 
   */
  Teuchos::SerialDenseMatrix<int,double> getConstraintGradient() const  { return constraintGradient;}
  /// Store the current gradients of the constraints 
  virtual void setConstraintGradient(const Teuchos::SerialDenseMatrix<int,double>& constraint_grad) 
                       { constraintGradient = constraint_grad;}

  /**
   * @return Value of merit function 
   */
  double getCost() const  { return cost;}
  /// Store current value of merit function
  void setCost(double value) { cost = value;}

  /**
   * @return Switch to turn on feasibility recovery method 
   */
  bool getFeasibilityRecoveryFlag() const {return feas_flag;}
  /// Set switch to turn on feasibility recovery method 
  void setFeasibilityRecoveryFlag(bool flag) {feas_flag = flag;}

  /**
   * @return Maximum number of iterations for feasibility recovery method 
   */
  int getMaxFeasIter() const {return max_feas_iter;}
  /// Set maximum number of iterations for feasibility recovery method of trust-region
  void setMaxFeasIter(int k) {max_feas_iter = k;}

  /**
   * @return Hessian matrix 
   */
  Teuchos::SerialSymDenseMatrix<int,double> getHessian() const {return Hessian;}
  /// Store current Hessian matrix 
  void setHessian(Teuchos::SerialSymDenseMatrix<int,double>& H) {Hessian = H;}
  
  /// Compute the Hessian of the Lagrangrian or its approximation at iteration k
  virtual Teuchos::SerialSymDenseMatrix<int,double> updateH(Teuchos::SerialSymDenseMatrix<int,double>& H, int k) = 0;

  /// Compute the active set according to Facchinei, Fischer, and Kanzow indicator 
  Teuchos::SerialDenseVector<int,double> computeFFK1Ind(const Teuchos::SerialDenseVector<int,double>& xc);
  /// Compute the active set according to Facchinei, Fischer, and Kanzow indicator 
  Teuchos::SerialDenseVector<int,double> computeFFK2Ind(const Teuchos::SerialDenseVector<int,double>& xc);
  /// Compute the Tapia indicators for the constrained problem
  Teuchos::SerialDenseVector<int,double> computeTapiaInd(const Teuchos::SerialDenseVector<int,double>& step);

  /// Print status of the constrained Newton's method
  void printStatus(char *);
  /// Output the Lagrangian multipliers to the screen 
  void printMultipliers(char *);
  /// Print the Lagrangian multipliers to a file 
  void fPrintMultipliers(std::ostream *nlpout, char *);
  /// Print second order sufficiency information to a file 
  void fPrintSecSuff(std::ostream *nlpout, Teuchos::SerialDenseVector<int,double>& info);

  friend int trustregion(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&,
			 Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&,
			 double&, double&, double stpmax, double stpmin);

  friend int trustpds(NLP1*, std::ostream*, Teuchos::SerialSymDenseMatrix<int,double>&,
		      Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&,
		      double&, double&, double stpmax, double stpmin, int);
};

/**
 * Constrained Newton classes that will accept either an NLP1 or NLP2
 */

class OptConstrNewton1Deriv: public OptConstrNewtonLike {
public:

  OptConstrNewton1Deriv() {}

  OptConstrNewton1Deriv(NLP1* p): 
    OptConstrNewtonLike(p->getDim()), mem_nlp(p) {;}
   
  OptConstrNewton1Deriv(NLP1* p, UPDATEFCN u): 
    OptConstrNewtonLike(p->getDim(),u), mem_nlp(p) {;}

  OptConstrNewton1Deriv(NLP1* p, TOLS t): 
    OptConstrNewtonLike(p->getDim(),t), mem_nlp(p) {;}
  
  virtual ~OptConstrNewton1Deriv(){}

private:
  NLP1* mem_nlp;
  
protected:
  /**
   * @ returns a pointer to an NLP1
   */
  NLP1* nlprob() const { return mem_nlp; }
};

/**
 * Constrained Newton classes that require an NLP2
 */

class OptConstrNewton2Deriv: public OptConstrNewtonLike {
public:

  OptConstrNewton2Deriv() {}

  OptConstrNewton2Deriv(NLP2* p): 
    OptConstrNewtonLike(p->getDim()), mem_nlp(p){;}

  OptConstrNewton2Deriv(NLP2* p, UPDATEFCN u): 
    OptConstrNewtonLike(p->getDim(),u), mem_nlp(p){;}

  OptConstrNewton2Deriv(NLP2* p, TOLS t): 
    OptConstrNewtonLike(p->getDim(),t), mem_nlp(p){;}
  
  virtual ~OptConstrNewton2Deriv(){}

private:
  NLP2* mem_nlp;
  
protected:
  /**
   * @ returns a pointer to an NLP1
   */
  NLP1* nlprob() const { return mem_nlp;}
  /**
   * @ returns a pointer to an NLP2
   */
  NLP2* nlprob2() const { return mem_nlp;}
};

} // namespace OPTPP

#endif
