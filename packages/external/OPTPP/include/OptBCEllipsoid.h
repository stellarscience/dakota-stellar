#ifndef OptBCEllipsoid_h
#define OptBCEllipsoid_h

#include "Opt.h"

namespace OPTPP{

/**
 * Bound Constrained Newton abstract data classes.
 * OptBCEllipsoid implements the bound constrained ellipsoid method.
 *
 * @author J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov 
 * @author Charles Tong
 * @note Modified by P.J. Williams to accept new constrained nlp objects 
 */

class OptBCEllipsoid: public OptimizeClass {
private:
  NLP1*	mem_nlp;			///< Pointer to an NLP1*
  int	nlp_type;			///< NLP identifier 
  
protected:
  /// Pointer to an NLP1*
  NLP1*   nlprob() const {return mem_nlp;} 
  /// Function lower/upper bounds 
  double  fval_lowbound, fval_upbound; 
  /// Radius of the initial sphere
  double  initial_radius;              
  /// Scaling indicator
  int     xscal_flag;                  
  /// Deep cut indicator
  int     deepcutflag;                 

public:

//----------------------------------------------------------------------
// Constructors & Destructors
//----------------------------------------------------------------------

  /**
    * Default constructor
    */
  OptBCEllipsoid() {}

  /**
    * @param p a pointer to an NLP1 object
    * @see OptBCEllipsoid(NLP1* p, UPDATEFCN u)
    * @see OptBCEllipsoid(NLP1* p, TOLS t)
    */
  OptBCEllipsoid(NLP1* p): OptimizeClass(p->getDim()), mem_nlp(p), 
      initial_radius(-1.0e0), xscal_flag(0), deepcutflag(0) {nlp_type = 1;}

  /**
    * @param p a pointer to an NLP1 object
    * @param u an update function 
    * @see OptBCEllipsoid(NLP1* p)
    * @see OptBCEllipsoid(NLP1* p, UPDATEFCN u)
    */
  OptBCEllipsoid(NLP1* p, UPDATEFCN u): OptimizeClass(p->getDim()), 
      mem_nlp(p), initial_radius(-1.0e0), xscal_flag(0), deepcutflag(0) 
      {update_fcn = u; nlp_type = 1;}

  /**
    * @param p a pointer to an NLP1 object
    * @param t a TOLS object
    * @see OptBCEllipsoid(NLP1* p)
    * @see OptBCEllipsoid(NLP1* p, UPDATEFCN u)
    */
  OptBCEllipsoid(NLP1* p, TOLS t): OptimizeClass(p->getDim(),t), 
      mem_nlp(p), initial_radius(-1.0e0), xscal_flag(0), deepcutflag(0) 
      {nlp_type = 1;}

  /**
    * Destructor
    */
  ~OptBCEllipsoid(){}

//----------------------------------------------------------------------------
// these functions have to be defined
//----------------------------------------------------------------------------

  virtual void         acceptStep(int k, int step_type);
  /// Check to see if algorithm satisfies the convergence criteria 
  virtual int          checkConvg();
  Teuchos::SerialDenseVector<int,double> computeSearch(Teuchos::SerialSymDenseMatrix<int,double>&) {exit(-1); return xprev;}
  virtual void         optimize(); 	///< Call the optimization method
  virtual void         readOptInput(); 	///< Read user-specified input options
  /// Print status of the bound constrained ellipsoidal method
  virtual void         printStatus(char *);

  virtual void         updateModel(int k, int ndim, Teuchos::SerialDenseVector<int,double> x) 
      {OptimizeClass::defaultUpdateModel(k, ndim, x);}
  virtual void         reset(); 

//----------------------------------------------------------------------------
// These are defined locally 
//----------------------------------------------------------------------------

  /// Sets up the optimization method 
  virtual void         initOpt();

  /// Sets the initial ellipsoid radius
  void         setInitialEllipsoid(double rad) {initial_radius = rad;}

  /// Computes feasibility of the constraints 
  double       computeFeasibility(Teuchos::SerialDenseVector<int,double>&);

  /// Picks the row corresponding to the most infeasible constraints
  Teuchos::SerialDenseVector<int,double> computeConstraintSubgradient(Teuchos::SerialDenseVector<int,double>&);

  /// Step taken if the current x is infeasible 
  int          infeasibilityStep(Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialSymDenseMatrix<int,double> &, 
                 double &);

  /// Deep cut step for upper bound
  int          halfSpaceStep(Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialSymDenseMatrix<int,double> &,double &);

  /// Set the scaling vector for x 
  void         setXScale(Teuchos::SerialDenseVector<int,double> &x) {OptimizeClass::setXScale(x); 
                 xscal_flag = 1;}
  /// Set deepcutflag = 1 
  void         setDeepCut()               {deepcutflag = 1;}

  /// Set deepcutflag = 0 
  void         resetDeepCut()             {deepcutflag = 0;}

  /// Given x, compute the gamma function 
  double       computeGamma(double);

};

} // namespace OPTPP

#endif
