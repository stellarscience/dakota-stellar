#ifndef OptBCNewton_h
#define OptBCNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "OptBCNewtonLike.h"

namespace OPTPP {

/**
 * OptBCNewton is a derived class of OptBCNewtonLike.
 * OptBCNewton implements a bound constrained Newton method. 
 * These methods will use the active set method.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 03/2007
 */

class OptBCNewton: public OptBCNewton2Deriv {
protected:

 public:
 /**
  * Default Constructor
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(): 
    OptBCNewton2Deriv()
    { std::cerr << "OptBCNewton :: instantiation \n";
      strcpy(method,"Bound constrained Newton");
    }

 /**
  * @param p a pointer to an NLP1.
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(NLP2* p): 
    OptBCNewton2Deriv(p)
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /**
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, TOLS t)
  */
  OptBCNewton(NLP2* p, UPDATEFCN u): 
    OptBCNewton2Deriv(p, u)
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  * @see OptBCNewton(NLP2* p)
  * @see OptBCNewton(NLP2* p, UPDATEFCN u)
  */
  OptBCNewton(NLP2* p, TOLS t): 
    OptBCNewton2Deriv(p, t)
    { strcpy(method,"Bound constrained Newton"); work_set = false; }

 /*
  * Destructor
  */
  virtual ~OptBCNewton(){;}

//----------------------------------
// These are defined elsewhere
//----------------------------------
  virtual int             checkDeriv();
  virtual void            initHessian();
  virtual void            printStatus(char *);
  virtual real            stepTolNorm() const;
  Teuchos::SerialSymDenseMatrix<int,double> updateH(Teuchos::SerialSymDenseMatrix<int,double>& H, int k);

protected:
//  NLP2*   nlprob2() const {return nlp; }
};

} // namespace OPTPP

#endif
