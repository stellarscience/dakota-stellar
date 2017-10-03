#ifndef OptBCFDNewton_h
#define OptBCFDNewton_h

/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.
 J.C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#include "OptBCNewtonLike.h"

namespace OPTPP {

/**
 * OptBCFDNewton is a derived class of OptBCNewtonLike.
 * OptBCFDNewton implements a Bound constrained Finite-Difference 
 * Newton method.  These methods will use the active set method.
 *
 * @author J.C. Meza, Lawrence Berkeley National Laboratory
 * @note Modified by P.J. Williams, pwillia@sandia.gov
 * @date Last modified 03/2007
 */

class OptBCFDNewton: public OptBCNewton1Deriv {
 protected:

 public:
 /**
  * Default constructor 
  */
  OptBCFDNewton(): 
    OptBCNewton1Deriv() 
    { std::cerr << "OptBCFDNewton :: instantiation \n";
      strcpy(method,"Bound constrained QNewton"); work_set = false; 
    }
 /**
  * Constructors
  * @param p a pointer to an NLP1.
  * @see OptBCFDNewton(NLP1* p, UPDATEFCN u)
  * @see OptBCFDNewton(NLP1* p, TOLS t)
  */
  OptBCFDNewton(NLP1* p): 
    OptBCNewton1Deriv(p)
    { strcpy(method,"Bound constrained FDNewton"); work_set = false; }
 /**
  * Constructors
  * @param p a pointer to an NLP1.
  * @param u a function pointer.
  * @see OptBCFDNewton(NLP1* p)
  * @see OptBCFDNewton(NLP1* p, TOLS t)
  */
  OptBCFDNewton(NLP1* p, UPDATEFCN u): 
    OptBCNewton1Deriv(p, u)
    { strcpy(method,"Bound constrained FDNewton"); work_set = false; }
 /**
  * @param p a pointer to an NLP1.
  * @param t tolerance class reference.
  * @see OptBCFDNewton(NLP1* p)
  * @see OptBCFDNewton(NLP1* p, UPDATEFCN u)
  */
  OptBCFDNewton(NLP1* p, TOLS t): 
    OptBCNewton1Deriv(p, t)
    { strcpy(method,"Bound constrained FDNewton"); work_set = false; }

 /**
  * Destructor
  */
  virtual ~OptBCFDNewton(){;}

//---------------------------------
// These are defined elsewhere
//---------------------------------

  virtual int   checkDeriv();
  Teuchos::SerialSymDenseMatrix<int,double> updateH(Teuchos::SerialSymDenseMatrix<int,double>& H, int k);

};

} // namespace OPTPP
#endif
