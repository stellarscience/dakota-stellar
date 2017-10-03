#ifndef LSQNLF_h
#define LSQNLF_h
/*----------------------------------------------------------------------
 Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
 DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 work by or on behalf of the U.S. Government.

 J. C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#include "NLP2.h"

extern "C" {
  double get_cpu_time();
  double get_wall_clock_time();
}

namespace OPTPP {

typedef void (*INITFCN)(int, Teuchos::SerialDenseVector<int,double>&);

typedef CompoundConstraint* (*INITCONFCN)(int);

typedef void (*USERFCNLSQ0)(int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, int&);

typedef void (*USERFCNLSQ0V)(int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, int&, void* v);

typedef void (*USERFCNLSQ1)(int, int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&, int&);

typedef void (*USERFCNLSQ1V)(int, int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&, int&, void* v);

typedef void (*USERNLNCON0)(int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, int&);

typedef void (*USERNLNCON1)(int, int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&, int&);

typedef void (*USERNLNCON2)(int, int, const Teuchos::SerialDenseVector<int,double>&, 
  Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&, 
  OptppArray<Teuchos::SerialSymDenseMatrix<int,double> >&, int&);


/**
 * LSQNLF is a derived class of NLP2.  The LSQNLF class implements the 
 * function,  gradient, and Hessian evaluator for a least square
 * objective.
 *
 * @author P. J. Williams, Sandia National Laboratories
 * @date Last modified 11/2005
 */

class LSQNLF: public NLP2 {

protected:
  USERFCNLSQ0  fcn0;	        ///< User-defined objective function
  USERFCNLSQ0V fcn0_v;	        ///< User-defined objective function w/ void ptr
  USERFCNLSQ1  fcn1;	        ///< User-defined objective function
  USERFCNLSQ1V fcn1_v;	        ///< User-defined objective function w/ void ptr
  USERNLNCON1 confcn;		///< User-defined nonlinear constraints
  INITFCN init_fcn;		///< User-defined initfcn for the obj. function
  INITCONFCN init_confcn;	///< User-defined initfcn for the constraints
  bool init_flag;		///< Has the function been initialized?
  bool Jacobian_current;	///< Has the function been updated?
  int lsqterms_;   		///< Number of least square terms in objective
  Teuchos::SerialDenseVector<int,double> fvector; ///< Vector of objective function values 
  Teuchos::SerialDenseMatrix<int,double> Jacobian_;	///< Jacobian_ of objective function residuals
  Teuchos::SerialDenseMatrix<int,double> partial_jac;
  void* vptr; 			///< Void pointer 

  static void f0_helper(int n, const Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& f, 
         int& result, void *v) 
        {LSQNLF *o = (LSQNLF*)v; (*o->fcn0)(n,xc,f,result);}
  static void f1_helper(int m, int n, const Teuchos::SerialDenseVector<int,double>& xc, 
         Teuchos::SerialDenseVector<int,double>& f, Teuchos::SerialDenseMatrix<int,double>& g, int& result, void *v) 
        {LSQNLF *o = (LSQNLF*)v; (*o->fcn1)(m,n,xc,f,g,result);}

private:
  Teuchos::SerialDenseVector<int,double> tempF;   ///< Vector of objective function residuals 
  Teuchos::SerialDenseVector<int,double> specLSQF;///< Vector of objective function residuals 

public:
  // Constructor

#ifdef OPTPP_HAVE_MPI

  LSQNLF(): 
     NLP2(){;}
  LSQNLF(int ndim, int lsqterms): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(0), fcn1_v(0), confcn(0),
    init_fcn(0), init_confcn(0), init_flag(false), Jacobian_current(false), 
    lsqterms_(lsqterms), fvector(lsqterms), 
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
     {
	 fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	 SpecFlag = Spec1;
     }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0 f, INITFCN i, 
	  CompoundConstraint* constraint = 0): 
    NLP2(ndim, constraint), fcn0(f), fcn0_v(f0_helper), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(0), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1 f, INITFCN i, 
	  CompoundConstraint* constraint = 0): 
    NLP2(ndim, constraint), fcn0(0), fcn0_v(0), fcn1(f), fcn1_v(f1_helper), 
    confcn(0), init_fcn(i), init_confcn(0), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0 f, INITFCN i, INITCONFCN c): 
    NLP2(ndim), fcn0(f), fcn0_v(f0_helper), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1 f, INITFCN i, INITCONFCN c): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(f), fcn1_v(f1_helper), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }
  /// Alternate function pointers with user-supplied void function pointer
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0V f, INITFCN i, INITCONFCN c, void* v): 
    NLP2(ndim), fcn0(0), fcn0_v(f), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(v),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1V f, INITFCN i, INITCONFCN c, void* v): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(0), fcn1_v(f), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(v),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = Spec1;
    }

#else

  LSQNLF(): 
     NLP2(){;}
  LSQNLF(int ndim, int lsqterms): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(0), fcn1_v(0), confcn(0),
    init_fcn(0), init_confcn(0), init_flag(false), Jacobian_current(false), 
    lsqterms_(lsqterms), fvector(lsqterms), 
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
     {
	 fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	 SpecFlag = NoSpec;
     }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0 f, INITFCN i, 
	  CompoundConstraint* constraint = 0): 
    NLP2(ndim, constraint), fcn0(f), fcn0_v(f0_helper), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(0), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1 f, INITFCN i, 
	  CompoundConstraint* constraint = 0): 
    NLP2(ndim, constraint), fcn0(0), fcn0_v(0), fcn1(f), fcn1_v(f1_helper), 
    confcn(0), init_fcn(i), init_confcn(0), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0 f, INITFCN i, INITCONFCN c): 
    NLP2(ndim), fcn0(f), fcn0_v(f0_helper), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1 f, INITFCN i, INITCONFCN c): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(f), fcn1_v(f1_helper), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(this),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }
  /// Alternate function pointers with user-supplied void function pointer
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ0V f, INITFCN i, INITCONFCN c, void* v): 
    NLP2(ndim), fcn0(0), fcn0_v(f), fcn1(0), fcn1_v(0), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(v),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }
  LSQNLF(int ndim, int lsqterms, USERFCNLSQ1V f, INITFCN i, INITCONFCN c, void* v): 
    NLP2(ndim), fcn0(0), fcn0_v(0), fcn1(0), fcn1_v(f), 
    confcn(0), init_fcn(i), init_confcn(c), init_flag(false), 
    Jacobian_current(false), lsqterms_(lsqterms), fvector(lsqterms),
    Jacobian_(lsqterms,ndim), partial_jac(lsqterms,ndim), vptr(v),
    tempF(lsqterms), specLSQF(lsqterms)
    { 
	fvector = 1.0e30;  Jacobian_ = 1.0e30; tempF = 1.0e30;
	SpecFlag = NoSpec;
    }

#endif

  // Destructor
  virtual ~LSQNLF() {;}                     

  void setFcnResidual(Teuchos::SerialDenseVector<int,double>& f) {tempF = f;}
  Teuchos::SerialDenseVector<int,double> getFcnResidual() const {return tempF;}

  /// Reset parameters 
  virtual void reset();          
  /// Initialize selected function
  virtual void initFcn();          
  /// Evaluate the function, gradient, and Hessian
  virtual void eval();            
  /// Evaluate the function 
  virtual real evalF();                         	
  /// Evaluate the function at x 
  virtual real evalF(const Teuchos::SerialDenseVector<int,double>& x);    	
  /// Evaluate a finite-difference gradient 
  virtual Teuchos::SerialDenseVector<int,double> evalG();              		
  /// Evaluate a finite-difference gradient at x 
  virtual Teuchos::SerialDenseVector<int,double> evalG(const Teuchos::SerialDenseVector<int,double>& x);  	
  /// Evaluate the Hessian 
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalH();              	

  /// Evaluate the Lagrangian at x
  virtual real evalLagrangian(const Teuchos::SerialDenseVector<int,double>& x, 
    Teuchos::SerialDenseVector<int,double>& mult, const Teuchos::SerialDenseVector<int,double>& type) ;
  /// Evaluate the gradient of the Lagrangian at x
  virtual Teuchos::SerialDenseVector<int,double> evalLagrangianGradient(const Teuchos::SerialDenseVector<int,double>& x, const Teuchos::SerialDenseVector<int,double>& mult, const Teuchos::SerialDenseVector<int,double>& type) ;

  /// Evaluate nonlinear constraints at x
  virtual Teuchos::SerialDenseVector<int,double> evalCF(const Teuchos::SerialDenseVector<int,double>& x);  	
  /// Evaluate gradient of nonlinear constraints at x
  virtual Teuchos::SerialDenseMatrix<int,double> evalCG(const Teuchos::SerialDenseVector<int,double>& x);  	
private:
  /// Evaluate hessian at x 
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalH(Teuchos::SerialDenseVector<int,double>& x); 	
  /// Evaluate constraint hessian at x 
  virtual Teuchos::SerialSymDenseMatrix<int,double> evalCH(Teuchos::SerialDenseVector<int,double> &x);	
  virtual OptppArray<Teuchos::SerialSymDenseMatrix<int,double> > evalCH(Teuchos::SerialDenseVector<int,double> &x,
    int darg);
  virtual void evalC(const Teuchos::SerialDenseVector<int,double>& x);  	

  /// Construct forward finite-difference Jacobian of objective funtion 
  Teuchos::SerialDenseMatrix<int,double> LSQFDJac(const Teuchos::SerialDenseVector<int,double>& sx, 
    const Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& fx, 
    Teuchos::SerialDenseMatrix<int,double>& partial_jac);
  /// Construct backward finite-difference Jacobian of objective funtion 
  Teuchos::SerialDenseMatrix<int,double> LSQBDJac(const Teuchos::SerialDenseVector<int,double>& sx, 
    const Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& fx, 
    Teuchos::SerialDenseMatrix<int,double>& partial_jac);
  /// Construct central finite-difference Jacobian of objective funtion 
  Teuchos::SerialDenseMatrix<int,double> LSQCDJac(const Teuchos::SerialDenseVector<int,double>& sx, 
    const Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& fx, 
    Teuchos::SerialDenseMatrix<int,double>& partial_jac);
};

} // namespace OPTPP
#endif
