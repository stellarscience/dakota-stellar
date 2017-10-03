#ifndef Appl_Data_h
#define Appl_Data_h

/*----------------------------------------------------------------------
  Copyright (c) 2001, Sandia Corporation.   Under the terms of Contract 
  DE-AC04-94AL85000, there is a non-exclusive license for use of this 
  work by or on behalf of the U.S. Government.
 ----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_STD
#include <cstring>
#else
#include <string.h>
#endif

#include "globals.h"
#include "OptppArray.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;

/**
 * @author J. C. Meza, Sandia National Laboratories, meza@ca.sandia.gov
 * @note Modified by P.J. Williams 02/2006
 */

namespace OPTPP {

class Appl_Data {
private :
  /// Dimension of the problem
  int             	dimension;		
  /// Current point
  Teuchos::SerialDenseVector<int,double>    *xparm;		
  /// Objective function value 
  double          function_value;		
  /// Gradient of the objective function 
  SerialDenseVector<int, double>    *gradient;		
  /// Hessian of the objective function 
  SerialSymDenseMatrix<int, double> *Hessian;		
  /// Constraint value 
  SerialDenseVector<int, double>    *constraint_value;	
  /// Gradient of the constraints 
  SerialDenseMatrix<int,double>         *constraint_gradient;	
  /// Hessian of the constraints 
  OptppArray<Teuchos::SerialSymDenseMatrix<int, double> > *constraint_Hessian; 
  /// Residuals of the least square objective function 
  SerialDenseVector<int, double>    *lsq_residuals;	
  /// Jacobian of the least square objective function 
  SerialDenseMatrix<int, double>          *lsq_jacobian;        
  /// Is the function value current? 
  bool            function_current;		
  /// Is the gradient current? 
  bool            gradient_current;		
  /// Is the Hessian current? 
  bool            Hessian_current;		

public:
  /**
   * Default Constructor
   */
  Appl_Data();

  /**
   * Destructor
   */
  ~Appl_Data();

  void reset();

  bool Compare(const Teuchos::SerialDenseVector<int,double>&);
  bool getF(const Teuchos::SerialDenseVector<int,double>&, double&);
  bool getGrad(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&);
  bool getHess(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialSymDenseMatrix<int,double>&);

  bool getCF(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&);
  bool getCGrad(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&);
  bool getCHess(const Teuchos::SerialDenseVector<int,double>&, OptppArray<Teuchos::SerialSymDenseMatrix<int,double> >&);

  bool getLSQF(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseVector<int,double>&);
  bool getLSQJac(const Teuchos::SerialDenseVector<int,double>&, Teuchos::SerialDenseMatrix<int,double>&);
  

  /// Update the objective function value
  void update(int,int,const Teuchos::SerialDenseVector<int,double>&,double);
  /// Update the objective function and gradient
  void update(int,int,const Teuchos::SerialDenseVector<int,double>&,double,Teuchos::SerialDenseVector<int,double>&);
  /// Update the objective function, gradient, and Hessian
  void update(int,int,const Teuchos::SerialDenseVector<int, double>&,double,
	      Teuchos::SerialDenseVector<int,double>&,Teuchos::SerialSymDenseMatrix<int,double>&);

  /// Update the nonlinear constraint functions
  void constraint_update(int,int,int,const Teuchos::SerialDenseVector<int,double>&,
			 Teuchos::SerialDenseVector<int,double>&);
  /// Update the nonlinear constraint functions and Jacobian
  void constraint_update(int,int,int,const Teuchos::SerialDenseVector<int,double>&,
			 Teuchos::SerialDenseVector<int,double>&,Teuchos::SerialDenseMatrix<int,double>&);
  /// Update the nonlinear constraint functions, Jacobian, and Hessians
  void constraint_update(int,int,int,const Teuchos::SerialDenseVector<int,double>&,
			 Teuchos::SerialDenseVector<int,double>&,Teuchos::SerialDenseMatrix<int,double>&,OptppArray<Teuchos::SerialSymDenseMatrix<int,double> >&);

  /// Update the least square residuals 
  void lsq_update(int,int,int,const Teuchos::SerialDenseVector<int,double>&,Teuchos::SerialDenseVector<int,double>&);

  /// Update the least square residuals and Jacobian
  void lsq_update(int,int,int,const Teuchos::SerialDenseVector<int,double>&,
		  Teuchos::SerialDenseVector<int,double>&,Teuchos::SerialDenseMatrix<int,double>&);
};

} // namespace OPTPP

#endif
