#ifndef GenSetBase_h
#define GenSetBase_h
//------------------------------------------------------------------------
// Generating Set Class - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#define WANT_MATH

#include <string>
#include <iostream>
#ifdef HAVE_STD
#include <cfloat>
#include <cstdlib>
#else
#include <float.h>
#include <stdlib.h>
#endif


#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"


namespace OPTPP {

class GenSetBase {  // Generating Set Base Class

 protected:
  int  Vdim;

  int  Size;
  int  nAct;

  Teuchos::SerialDenseVector<int,double> ActiveIDs;
  Teuchos::SerialDenseVector<int,double> InactiveIDs;

 public:
  virtual std::string classnm() { return "GenSetBase";}; 

  // default constructor;
  GenSetBase() : Vdim(0),  Size(0), nAct(0) {};

  // Constructor with specific vector-size
  GenSetBase(int n) : Vdim(n), Size(0), nAct(0) {};

  /// Destructor
  virtual ~GenSetBase() {;}  

  // Basic init method -- call after default constr.
  void init(int vd) { Vdim = vd; }

  // Basic set/get methods
  void setSize(int s) { Size = s; }
  void setVdim(int n) { Vdim = n; }
  int size() { return Size; }
  int vdim() { return Vdim; }

  //--
  // Generating Methods
  //--

  // -- wrt ALL Directions --

  //  virtual Teuchos::SerialDenseVector<int,double> generate(int i);
  ///< Returns  d_i = ith element of D

  //  virtual void generate(int i, Teuchos::SerialDenseVector<int,double>& y);
  ///< Stores d_i in y

  //  virtual Teuchos::SerialDenseVector<int,double> generate(int i, double a, Teuchos::SerialDenseVector<int,double>& x);
  ///< Returns the vector y_i =  x + a*d_i
  
  virtual void generate(int i, double a, Teuchos::SerialDenseVector<int,double> &x, 
    Teuchos::SerialDenseVector<int,double> &y) = 0;  
  ///< Stores in y the vector  x + a*d_i

  // -- wrt ACTIVE Directions --
  /*
  virtual 
    Teuchos::SerialDenseVector<int,double> generateActive(int i)
    { return generate(activeID(i)); }

  virtual 
    void generateActive(int i, Teuchos::SerialDenseVector<int,double>& y) 
    { generate(activeID(i), y); }

  virtual 
    Teuchos::SerialDenseVector<int,double> generateActive(int i, double s, Teuchos::SerialDenseVector<int,double>& x)
    { return generate(activeID(i), s, x); }
  */
  virtual void generateActive(int i, double s, Teuchos::SerialDenseVector<int,double> &x, 
    Teuchos::SerialDenseVector<int,double> &y)
  {generate(activeID(i), s, x, y); }

  // -- wrt INACTIVE Directions --
  /*
  virtual 
    Teuchos::SerialDenseVector<int,double> generateInactive(int i)
    ///< Returns b_i, the ith element of the INACTIVE subset of D
    { 
      Teuchos::SerialDenseVector<int,double> v; 
      v = generate(inactiveID(i)); 
      return v; 
    }

  virtual 
    void generateInactive(int i, Teuchos::SerialDenseVector<int,double>& y)
    ///< Stores b_i in y 
    { generate(inactiveID(i), y); }

  virtual 
    Teuchos::SerialDenseVector<int,double> generateInactive(int i, double s, Teuchos::SerialDenseVector<int,double>& x)
    ///< Returns the vector x + s*b_i,
    { return generate(inactiveID(i), s, x); }
  */
  virtual 
    void generateInactive(int i, double s, Teuchos::SerialDenseVector<int,double> &x, Teuchos::SerialDenseVector<int,double> &y)
    ///< Stores in y the vector x + s*b_i,
  {generate(inactiveID(i), s, x, y); }

  //--
  // Pruning methods
  //--
  virtual void initActive() {  
    // call this in constructor of derived class 
    // after size of derived class has been set
    if (Size==0) { 
      std::cerr << "!!! ERROR: GenSetBase::initActive() called when size==0\n";
      return;
    }
    nAct = Size;
    ActiveIDs.resize(Size);
    for (int i=0; i<Size; i++) ActiveIDs(i) = i; 
    InactiveIDs.resize(Size); 
    InactiveIDs = 0;
  }

  virtual int nActive() { return nAct; }
  virtual int nInactive() { return (Size - nAct); }
  virtual int activeID(int j) { return static_cast<int>(ActiveIDs(j)); }
  virtual int inactiveID(int j) { return static_cast<int>(InactiveIDs(j)); }

  virtual int init(){ return 0;}    ///< Computes initial generating set D
  virtual int init(Teuchos::SerialDenseVector<int,double>& pV){ return 0;}    

  virtual int update(){ return 0;}    ///< Updates D on each iteration
  virtual int update(Teuchos::SerialDenseVector<int,double>& pV){return 0;}        

  virtual bool prunes(){return false;} 
  ///< switch to true if implementing pruning in derived class

  bool   generateAll(Teuchos::SerialDenseMatrix<int,double>& M, Teuchos::SerialDenseVector<int,double>& X, double Delta=1.0);
  Teuchos::SerialDenseMatrix<int,double> generateAll(Teuchos::SerialDenseVector<int,double>& X, double D=1.0) {
    Teuchos::SerialDenseMatrix<int,double> M(Vdim,Size);
    generateAll(M,X,D);
    return M;
  }
  Teuchos::SerialDenseMatrix<int,double> generateAll(double Delta=1.0) {
    Teuchos::SerialDenseVector<int,double> X(Vdim);
    X = 0;
    return generateAll(X,Delta); 
  }

  bool   generateAllActive(Teuchos::SerialDenseMatrix<int,double>& M, Teuchos::SerialDenseVector<int,double>& X, double Delta=1.0);
  Teuchos::SerialDenseMatrix<int,double> generateAllActive(Teuchos::SerialDenseVector<int,double>& X, double D=1.0) {
    int n = nActive();
    int m = Vdim;
    Teuchos::SerialDenseMatrix<int,double> M(m,n);
    generateAllActive(M,X,D);
    return M;
  }
  Teuchos::SerialDenseMatrix<int,double> generateAllActive(double Delta=1.0) {
    Teuchos::SerialDenseVector<int,double> X(Vdim);
    X = 0;
    return generateAllActive(X,Delta); 
  }

  Teuchos::SerialDenseMatrix<int,double> pllMesh(int P, Teuchos::SerialDenseVector<int,double>& xc, Teuchos::SerialDenseVector<int,double>& xn, double d=0.0);

}; // end of GenSetBase class

} // namespace OPTPP

#endif
