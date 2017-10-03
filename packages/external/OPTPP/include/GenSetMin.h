#ifndef GenSetMin_h
#define GenSetMin_h
//------------------------------------------------------------------------
// Generating Set Classes - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBase.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"

namespace OPTPP {

//------------------------------------------------------------------------
// Standard Basis with n+1 vectors
//------------------------------------------------------------------------
class GenSetMin : public GenSetBase { 

 public:
  virtual std::string classnm() { return "GenSetMin";}

  // default constructor;
  GenSetMin(){};

  // Size specific constructor
  GenSetMin(int n) : GenSetBase(n) { 
    setSize(n+1); 
    initActive();
  };

  // iteration methods - virtual in base

  void generate(int i, double a, Teuchos::SerialDenseVector<int,double> &x, Teuchos::SerialDenseVector<int,double> &y);
  ///< Stores the search direction in the vector y

  //
  // pruning methods - virtual in base
  //
  int init(){ return 0;}    ///< Computes initial generating set D
  int init(Teuchos::SerialDenseVector<int,double>& pV);     ///< Computes initial generating set D
  int update(){ return 0;}    ///< Updates D on each iteration
  int update(Teuchos::SerialDenseVector<int,double>& pV);   ///< Updates D on each iteration
  bool prunes() { return true;}

 
 ///< Destructor
  virtual ~GenSetMin() {;}  
};

} // namespace OPTPP

#endif
