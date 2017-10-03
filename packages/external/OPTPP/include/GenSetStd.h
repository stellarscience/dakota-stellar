#ifndef GenSetStd_h
#define GenSetStd_h
//------------------------------------------------------------------------
// Generating Set Classes - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBase.h"

namespace OPTPP {

//------------------------------------------------------------------------
// Standard Basis with 2*n vectors
//------------------------------------------------------------------------
class GenSetStd : public GenSetBase {  // Standard with 2*n vectors

 public:
  virtual std::string classnm() { return "GenSetStd";}

  // default constructor;
  GenSetStd(){};

  // Vsize specific constructor
  GenSetStd(int n) : GenSetBase(n) { 
    setSize(2*n);
    initActive();
  };

  void generate(int i, double a, Teuchos::SerialDenseVector<int,double> &x, Teuchos::SerialDenseVector<int,double> &y);
  ///< Stores the search direction in the vector y

  // overloaded pruning methods - virtual in base
  int init(){ return 0;}    ///< Computes initial generating set D
  int init(Teuchos::SerialDenseVector<int,double>& pV);  
  int update(){ return 0;}    ///< Updates D on each iteration
  int update(Teuchos::SerialDenseVector<int,double>& pV);   
  bool prunes() { return true; }

  ///< Destructor
  virtual ~GenSetStd() {;}  
};

} // namespace OPTPP

#endif
