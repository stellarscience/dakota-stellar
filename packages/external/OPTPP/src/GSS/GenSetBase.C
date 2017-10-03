//------------------------------------------------------------------------
// Generating Set Class - for use with OptGSS
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBase.h"

using Teuchos::SerialDenseVector;
using Teuchos::SerialDenseMatrix;

using std::cerr;
using std::endl;

namespace OPTPP {
//
// Methods for Generating Directions
//


bool GenSetBase::generateAll(SerialDenseMatrix<int,double>& M, SerialDenseVector<int,double>& X, double D){ 
  if (Size<=0 || Vdim<=0) {
    cerr << "***ERROR: GenSetBase::generateAll(SerialDenseMatrix<int,double>,...) "
	 << "called with size=" << Size << ", vdim=" << Vdim << endl;
    return false;
  }
  if (M.numCols() != Size || M.numRows() != Vdim) {
    cerr << "***ERROR: GenSetBase::generateAll(SerialDenseMatrix<int,double>,...) "
	 << "dimesion of M expected to be "
	 << Vdim << "-by-" << Size 
	 << " but is " << M.numRows() << "-by-" << M.numCols()
	 << endl;
    return false;
  }
  SerialDenseVector<int,double> xi(Vdim);
  for (int i=0; i<Size; i++) {
    std::cout<<"Calling generate from GenSetBase.C"<<std::endl;
    generate(i+1, D, X, xi);
    // M.Column(i) = xi;
    for(int j=0; j< xi.length(); j++)
      {M(j,i) = xi(j);}
  
  }
  return true;
}

bool GenSetBase::generateAllActive(SerialDenseMatrix<int,double>& M, SerialDenseVector<int,double>& X, double D){ 
  if (Size<=0 || Vdim<=0 || nActive()<=0) {
    cerr << "***ERROR: GenSetBase::generateAllActive(SerialDenseMatrix<int,double>,...) "
	 << "called with size=" << Size << ", vdim=" << Vdim 
	 << " nActive = " << nActive() 
	 << endl;
    return false;
  }
  if (M.numCols() != nActive() || M.numRows() != Vdim ) {
    cerr << "***ERROR: GenSetBase::generateAllActive(SerialDenseMatrix<int,double>,...) "
	 << "dimesion of M expected to be "
	 << Vdim << "-by-" << nActive()
	 << " but is " << M.numRows() << "-by-" << M.numCols()
	 << endl;
    return false;
  }
  SerialDenseVector<int,double> xi(Vdim);
  for (int i=0; i<nActive(); i++) {
    generateActive(i, D, X, xi);
    //M.Column(i) = xi;
for(int j=0; j< xi.length(); j++)
      {M(j,i) = xi(j);}
  }
  return true;
}
/*
SerialDenseVector<int,double> GenSetBase::generate(int i) { 
  ///< returns d_i, the ith basis element
  SerialDenseVector<int,double> y(Vdim);
  y = 0;
  generate(i, 0, y, y); 
  return y;
}

void GenSetBase::generate(int i, SerialDenseVector<int,double>& y) {
  ///< Stores d_i in y
  y = 0;
  generate(i, 0, y, y);
}

SerialDenseVector<int,double> GenSetBase::generate(int i, double a, SerialDenseVector<int,double>& x) {
  ///< Returns  x + a * d_i, with d_i = ith basis element
  SerialDenseVector<int,double> y(Vdim);
  generate(i, a, x, y);
  return y;
}
*/
// -- purely virtual--  to be defined in  each derived class:
//
// generate(int i, double a, SerialDenseVector<int,double> &x, SerialDenseVector<int,double> &y)
// 
// ///< Set  y = x + a * d_i; should allow y and x to be same vector.
//


SerialDenseMatrix<int,double> GenSetBase::pllMesh(int P, SerialDenseVector<int,double>& xc, SerialDenseVector<int,double>& xn, double r) 
  // P : num points we want ( ~ num of processors)
  // xc: the current point; 
  // xn: the "newton" point, tells the dir in which the mesh is grown
  //  r: ||xc-xn||, if known.
  //
  // return matrix M with genset generated at y_k = xc+k(xn-xc)
  // with radius r_k = k^n*r_0, r_0=||xc-xn||/2
  //
  // *** current implementation not optimized for efficiency ***
  //
{

    int k = 0;
    SerialDenseVector<int,double> xk(xn.length()); 
    double       rk; 
    SerialDenseMatrix<int,double> M, A, B;

    // SerialDenseVector<int,double> ns = xn - xc;  // newton step
    SerialDenseVector<int,double> ns(xn.length());
    ns = xn;
    ns -= xc;

    int m = Vdim;
    int n = Size;

    // return at least the base point
    M.reshape(xn.length(),1);
    M = xn;
    int nump = P-1; 

    //--
    // main loop:
    //--
    if (r <= 0) r = std::sqrt(ns.dot(ns));
    double r0 = 0.25*r;
    double pert = r0*.05;
    while (nump > 0) {

      // generate points xc + k*newton_step + r^k*d(i), i=1..size()
      ++k;
      xk = xc;
      SerialDenseVector<int,double> AnotherTemp(ns.length());
      AnotherTemp = ns.scale(k);
      xk +=  AnotherTemp;
      rk = k*r0;
      A.reshape(n,m);
      A = generateAll(xk,rk);

      // perturbation to avoid potential overlaps
      int RAND_HALF = RAND_MAX / 2;
      for (int i=0; i<n; i++)
	for (int j=0; j<m; j++) {
	  int sig = (std::rand() > RAND_HALF)? 1 : -1;
	  double amp = std::rand() / RAND_MAX * pert;	
	  A(j,i) = A(j,i) + sig * amp ;
	}

      // after the first iteration want to add xk to M
      if (k>1) {
	B.reshape(M.numRows(),M.numCols());
	B=M;
	B.reshape(M.numRows(), M.numCols()+1);
	//B = M | xk;
	for(int i=0; i<xk.length(); i++)
	  B(i,M.numCols()) = xk(i);
	M.reshape(B.numRows(),B.numCols());
	M = B;
	--nump;
      }

      // addd to M as many columns as we can afford
      if (nump > n) {
	B=M;
	B.reshape(M.numRows(), M.numCols()+A.numCols());
	//B = M | A;
	for(int i = 0; i<M.numRows(); i++)
	  for(int j=M.numCols(); j< M.numCols() + A.numCols(); j++)
	    {B(i,j) = A(i,j);}

      }
      else if (nump>0) { 
	//C = A.SubMatrix(1,m,1,nump);
	//B = M | C;
	B = M;
	B.reshape(M.numRows(), M.numCols()+nump);
	for(int i=0; i<M.numRows();i++)
	  for(int j=M.numRows(); j<B.numCols(); j++)
	    {B(i,j) = A(i,j-M.numCols());}
      }
      
      M = B;
      nump = nump - n;

    }
    return M;
}

} // namespace OPTPP
