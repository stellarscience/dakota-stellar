//------------------------------------------------------------------------
// Generating Set Basis Class D = [I -I (1;1) (1;-1) (-1;1) (-1;-1)]
// for 2-d only!!!!
//-----------------------------------------------------------------------
// Basis vectors represented by integers
//------------------------------------------------------------------------

/*------------------------------------------------------------------------
 Copyright (c) 2003,
 Ricardo Oliva (raoliva@lbl.gov)
 Lawrence Berkeley National Laboratory
 ------------------------------------------------------------------------*/

#include "GenSetBox2d.h" 

using Teuchos::SerialDenseVector;
using std::cerr;

namespace OPTPP {

///> Stores the search direction in the vector y
void GenSetBox2d::generate(int i, double a, SerialDenseVector<int,double> &x, SerialDenseVector<int,double> &y)
{
  //  sets y = x + a * d[i] 
  if (i<1 || i>Size) {
    cerr << "Gen_Set_Box2d: Basis index out of range: " << i << "\n";
    return;
  }

  y == x;

  if (i<=Vdim)
    y(i) += a;  
  else if (i<=2*Vdim)
    y(i-Vdim) -= a;
  else {
    double w = a / std::sqrt(2.0);

    switch (i-2*Vdim) {
    case 1:
      y(1) += w;   y(2) += w;  break;
    case 2:
      y(1) += w;   y(2) -= w;  break;
    case 3:
      y(1) -= w;   y(2) += w;  break;
    case 4:
      y(1) -= w;   y(2) -= w;  break;
    }
  }
}

//--
// the pruning methods
//--

int GenSetBox2d::init(SerialDenseVector<int,double>& gX)  {

  ActiveIDs.resize(Size);
  for (int i=0; i<Size; i++) ActiveIDs(i) = i; 

  return update(gX);
}

int GenSetBox2d::update(SerialDenseVector<int,double>& gX)  {
  if (Size<1) {
    cerr << "GenSetBox2d Error: update() called on an empty GenSet\n";
    return -1;
  }

  //--
  // Update == Pruning
  //--
  // determine which search directions are descending
  // and sets Active accordingly
  //--
  int nIna = 0;
  nAct = 0; 
  ActiveIDs = 0; 
  InactiveIDs = 0;
  double gradangle = 0.0;

  // {I} ==> gX*d = gX(i)
  for (int i=0; i<Vdim; i++) {
    if (gX(i) <= gradangle) {
      ActiveIDs(++nAct) = i;
    } 
    else {
      InactiveIDs(++nIna) = i;
    }
  }

  // { -I }
  for (int i=Vdim; i<2*Vdim; i++) {
    if (gX(i-Vdim) >= gradangle) {
      ActiveIDs (++nAct) = i;
    }
    else {
      InactiveIDs(++nIna) = i;
    }
  }

  // { Corner Directions }
  for (int i=2*Vdim; i<Size; i++) {
    int s = i-2*Vdim;
    double dot;
    switch (s) { 
    case 1:  dot =  gX(0) + gX(1);  break;
    case 2:  dot =  gX(0) - gX(1);  break;
    case 3:  dot = -gX(0) + gX(1);  break;
    case 4:  dot = -gX(0) - gX(1);  break;
    default: dot = 0.0;
    }
    if (gradangle != 0.0) dot /= std::sqrt(2.0);

    if (dot < gradangle) {
      ActiveIDs (++nAct) = i;
    }
    else 
      InactiveIDs(++nIna) = i;

  } // for

  return 0;
}

} // namespace OPTPP
