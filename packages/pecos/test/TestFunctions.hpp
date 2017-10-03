#include "pecos_data_types.hpp"

using namespace Pecos;

void get_genz_coefficients( int num_dims, Real factor, int c_type, RealVector &c, RealVector &w )
{
  c.resize( num_dims );
  w.resize( num_dims );
  switch ( c_type )
    {
    case 0: // Linear decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.0;
	    c[d] = ( (Real)d + 4.5 ) / (Real)num_dims;
	    csum += c[d]; 
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    case 1: // Quadratic decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.0;
	    c[d] = 1.0 / (Real)( ( d + 1 ) * ( d + 1 ) );
	    csum += c[d]; 
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    case 2: // Exponential decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.;
	    c[d] = std::exp( (Real)(d+1)*std::log( 1.e-8 ) / (Real)num_dims );
	    csum += c[d];
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    default:
      throw(std::runtime_error("GetCoefficients() ensure type in [0,1,2]"));
    };
}

double genz(String an_comp, RealVector &xC)
{
  int coeff_type, fn_type;
  Real factor;
  if (an_comp == "os1")
    { coeff_type = 0; fn_type = 0; factor = 4.5; }
  else if (an_comp == "os2")
    { coeff_type = 1; fn_type = 0; factor = 4.5; }
  else if (an_comp == "os3")
    { coeff_type = 2; fn_type = 0; factor = 4.5; }
  else if (an_comp == "cp1")
    { coeff_type = 0; fn_type = 1; factor = .25; }
  else if (an_comp == "cp2")
    { coeff_type = 1; fn_type = 1; factor = .25; }
  else if (an_comp == "cp3")
    { coeff_type = 2; fn_type = 1; factor = .25; }
  else {
    std::cerr << "Error: option for genz() is not yet implemented : "<<an_comp
	 << std::endl;
    std::terminate();
  } 

  RealVector c, w;
  int numVars = xC.length();
  get_genz_coefficients( numVars, factor, coeff_type, c, w );
  Real fnVals, pi = 4.0 * std::atan( 1.0 );

  // **** f:
  switch (fn_type) {
    case 0: {
      fnVals = 2.0 * pi * w[0];
      for ( int d = 0; d < numVars; d++ ){
	fnVals += c[d] * xC[d];
      }
      fnVals = std::cos( fnVals );
      break;
    }
    case 1:{
      fnVals = 1.0;
      for ( int d = 0; d < numVars; d++ )
	{
	  fnVals += c[d]* xC[d];
	}
      fnVals = 1.0 / std::pow( fnVals, (Real)(numVars+1) );
      break;
    }
  }

  return fnVals; // no failure
}

double gerstner(String an_comp, RealVector &xC)
{

  if (xC.length() != 2 )
    throw(std::runtime_error("Pecos::gerstner() needs 2 variables\n"));

  short test_fn;
  Real x_coeff, y_coeff, xy_coeff;
  if (an_comp        == "iso1") {
    test_fn = 1; x_coeff = y_coeff = 10.;
  }
  else if (an_comp   == "iso2") {
    test_fn = 2; x_coeff = y_coeff = xy_coeff = 1.;
  }
  else if (an_comp   == "iso3") {
    test_fn = 3; x_coeff = y_coeff = 10.;
  }
  else if (an_comp == "aniso1") {
    test_fn = 1; x_coeff =  1.; y_coeff = 10.;
  }
  else if (an_comp == "aniso2") {
    test_fn = 2; x_coeff =  1.; y_coeff = xy_coeff = 10.;
  }
  else if (an_comp == "aniso3") {
    test_fn = 3; x_coeff = 10.; y_coeff = 5.;
  }
  else {
    throw(std::runtime_error("Pecos::gerstner() ensure type is iso1...3/aniso1...3\n"));
  }

  Real fnVals, x=xC[0], y=xC[1];
  switch (test_fn) {
    case 1:
      fnVals = x_coeff * std::exp(-x*x) + y_coeff * std::exp(-y*y);
      break;
    case 2:
      fnVals = x_coeff * std::exp(x) + y_coeff * std::exp(y)
              + xy_coeff * std::exp(x*y);
      break;
    case 3:
      fnVals = std::exp(-x_coeff * x * x - y_coeff * y * y);
      break;
  }

  return fnVals;

}

double custPol(String an_comp, RealVector &xC)
{
  int fn_type;
  RealVector cMain, cCross;
  int numVars = xC.length();

  cMain.resize( 2*numVars );
  cCross.resize( numVars*(numVars-1)/2 );

  if (an_comp == "pol2") { 
    /* 2-nd order polynomial */
    fn_type = 0;
    for ( int d = 0; d < numVars; d++ ) {
      cMain[d        ] = ( (Real) d + 1.0 ) / (Real) numVars ;
      //cMain[d        ] = 0.0;
      cMain[d+numVars] = ( (Real) d + 5.0 ) / (Real) numVars ;
    }
    int idx=0;
    for ( int d1 = 1; d1 < numVars; d1++ ) {
      for ( int d2 = 0; d2 < d1; d2++,idx++ ) {
        cCross[idx] = ( (Real) (d1+1.0)*(d2+1.0) + 1.0 ) / ((Real) numVars*numVars) ;
        //cCross[idx] = 0.0;
      }
    }
  }
  else {
    std::cerr << "Error: option for custPol() is not yet implemented : "<<an_comp
	 << std::endl;
    std::terminate();
  } 

  // **** f:
  Real fnVals;
  switch (fn_type) {
    case 0: {
      fnVals = 0.0;
      for ( int d = 0; d < numVars; d++ ){
	fnVals += cMain[d        ] * xC[d] * xC[d];
	fnVals += cMain[d+numVars] * xC[d] ;
        //printf("%d: %e x %e\n",d,cMain[d+numVars],xC[d]);
      }
      int idx=0;
      for ( int d1 = 1; d1 < numVars; d1++ ) {
        for ( int d2 = 0; d2 < d1; d2++,idx++ ) {
          fnVals += cCross[idx] * xC[d1] * xC[d2] ;
          //printf("%d,%d: %e x %e x %e\n",d1,d2,cCross[idx],xC[d1],xC[d2]);
        }
      }
      //printf("fnct=%e\n----------------\n",fnVals);
      break;
    }
  }

  return fnVals; // no failure

}
