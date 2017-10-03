
#include <ctype.h>
#include <string>

#include <Teuchos_UnitTestHarness.hpp> 

#include "pecos_data_types.hpp"
#include "BasisPolynomial.hpp"

using namespace Pecos;

//----------------------------------------------------------------

TEUCHOS_UNIT_TEST(pochhammer, simple1)
{

  Real testval = -5.0;
  std::vector<Real> values;
  for( unsigned short n = 0; n < 10; ++n ) {
    values.push_back(BasisPolynomial::pochhammer(testval, n));
    //std::cout << n << "\t" << values[n] << std::endl;
  }

  TEST_EQUALITY(    1,  values[0] );
  TEST_EQUALITY(   -5,  values[1] );
  TEST_EQUALITY(   20,  values[2] );
  TEST_EQUALITY(  -60,  values[3] );
  TEST_EQUALITY(  120,  values[4] );
  TEST_EQUALITY( -120,  values[5] );
  TEST_EQUALITY(    0,  values[6] );
  TEST_EQUALITY(    0,  values[7] );
  TEST_EQUALITY(    0,  values[8] );
  TEST_EQUALITY(    0,  values[9] );
}

//----------------------------------------------------------------
