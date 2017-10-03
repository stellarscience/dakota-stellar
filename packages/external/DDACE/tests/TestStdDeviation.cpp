#include "TestStdDeviation.h"
#include <cmath>
using namespace std;

TestStdDeviation::TestStdDeviation()
{
}

TestStdDeviation::~TestStdDeviation()
{
}

void TestStdDeviation::run()
{
  testStdDeviationDouble();
  testStdDeviationArray();
}

void TestStdDeviation::testStdDeviationDouble()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the double constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // create a std deviation
  StdDeviation double_stdd( 0.5 );
  _test( double_stdd.value() == 0.5 );
}

void TestStdDeviation::testStdDeviationArray()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the array constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // create an array of values and use to
  // create a std deviation then check
  // the calculated std deviation
  std::vector<double> numbers;
  for( int i = 1; i <= 10; i++ ) {
 
   numbers.push_back( (double)i );
  }
  StdDeviation array_stdd( numbers );
  _test( closeEnough( array_stdd.value(), 3.02765, 0.0001 ) );
}

  
