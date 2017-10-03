#include "TestMean.h"

TestMean::TestMean()
{
}

TestMean::~TestMean()
{
}

void TestMean::run()
{
  testMean();
  testMeanDouble();
  testMeanArray();
}

void TestMean::testMean()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the default constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // create a mean and test the value
  Mean default_mean;
  _test( default_mean.value() == 0 );
}

void TestMean::testMeanDouble()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the double constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Mean double_mean( 0.5 );
  _test( double_mean.value() == 0.5 );
}

void TestMean::testMeanArray()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the array constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  std::vector<double> numbers;
  for( int i = 1; i <= 10; i++ ) {
    numbers.push_back( (double)i );
  }
  Mean number_mean( numbers );
  _test( number_mean.value() == 5.5 );
}

