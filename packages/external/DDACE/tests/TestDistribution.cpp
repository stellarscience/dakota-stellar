#include "TestDistribution.h"
#include "UniformDistribution.h"
#include "arrcmp.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h>

using namespace std;

TestDistribution::TestDistribution()
{
   // set the seed
   seed = 779;
   DistributionBase::setSeed( seed );

   // set lower and upper bounds
   lb = 5;
   ub = 10;

   //get the first deviate
   Distribution distArg( UniformDistribution( lb, ub ) );
   this->deviate = distArg.getDeviate();
}

TestDistribution::~TestDistribution()
{
}

void TestDistribution::run()
{
    testDistributionDefault(); 
    testDistributionArg(); 
    testGetDeviateNoProb();
    testGetDeviateWithProb();
    testGetCDF();
    testLowerBound();
    testUpperBound();
    testMean();
    testStdDev();
    testPrint();
    testPrintAttributes();
    testTypeName();
    testUsePseudoRandom();
    
}


void TestDistribution::testUsePseudoRandom() {
	
	DistributionBase::usePseudoRandom(true);
	
	DistributionBase::setSeed(0);
	double number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.000)<0.0001);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.001)<0.0001);
	
	DistributionBase::setSeed(500);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.500)<0.0001);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.501)<0.0001);
	
	DistributionBase::setSeed(998);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.998)<0.0001);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.999)<0.0001);
	number = DistributionBase::uniformUnitDeviate();
	_test(fabs(number - 0.0)<0.0001);	
	
	DistributionBase::usePseudoRandom(false);
	DistributionBase::setSeed(0);
	number = DistributionBase::uniformUnitDeviate();
	_test(number != 0);
	
	DistributionBase::setSeed(0);
		
}


void TestDistribution::testDistributionDefault()
{
   // use constructor Distribution();
   Distribution distEmpty;

   try {
     _test( distEmpty.lowerBound() == 0 );
  
     // if here then no exception was thrown, record failure.
     _test( false );
   } catch( std::exception& e ) {
     // if here then exception was thrown, record success
     _test( true );
   }
}

void TestDistribution::testDistributionArg()
{
   // use constructor Distribution( DistributionBase& base );
   Distribution distArg( UniformDistribution( lb, ub ) );

   _test( distArg.lowerBound() == lb );
   _test( distArg.upperBound() == ub );
}

void TestDistribution::testGetDeviateNoProb()
{
   DistributionBase::setSeed( seed );
   Distribution distArg( UniformDistribution( lb, ub ) );

   // test the returned deviate value
   _test( fabs(distArg.getDeviate() - this->deviate)  < 0.001 );
}

void TestDistribution::testGetDeviateWithProb()
{
   DistributionBase::setSeed( seed );
   Distribution distArg( UniformDistribution( lb, ub ) );

   // test the returned deviate value
   _test( distArg.getDeviate( 0.5 ) == (lb + (ub-lb)*0.5) );
}

void TestDistribution::testGetCDF()
{
   Distribution distArg( UniformDistribution( lb, ub ) );

   // test the returned CDF value

   double a = distArg.getCDF( 0.5 );
   double b = ((0.5-lb)/(ub-lb));
   _test( a == b );
   // _test( distArg.getCDF( 0.5 ) == ((0.5-lb)/(ub-lb)) );
}

void TestDistribution::testLowerBound()
{
   Distribution distArg( UniformDistribution( lb, ub ) );
   
   _test( distArg.lowerBound() == lb );
}

void TestDistribution::testUpperBound()
{
   Distribution distArg( UniformDistribution( lb, ub ) );

   _test( distArg.upperBound() == ub );
}

void TestDistribution::testMean()
{
   Distribution distArg( UniformDistribution( lb, ub ) );

   _test( distArg.mean() == lb+(ub-lb)/2 );
}

void TestDistribution::testStdDev()
{
   Distribution distArg( UniformDistribution( lb, ub ) );
   
   _test( fabs(distArg.stdDev() - sqrt(pow(ub-lb,2)/12.0)) < 1.0e-10 );
}

void TestDistribution::testPrint()
{
   char buf[128];

   Distribution distArg( UniformDistribution( lb, ub ) );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   memset( buf, 0, 128 );
#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "UNIFORM %d %d", (int)lb, (int)ub );
#else
   int numberOfCharsRead = sprintf( buf, "UNIFORM %d %d", (int)lb, (int)ub );
   assert(numberOfCharsRead <= 128);
#endif

   test_str = buf; 

   fout.open( "TestDistribution_Log" );
   if( fout )
   {
      distArg.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      _test( data == test_str );
      fin.close();
   }
   else
   {
      _test( false );
   }
}
  

void TestDistribution::testPrintAttributes()
{
   char buf[128];

   Distribution distArg( UniformDistribution( lb, ub ) );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   memset( buf, 0, 128 );
#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "distribution=\"uniform\" lower=\"%d\" upper=\"%d\"", (int)lb, (int)ub );
#else
   int numberOfCharsRead = 
     sprintf( buf, "distribution=\"uniform\" lower=\"%d\" upper=\"%d\"", (int)lb, (int)ub );
   assert(numberOfCharsRead <= 128);
#endif
   test_str = buf;

   fout.open( "TestDistribution_Log" );
   if( fout )
   {
      distArg.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      _test( data == test_str );
      fin.close();
   }
   else
   {
      _test( false );
   }
}

void TestDistribution::testTypeName()
{
   Distribution distArg( UniformDistribution( lb, ub ) );

   _test( distArg.typeName() == std::string( "UniformDistribution" ) );
}
