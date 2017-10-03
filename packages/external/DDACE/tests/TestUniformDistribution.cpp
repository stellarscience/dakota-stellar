#include "TestUniformDistribution.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

TestUniformDistribution::TestUniformDistribution()
{
   // set the seed
   seed = 779;
   DistributionBase::setSeed( seed );

   // set lower and upper bounds
   lb = 5;
   ub = 10;

   /* get some deviates */
   UniformDistribution distnb;
   DistributionBase::setSeed( seed );
   this->deviates[0] = distnb.getDeviate();
   DistributionBase::setSeed( seed );

   UniformDistribution distwb( lb, ub );
   DistributionBase::setSeed( seed );
   this->deviates[1] = distwb.getDeviate();
   DistributionBase::setSeed( seed );



}

TestUniformDistribution::~TestUniformDistribution()
{
}

void TestUniformDistribution::run()
{
    testUniformDistributionNoBounds(); 
    testUniformDistributionWithBounds(); 
    testClone();
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
}
void TestUniformDistribution::testUniformDistributionNoBounds()
{
   // use constructor UniformDistribution();
   UniformDistribution distnb;

   _test( distnb.lowerBound() == 0 );
   _test( distnb.upperBound() == 1 );
}

void TestUniformDistribution::testUniformDistributionWithBounds()
{
   // use constructor UniformDistribution(double lowerBound, double upperBound);
   UniformDistribution distwb( lb, ub );

   _test( distwb.lowerBound() == lb );
   _test( distwb.upperBound() == ub );
}
 
void TestUniformDistribution::testClone()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );
   
   // call clone and test results
   UniformDistribution* rtn = (UniformDistribution*)distnb.clone();

   _test( rtn->lowerBound() == distnb.lowerBound() );
   _test( rtn->upperBound() == distnb.upperBound() );
   _test( rtn->mean() == distnb.mean() );
   _test( rtn->stdDev() == distnb.stdDev() );
   _test( rtn->typeName() == distnb.typeName() );

   // clean up memory
   delete rtn;

   // call clone and test results
   rtn = (UniformDistribution*)distwb.clone();

   _test( rtn->lowerBound() == distwb.lowerBound() );
   _test( rtn->upperBound() == distwb.upperBound() );
   _test( rtn->mean() == distwb.mean() );
   _test( rtn->stdDev() == distwb.stdDev() );
   _test( rtn->typeName() == distwb.typeName() );

   // clean up memory
   delete rtn;
}

void TestUniformDistribution::testGetDeviateNoProb()
{
   DistributionBase::setSeed( seed );
   UniformDistribution distnb;
   DistributionBase::setSeed( seed );
   double deviate = distnb.getDeviate();
   _test( fabs(deviate - this->deviates[0]) < 0.001 );
   DistributionBase::setSeed( seed );


   UniformDistribution distwb( lb, ub );
   DistributionBase::setSeed( seed );
   deviate = distwb.getDeviate();
   _test( fabs(deviate - this->deviates[1])  < 0.001 );
   DistributionBase::setSeed( seed );
}

void TestUniformDistribution::testGetDeviateWithProb()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   // test the returned deviate value
   _test( distnb.getDeviate( 0.5 ) == 0.5 );
   _test( distwb.getDeviate( 0.5 ) == (lb + (ub-lb)*0.5) );
}

void TestUniformDistribution::testGetCDF()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   // test the returned CDF value
   _test( distnb.getCDF( 0.5 ) == 0.5 );

   double a = distwb.getCDF( 0.5 );
   double b = ((0.5-lb)/(ub-lb));
   _test( a == b );
   // _test( distwb.getCDF( 0.5 ) == ((0.5-lb)/(ub-lb)) );
}

void TestUniformDistribution::testLowerBound()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   _test( distnb.lowerBound() == 0 );
   
   _test( distwb.lowerBound() == lb );
}

void TestUniformDistribution::testUpperBound()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   _test( distnb.upperBound() == 1 );

   _test( distwb.upperBound() == ub );
}

void TestUniformDistribution::testMean()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   _test( distnb.mean() == 0.5 );

   _test( distwb.mean() == lb+(ub-lb)/2 );
}

void TestUniformDistribution::testStdDev()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   _test( fabs(distnb.stdDev() - sqrt(1/12.0)) < 1.0e-10);
   
   _test( fabs(distwb.stdDev() - sqrt(pow(ub-lb,2)/12.0)) < 1.0e-10 );
}

void TestUniformDistribution::testPrint()
{
   char buf[128];

   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "UNIFORM 0 1" );
#else
   int numberOfChars = sprintf( buf, "UNIFORM 0 1" );
   assert(numberOfChars < 128);
#endif // HAVE_SNPRINTF

   test_str = buf; 
   
   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   fout.open( "TestUniformDistribution_Log" );
   if( fout )
   {
      distnb.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestUniformDistribution_Log" );
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

   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   memset( buf, 0, 128 );
#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "UNIFORM %d %d", (int)lb, (int)ub );
#else
   numberOfChars = sprintf( buf, "UNIFORM %d %d", (int)lb, (int)ub );
   assert(numberOfChars < 128);
#endif // HAVE_SNPRINTF
   test_str = buf; 

   fout.open( "TestUniformDistribution_Log" );
   if( fout )
   {
      distwb.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestUniformDistribution_Log" );
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
  

void TestUniformDistribution::testPrintAttributes()
{
   char buf[128];

   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "distribution=\"uniform\" lower=\"0\" upper=\"1\"");
#else
   int n = sprintf( buf, "distribution=\"uniform\" lower=\"0\" upper=\"1\"");
   assert(n < 128);
#endif // HAVE_SNPRINTF
   test_str = buf;
   
   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   fout.open( "TestUniformDistribution_Log" );
   if( fout )
   {
      distnb.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestUniformDistribution_Log" );
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

   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   memset( buf, 0, 128 );
#ifdef HAVE_SNPRINTF
   snprintf( buf, 128, "distribution=\"uniform\" lower=\"%d\" upper=\"%d\"", (int)lb, (int)ub );
#else
   n = sprintf( buf, "distribution=\"uniform\" lower=\"%d\" upper=\"%d\"", (int)lb, (int)ub );
   assert(n < 128);
#endif // HAVE_SNPRINTF
   test_str = buf;

   fout.open( "TestUniformDistribution_Log" );
   if( fout )
   {
      distwb.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestUniformDistribution_Log" );
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

void TestUniformDistribution::testTypeName()
{
   UniformDistribution distnb;
   UniformDistribution distwb( lb, ub );

   _test( distnb.typeName() == std::string( "UniformDistribution" ) );
   _test( distwb.typeName() == std::string( "UniformDistribution" ) );
}




































