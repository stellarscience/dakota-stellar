#include "TestNormalDistribution.h"
#include "Distribution.h"
#include "StdDeviation.h"
#include "Mean.h"
#include "arrcmp.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string.h>

using namespace std;

TestNormalDistribution::TestNormalDistribution()
{
   // set the seed
   seed = 779;
   DistributionBase::setSeed( seed );

   // set the mean, standard deviation and data set
   data.resize( 10 );
   for( int i = 0; i < 10; i++ )
     {
       data[i] = i+1;
     }
   mean = new Mean( 5.5 );
   stdDev = new StdDeviation( data );

   // set lower and upper bounds and number of deviations
   numDevs = 2;
   lb = 5;
   ub = 10;
   lb_est = mean->value() - numDevs*stdDev->value();
   ub_est = mean->value() + numDevs*stdDev->value();




  DistributionBase::setSeed( seed );
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  DistributionBase::setSeed( seed );
  this->deviates[0] = distSigmaWithNumDeviations.getDeviate();
  DistributionBase::setSeed( seed );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  DistributionBase::setSeed( seed );
  this->deviates[1] = distSigmaWithoutNumDeviations.getDeviate();
  DistributionBase::setSeed( seed );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  DistributionBase::setSeed( seed );
  this->deviates[2] = distBoundsWithoutNumDeviations.getDeviate();
  DistributionBase::setSeed( seed );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );
  DistributionBase::setSeed( seed );
  this->deviates[3] = distBoundsWithNumDeviations.getDeviate();
  DistributionBase::setSeed( seed );


}

TestNormalDistribution::~TestNormalDistribution()
{
  // clean up memory
  delete mean;
  delete stdDev;
}

void TestNormalDistribution::run()
{
    testNormalDistributionSigmaWithNumDeviations(); 
    testNormalDistributionSigmaWithoutNumDeviations();
    testNormalDistributionBoundsWithoutNumDeviations();
    testNormalDistributionBoundsWithNumDeviations();
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

void TestNormalDistribution::testNormalDistributionSigmaWithNumDeviations()
{
  // use constructor NormalDistribution(const mean& m, 
  //    		       const stdDeviation& sigma, 
  //			       double numDeviations)
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );

  _test( distSigmaWithNumDeviations.lowerBound() == lb_est );
  _test( distSigmaWithNumDeviations.upperBound() == ub_est );
  _test( distSigmaWithNumDeviations.mean() == mean->value() );
  _test( distSigmaWithNumDeviations.stdDev() == stdDev->value() );
}

void TestNormalDistribution::testNormalDistributionSigmaWithoutNumDeviations()
{
  // use constructor NormalDistribution(const mean& m, 
  //			       const stdDeviation& sigma)
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );

  _test( distSigmaWithoutNumDeviations.lowerBound() == lb_est );
  _test( distSigmaWithoutNumDeviations.upperBound() == ub_est );
  _test( distSigmaWithoutNumDeviations.mean() == mean->value() );
  _test( distSigmaWithoutNumDeviations.stdDev() == stdDev->value() );
}

void TestNormalDistribution::testNormalDistributionBoundsWithoutNumDeviations()
{
  // use constructor NormalDistribution(double lower, double upper)
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );

  _test( distBoundsWithoutNumDeviations.lowerBound() == lb );
  _test( distBoundsWithoutNumDeviations.upperBound() == ub );
  _test( distBoundsWithoutNumDeviations.mean() == 7.5 );
  _test( fabs(distBoundsWithoutNumDeviations.stdDev()) - 1.275533 < 0.001 );
}

void TestNormalDistribution::testNormalDistributionBoundsWithNumDeviations()
{
  // use constructor NormalDistribution(double lower, 
  //			     double upper, 
  //  			     double numDeviations)
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  _test( distBoundsWithNumDeviations.lowerBound() == lb );
  _test( distBoundsWithNumDeviations.upperBound() == ub );
  _test( distBoundsWithNumDeviations.mean() == (lb + ub)/2 );
  _test( distBoundsWithNumDeviations.stdDev() == (ub-lb)/(2*numDevs) );
}
 
void TestNormalDistribution::testClone()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );
   
   // call clone and test results
   NormalDistribution* rtn = (NormalDistribution*)distSigmaWithNumDeviations.clone();

   _test( rtn->lowerBound() == distSigmaWithNumDeviations.lowerBound() );
   _test( rtn->upperBound() == distSigmaWithNumDeviations.upperBound() );
   _test( rtn->mean() == distSigmaWithNumDeviations.mean() );
   _test( rtn->stdDev() == distSigmaWithNumDeviations.stdDev() );
   _test( rtn->typeName() == distSigmaWithNumDeviations.typeName() );

   // clean up memory
   delete rtn;

   // call clone and test results
   rtn = (NormalDistribution*)distSigmaWithoutNumDeviations.clone();

   _test( rtn->lowerBound() == distSigmaWithoutNumDeviations.lowerBound() );
   _test( rtn->upperBound() == distSigmaWithoutNumDeviations.upperBound() );
   _test( rtn->mean() == distSigmaWithoutNumDeviations.mean() );
   _test( rtn->stdDev() == distSigmaWithoutNumDeviations.stdDev() );
   _test( rtn->typeName() == distSigmaWithoutNumDeviations.typeName() );

   // clean up memory
   delete rtn;

   // call clone and test results
   rtn = (NormalDistribution*)distBoundsWithoutNumDeviations.clone();

   _test( rtn->lowerBound() == distBoundsWithoutNumDeviations.lowerBound() );
   _test( rtn->upperBound() == distBoundsWithoutNumDeviations.upperBound() );
   _test( rtn->mean() == distBoundsWithoutNumDeviations.mean() );
   _test( rtn->stdDev() == distBoundsWithoutNumDeviations.stdDev() );
   _test( rtn->typeName() == distBoundsWithoutNumDeviations.typeName() );

   // clean up memory
   delete rtn;

   // call clone and test results
   rtn = (NormalDistribution*)distBoundsWithNumDeviations.clone();

   _test( rtn->lowerBound() ==distBoundsWithNumDeviations.lowerBound() );
   _test( rtn->upperBound() ==distBoundsWithNumDeviations.upperBound() );
   _test( rtn->mean() ==distBoundsWithNumDeviations.mean() );
   _test( rtn->stdDev() ==distBoundsWithNumDeviations.stdDev() );
   _test( rtn->typeName() ==distBoundsWithNumDeviations.typeName() );

   // clean up memory
   delete rtn;
}

void TestNormalDistribution::testGetDeviateNoProb()
{
  DistributionBase::setSeed( seed );
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  DistributionBase::setSeed( seed );
  _test( fabs(distSigmaWithNumDeviations.getDeviate()) - this->deviates[0] < 0.001 );
  DistributionBase::setSeed( seed );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  DistributionBase::setSeed( seed );
   _test( fabs(distSigmaWithoutNumDeviations.getDeviate()) - this->deviates[1] < 0.001 );
  DistributionBase::setSeed( seed );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  DistributionBase::setSeed( seed );
   _test( fabs(distBoundsWithoutNumDeviations.getDeviate()) - this->deviates[2] < 0.001 );
  DistributionBase::setSeed( seed );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );
  DistributionBase::setSeed( seed );
   _test( fabs(distBoundsWithNumDeviations.getDeviate()) - this->deviates[3] < 0.001 );
  DistributionBase::setSeed( seed );

}

void TestNormalDistribution::testGetDeviateWithProb()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  // test the returned deviate value
  _test( distSigmaWithNumDeviations.getDeviate( 0.5 ) == 5.5 );
  _test( distSigmaWithoutNumDeviations.getDeviate( 0.5 ) == 5.5 );
  _test( distBoundsWithoutNumDeviations.getDeviate( 0.5 ) == 7.5 );
  _test( distBoundsWithNumDeviations.getDeviate( 0.5 ) == 7.5 );
}

void TestNormalDistribution::testGetCDF()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

   // test the returned CDF value
   _test( fabs(distSigmaWithNumDeviations.getCDF( 0.5 )) - 0.0568407 < 0.001 );
   _test( fabs(distSigmaWithoutNumDeviations.getCDF( 0.5 )) - 0.0568407 < 0.001 );
   _test( fabs(distBoundsWithoutNumDeviations.getCDF( 0.5 )) - 0.0263157 < 0.001 );
   _test( fabs(distBoundsWithNumDeviations.getCDF( 0.5 )) - 0.0263157 < 0.001 );
}

void TestNormalDistribution::testLowerBound()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  // call lowerBound() and test the results
  _test( distSigmaWithNumDeviations.lowerBound() == (mean->value() - numDevs*(stdDev->value())) );
  _test( distSigmaWithoutNumDeviations.lowerBound() == (mean->value() - 2*(stdDev->value())) );
  _test( distBoundsWithoutNumDeviations.lowerBound() == lb );
  _test( distBoundsWithNumDeviations.lowerBound() == lb );
}

void TestNormalDistribution::testUpperBound()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  // call upperBound() and test the results
  _test( distSigmaWithNumDeviations.upperBound() - (mean->value() + numDevs*(stdDev->value())) < 0.001 );
  _test( distSigmaWithoutNumDeviations.upperBound() - (mean->value() + 2*(stdDev->value())) < 0.001 );
  _test( distBoundsWithoutNumDeviations.upperBound() == ub );
  _test( distBoundsWithNumDeviations.upperBound() == ub );
}

void TestNormalDistribution::testMean()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  // call mean() and tests result
  _test( distSigmaWithNumDeviations.mean() == mean->value() );
  _test( distSigmaWithoutNumDeviations.mean() == mean->value() );
  _test( distBoundsWithoutNumDeviations.mean() == 7.5 );
  _test( distBoundsWithNumDeviations.mean() == (lb + ub)/2 );
}

void TestNormalDistribution::testStdDev()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  // call stdDev() and test results
  _test( distSigmaWithNumDeviations.stdDev() == stdDev->value() );
  _test( distSigmaWithoutNumDeviations.stdDev() == stdDev->value() );
  _test( fabs(distBoundsWithoutNumDeviations.stdDev()) - 1.275533 < 0.001 );
  _test( distBoundsWithNumDeviations.stdDev() == (ub-lb)/(2*numDevs) ); 
}

void TestNormalDistribution::testPrint()
{
   char buf[128];

  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

#ifdef HAVE_SNPRINTF
   snprintf( buf, 127, "NORMAL MEAN %.1f DEV %.5f CUTOFF %d", distSigmaWithNumDeviations.mean(),
                                                            distSigmaWithNumDeviations.stdDev(),
                                                            (int)floor((distSigmaWithNumDeviations.upperBound()-distSigmaWithNumDeviations.lowerBound())/
                                                            2.0/distSigmaWithNumDeviations.stdDev()) );
#else
       int numberOfChars =  sprintf( buf, "NORMAL MEAN %.1f DEV %.5f CUTOFF %d", distSigmaWithNumDeviations.mean(),
                                                            distSigmaWithNumDeviations.stdDev(),
                                                            (int)floor((distSigmaWithNumDeviations.upperBound()-distSigmaWithNumDeviations.lowerBound())/
                                                            2.0/distSigmaWithNumDeviations.stdDev()) );
    assert(numberOfChars<127);
#endif // HAVE_SNPRINTF
    
   test_str = buf; 
   
   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distSigmaWithNumDeviations.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
      fin.close();
   }
   else
   {
      _test( false );
   }

   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   memset( buf, 0, 128 );

   //
   // HERE!!! add formatting info to format string
   //

#ifdef HAVE_SNPRINTF
   snprintf( buf, 127, "NORMAL MEAN %.1f DEV %.5f CUTOFF %d", distSigmaWithoutNumDeviations.mean(),
                                                            distSigmaWithoutNumDeviations.stdDev(),
                                                            (int)floor((distSigmaWithoutNumDeviations.upperBound()-distSigmaWithoutNumDeviations.lowerBound())/
                                                            2.0/distSigmaWithoutNumDeviations.stdDev()) );
#else
   numberOfChars = sprintf( buf, "NORMAL MEAN %.1f DEV %.5f CUTOFF %d", distSigmaWithoutNumDeviations.mean(),
                                                            distSigmaWithoutNumDeviations.stdDev(),
                                                            (int)floor((distSigmaWithoutNumDeviations.upperBound()-distSigmaWithoutNumDeviations.lowerBound())/
                                                            2.0/distSigmaWithoutNumDeviations.stdDev()) );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF

   test_str = buf; 

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distSigmaWithoutNumDeviations.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

    if( _test( data == test_str ) == false )
    {
      cerr << "data:     \"" << data << "\"" << endl
          << "test_str: \"" << test_str << "\"" << endl;
    }
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
   snprintf( buf, 127, "NORMAL MEAN %.1f DEV %.5f CUTOFF %.5f", distBoundsWithoutNumDeviations.mean(),
                                                            distBoundsWithoutNumDeviations.stdDev(),
                                                            (distBoundsWithoutNumDeviations.upperBound()-distBoundsWithoutNumDeviations.lowerBound())/
                                                            2.0/distBoundsWithoutNumDeviations.stdDev() );

#else
   numberOfChars = sprintf( buf, "NORMAL MEAN %.1f DEV %.5f CUTOFF %.5f", distBoundsWithoutNumDeviations.mean(),
                                                            distBoundsWithoutNumDeviations.stdDev(),
                                                            (distBoundsWithoutNumDeviations.upperBound()-distBoundsWithoutNumDeviations.lowerBound())/
                                                            2.0/distBoundsWithoutNumDeviations.stdDev() );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF
   test_str = buf; 

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distBoundsWithoutNumDeviations.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
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
   snprintf( buf, 127, "NORMAL MEAN %.1f DEV %.2f CUTOFF %d", distBoundsWithNumDeviations.mean(),
                                                            distBoundsWithNumDeviations.stdDev(),
                                                            (int)floor((distBoundsWithNumDeviations.upperBound()-distBoundsWithNumDeviations.lowerBound())/
                                                            2.0/distBoundsWithNumDeviations.stdDev()) );

#else
   numberOfChars = sprintf( buf, "NORMAL MEAN %.1f DEV %.2f CUTOFF %d", distBoundsWithNumDeviations.mean(),
                                                            distBoundsWithNumDeviations.stdDev(),
                                                            (int)floor((distBoundsWithNumDeviations.upperBound()-distBoundsWithNumDeviations.lowerBound())/
                                                            2.0/distBoundsWithNumDeviations.stdDev()) );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF

   test_str = buf; 

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distBoundsWithNumDeviations.print( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
      fin.close();
   }
   else
   {
      _test( false );
   }
}
  

void TestNormalDistribution::testPrintAttributes()
{
   char buf[128];

  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

   ofstream fout;
   ifstream fin;

   string data;
   string test_str;

#ifdef HAVE_SNPRINTF
   snprintf( buf, 127, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%d\"", 
                       distSigmaWithNumDeviations.mean(), distSigmaWithNumDeviations.stdDev(),
                       (int)floor((distSigmaWithNumDeviations.upperBound()-distSigmaWithNumDeviations.lowerBound())/
                       2.0/distSigmaWithNumDeviations.stdDev()) );
#else
   int numberOfChars = sprintf( buf, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%d\"", 
                       distSigmaWithNumDeviations.mean(), distSigmaWithNumDeviations.stdDev(),
                       (int)floor((distSigmaWithNumDeviations.upperBound()-distSigmaWithNumDeviations.lowerBound())/
                       2.0/distSigmaWithNumDeviations.stdDev()) );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF

   test_str = buf;
   
   // write the data to a file, then reopen the
   // file and verify the data was written correctly
   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distSigmaWithNumDeviations.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
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
   snprintf( buf, 127, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%d\"", 
                       distSigmaWithoutNumDeviations.mean(), distSigmaWithoutNumDeviations.stdDev(),
                       (int)floor((distSigmaWithoutNumDeviations.upperBound()-distSigmaWithoutNumDeviations.lowerBound())/
                       2.0/distSigmaWithoutNumDeviations.stdDev()) );
#else
   numberOfChars = sprintf( buf, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%d\"", 
                       distSigmaWithoutNumDeviations.mean(), distSigmaWithoutNumDeviations.stdDev(),
                       (int)floor((distSigmaWithoutNumDeviations.upperBound()-distSigmaWithoutNumDeviations.lowerBound())/
                       2.0/distSigmaWithoutNumDeviations.stdDev()) );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF
   test_str = buf;

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distSigmaWithoutNumDeviations.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
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
   snprintf( buf, 127, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%.5f\"", 
                       distBoundsWithoutNumDeviations.mean(), distBoundsWithoutNumDeviations.stdDev(),
                       (distBoundsWithoutNumDeviations.upperBound()-distBoundsWithoutNumDeviations.lowerBound())/
                       2.0/distBoundsWithoutNumDeviations.stdDev() );
#else
   numberOfChars = sprintf( buf, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.5f\" cutoff=\"%.5f\"", 
                       distBoundsWithoutNumDeviations.mean(), distBoundsWithoutNumDeviations.stdDev(),
                       (distBoundsWithoutNumDeviations.upperBound()-distBoundsWithoutNumDeviations.lowerBound())/
                       2.0/distBoundsWithoutNumDeviations.stdDev() );
   assert(numberOfChars<127);
#endif // HAVE_SNPRINTF
   test_str = buf;

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distBoundsWithoutNumDeviations.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
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
   snprintf( buf, 127, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.2f\" cutoff=\"%d\"", 
                       distBoundsWithNumDeviations.mean(), distBoundsWithNumDeviations.stdDev(),
                       (int)floor((distBoundsWithNumDeviations.upperBound()-distBoundsWithNumDeviations.lowerBound())/
                       2.0/distBoundsWithNumDeviations.stdDev()) );
#else
   numberOfChars = sprintf( buf, "distribution=\"normal\" mean=\"%.1f\" sigma=\"%.2f\" cutoff=\"%d\"", 
                       distBoundsWithNumDeviations.mean(), distBoundsWithNumDeviations.stdDev(),
                       (int)floor((distBoundsWithNumDeviations.upperBound()-distBoundsWithNumDeviations.lowerBound())/
                       2.0/distBoundsWithNumDeviations.stdDev()) );
    assert(numberOfChars<127);
#endif // HAVE_SNPRINTF
   test_str = buf;

   fout.open( "TestNormalDistribution_Log" );
   if( fout )
   {
      distBoundsWithNumDeviations.printAttributes( fout );
      fout << endl;
      fout.close();
   }
   else
   {
      _test( false );
   }

   fin.open( "TestNormalDistribution_Log" );
   if( fin )
   {
      fin.getline( buf, 128 );
      data = buf;

      if( _test( data == test_str ) == false )
      {
         cerr << "data:     \"" << data << "\"" << endl
              << "test_str: \"" << test_str << "\"" << endl;
      }
      fin.close();
   }
   else
   {
      _test( false );
   }
}

void TestNormalDistribution::testTypeName()
{
  NormalDistribution distSigmaWithNumDeviations( *mean, *stdDev, numDevs );
  NormalDistribution distSigmaWithoutNumDeviations( *mean, *stdDev );
  NormalDistribution distBoundsWithoutNumDeviations( lb, ub );
  NormalDistribution distBoundsWithNumDeviations( lb, ub, numDevs );

  std::string normal( "NormalDistribution" );

   _test( distSigmaWithNumDeviations.typeName() == normal );
   _test( distSigmaWithoutNumDeviations.typeName() == normal );
   _test( distBoundsWithoutNumDeviations.typeName() == normal );
   _test( distBoundsWithNumDeviations.typeName() == normal );
}

