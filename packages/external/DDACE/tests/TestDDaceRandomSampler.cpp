#include "TestDDaceRandomSampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

TestDDaceRandomSampler::TestDDaceRandomSampler() : seed( 779 )
{
	
   /* For unit tests, we do NOT want a random number generator.              */
   /* We want a number generator that pushes out a known sequence of numbers.*/	
   DistributionBase::usePseudoRandom(true);	
		
	
  // set the seed to ensure the data matches
  DistributionBase::setSeed( seed );
 
  // create distributions need by sampler
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );

  // set up the getSamples test data
  int i,j;
  std::vector<double> tmp( 0 );

  
  double   pts[3][3] = { {7.79, 7.8, 7.81},
                         {7.82,  7.83, 7.84},
                         {7.85, 7.86, 7.87} };

  test_data.resize( 0 );

  for( i = 0; i < 3; i++ )
  {
     test_data.push_back( tmp );
     for( j = 0; j < 3; j++ )
     {
        test_data[i].push_back( pts[i][j] );
     }
  }
}

TestDDaceRandomSampler::~TestDDaceRandomSampler()
{
   /* Restore the Random Number Generator */	
   DistributionBase::usePseudoRandom(false);		
	
  // reseed the Random Number Generator
  DistributionBase::setSeed( seed );
}

void TestDDaceRandomSampler::run()
{
   testDDaceRandomSamplerWithDist();
   testDDaceRandomSamplerWithoutDist();
   testGetSamples();
   testClone();
   testPrint();
   testTypeName();
   testGetParameter();
}

void TestDDaceRandomSampler::testDDaceRandomSamplerWithDist()
{  
   // call constructor DDaceRandomSampler(int nSamples, const Array<Distribution>& dist);
   DDaceRandomSampler samplerwd( 3, dists );
   
   _test( samplerwd.nSamples()      == 3 );
   _test( samplerwd.nInputs()       == 3 );
   _test( samplerwd.dist().size() == 3 );

   for( int i = 0; i < (int) samplerwd.dist().size(); i++ )
   {
     _test( dists[i].lowerBound() == samplerwd.dist()[i].lowerBound() );
     _test( dists[i].upperBound() == samplerwd.dist()[i].upperBound() );
   } 
}

void TestDDaceRandomSampler::testDDaceRandomSamplerWithoutDist()
{
   // call constructor DDaceRandomSampler(int nSamples, int nInputs);
   DDaceRandomSampler samplernd( 3, 3 );

   _test( samplernd.nSamples()      == 3 );
   _test( samplernd.nInputs()       == 3 );
   _test( samplernd.dist().size() == 0 );

   for( int i = 0; i < (int) samplernd.dist().size(); i++ )
   {
     _test( dists[i].lowerBound() == samplernd.dist()[i].lowerBound() );
     _test( dists[i].upperBound() == samplernd.dist()[i].upperBound() );
   }
}

void TestDDaceRandomSampler::testGetSamples()
{
  int i,j;
  
  /* reset the seed */
  DistributionBase::setSeed( seed );  

  /**
    * dists[] is an array that contains 3 different uniform distributions.
    * From each distribution, use a random # generator to generate one value.
    */
  DDaceRandomSampler samplerwd( 3, dists );


  DDaceRandomSampler samplernd( 3, 3 );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  // the array of sample data points is passed into "pass"
  //ret = samplerwd.getSamples( pass );
  samplerwd.getSamples( pass );
  ret = pass;


  if( _test( Arrcmp_ad_est( pass, test_data, 0.001) == 0 ) == false )
  {
    for( i = 0; i < 3; i++ )
    {
       cout << "[";
       for( j = 0; j < 3; j++ )
       {
         cout << "(" << pass[i][j] << "|" << test_data[i][j] << "), ";
       }
       cout << "]" << endl;
    }
  }
  _test( Arrcmp_ad_est( ret, test_data, 0.001 ) == 0 );

  //try {
    //ret = samplernd.getSamples( pass );
    //samplernd.getSamples( pass );
    //ret = pass;

    // if here then exception was thrown, record failure
    //_test( false );
  //} catch( ArgumentMisMatchException e ) {
  //} catch (...) {
    
    // if here exception was thrown, record success
    //_test( true );
  //}
}

void TestDDaceRandomSampler::testClone()
{
  DDaceRandomSampler samplerwd( 3, dists );
  DDaceRandomSampler samplernd( 3, 3 );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceRandomSampler* rtn = (DDaceRandomSampler*)samplerwd.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == samplerwd.typeName() );
   _test( rtn->nSamples() == samplerwd.nSamples() );
   _test( rtn->nInputs() == samplerwd.nInputs() );
   
   if( _test( (rtn->dist()).size() == samplerwd.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == samplerwd.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == samplerwd.dist()[i].upperBound() );
     }
   }
      
 
   a = rtn->lowerBounds();
   b = samplerwd.lowerBounds();
   _test( Arrcmp_d( a, b ) == 0 );
   a.resize( 0 );
   b.resize( 0 );

   a = rtn->upperBounds();
   b = samplerwd.upperBounds();
   _test( Arrcmp_d( a, b ) == 0 );

   // clean up the dynamic memory
   delete rtn;

   // call clone and test the results
   rtn = (DDaceRandomSampler*)samplernd.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == samplernd.typeName() );
   _test( rtn->nSamples() == samplernd.nSamples() );
   _test( rtn->nInputs() == samplernd.nInputs() );
   
   if( _test( (rtn->dist()).size() == samplernd.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == samplernd.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == samplernd.dist()[i].upperBound() );
     }
   }
      
 
   a = rtn->lowerBounds();
   b = samplernd.lowerBounds();
   _test( Arrcmp_d( a, b ) == 0 );
   a.resize( 0 );
   b.resize( 0 );

   a = rtn->upperBounds();
   b = samplernd.upperBounds();
   _test( Arrcmp_d( a, b ) == 0 );

   // clean up the dynamic memory
   delete rtn;
}

void TestDDaceRandomSampler::testPrint()
{  
  DDaceRandomSampler samplerwd( 3, dists );
  DDaceRandomSampler samplernd( 3, 3 );

  char    buf[256];
  char    num[16];
  string  data;
  string  test_str;

  // convert seed from 'long' to 'char*'
  itoa( DistributionBase::seed(), num, 16 );
  
  test_str = "<Random samples=\"3\" seed=\"";
  test_str += string( num );
  test_str += "\"/>";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceRandomSamplerWD_Log" );
  if( fout )
  {
     samplerwd.print( fout );
     fout << endl;
     fout.close();
  }
  else
  {
     _test( false );
  }

  fin.open( "TestDDaceRandomSamplerWD_Log" );
  if( fin )
  {
     fin.getline( buf, 255 );
     data = buf;

     _test( data == test_str );

     fin.close();
  }
  else
  {
     _test( false );
  }

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceRandomSamplerND_Log" );
  if( fout )
  {
     samplernd.print( fout );
     fout << endl;
     fout.close();
  }
  else
  {
     _test( false );
  }

  fin.open( "TestDDaceRandomSamplerND_Log" );
  if( fin )
  {
     fin.getline( buf, 255 );
     data = buf;

     _test( data == test_str );

     fin.close();
  }
  else
  {
     _test( false );
  }
}

void TestDDaceRandomSampler::testTypeName()
{
  DDaceRandomSampler samplerwd( 3, dists );
  DDaceRandomSampler samplernd( 3, 3 );

  // call typeName and test the results
  _test( samplerwd.typeName() == std::string( "DDaceRandomSampler" ) );
  _test( samplernd.typeName() == std::string( "DDaceRandomSampler" ) );
}

void TestDDaceRandomSampler::testGetParameter()
{
  DDaceRandomSampler samplerwd( 3, dists );
  DDaceRandomSampler samplernd( 3, 3 );

  // call getParamter, should throw an exception
  try {
     samplerwd.getParameter( std::string( "INPUT" ) );

     // if here then exception was not thrown record failure
     _test(false);
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }

  try {
    samplernd.getParameter( std::string( "INPUT" ) );

    // if here then exception was not thrown, record failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }
}

