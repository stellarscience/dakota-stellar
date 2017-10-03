#include "TestDDaceFactorialSampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

TestDDaceFactorialSampler::TestDDaceFactorialSampler() : seed( 779 )
{
	
   /* For unit tests, we do NOT want a random number generator.              */
   /* We want a number generator that pushes out a known sequence of numbers.*/	
   DistributionBase::usePseudoRandom(true);	
	
  // set random number generator seed
  DistributionBase::setSeed( seed );

  // set up the getSamples test data
  int i,j;
  std::vector<double> tmp( 0 );
  double   pts[16][2] = { {12.5,12.5},  // correct data for nSamples=16, nSymbols=4, nInputs=2
                          {37.5,12.5},
                          {62.5,12.5},
                          {87.5,12.5},
                          {12.5,37.5},
                          {37.5,37.5},
                          {62.5,37.5},
                          {87.5,37.5},
                          {12.5,62.5},
                          {37.5,62.5},
                          {62.5,62.5},
                          {87.5,62.5},
                          {12.5,87.5},
                          {37.5,87.5},
                          {62.5,87.5},
                          {87.5,87.5} };

  double pts_wn[16][2] = {{19.475,  19.5},
                          {44.525, 19.55},
                          {69.575, 19.6},
                          {94.625,19.65},
                          {19.675,44.7},
                          {44.725,44.75},
                          {69.775,44.8},
                          {94.825,44.85},
                          {19.875,69.9},
                          {44.925,69.95},
                          {69.975,70},
                          {95.025,70.05},
                          {20.075,95.1},
                          {45.125,95.15},
                          {70.175,95.2},
                          {95.225,95.25}};

  DistributionBase::setSeed( seed );

  // create distributions need by sampler
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );  
  std::vector<DDaceSamplePoint> data( 0 );
  DistributionBase::setSeed( seed );

  
                          
  test_data.resize( 0 );
  test_data_wn.resize( 0 );

  for( i = 0; i < 16; i++ )
  {
     test_data.push_back( tmp );
     test_data_wn.push_back( tmp );
     for( j = 0; j < 2; j++ )
     {
        test_data[i].push_back( pts[i][j] );
        test_data_wn[i].push_back( pts_wn[i][j] );
     }
  }
}

TestDDaceFactorialSampler::~TestDDaceFactorialSampler()
{
   /* Restore the Random Number Generator */	
   DistributionBase::usePseudoRandom(false);		
	
  // reseed the Random Number Generator
  DistributionBase::setSeed( seed );
}

void TestDDaceFactorialSampler::run()
{
   testDDaceFactorialSampler2();
   testDDaceFactorialSampler4();
   testGetSamplesWithoutNoise();
   testGetSamplesWithNoise();
   testClone();
   testPrint();
   testTypeName();
   testGetParameter();
}

void TestDDaceFactorialSampler::testDDaceFactorialSampler4()
{  
   // create a factorial sampler using the constructor 
   // DDaceFactorialSampler(int nSamples, int nSymbols, bool noise,
   //                       const Array<Distribution>& dist);
   try {
     DDaceFactorialSampler sampler4( 16, 4, false, dists );
     
     // if here then no exception was thrown, record success
     _test( true );
     _test( sampler4.nSamples()      == 16 );
     _test( sampler4.nInputs()       == 2 );
     //_test( sampler4.nSymbols()      == 4 );
     _test( sampler4.noise()         == false );
     _test( sampler4.dist().size() == 2 );

     for( int i = 0; i < (int) sampler4.dist().size(); i++ )
     {
       _test( dists[i].lowerBound() == sampler4.dist()[i].lowerBound() );
       _test( dists[i].upperBound() == sampler4.dist()[i].upperBound() );
     } 
   } catch( std::exception& e ) {
    
     // if here then exception was thrown, record failure
     _test( false );
   }
}

void TestDDaceFactorialSampler::testDDaceFactorialSampler2()
{
  // create a factorial sampler using the constructor
  // DDaceFactorialSampler(int nSamples, int nSymbols);
  try {
     DDaceFactorialSampler sampler2( 16, 4 );

     // if here then no exception was thrown, record success
     _test( true );
     _test( sampler2.nSamples()      == 16 );
     _test( sampler2.nInputs()       == 2 );
     //_test( sampler2.nSymbols()      == 4 );
     _test( sampler2.noise()         == false );
     _test( sampler2.dist().size() == 0 );

     for( int i = 0; i < (int) sampler2.dist().size(); i++ )
     {
       _test( dists[i].lowerBound() == sampler2.dist()[i].lowerBound() );
       _test( dists[i].upperBound() == sampler2.dist()[i].upperBound() );
     } 
  } catch( std::exception& e ) {
    
    // if here then exception was thrown, record failure
    _test( false );
  }
}

void TestDDaceFactorialSampler::testGetSamplesWithoutNoise()
{
  DistributionBase::setSeed( seed );
  DDaceFactorialSampler sampler4( 16, 4, false, dists );
  DDaceFactorialSampler sampler2( 16, 4 );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = sampler4.getSamples( pass );
  DistributionBase::setSeed( seed );  
  sampler4.getSamples( pass );
  ret = pass;

  _test( Arrcmp_ad( pass, test_data ) == 0 );
  _test( Arrcmp_ad( ret, test_data ) == 0 );

  try {
    //ret = sampler2.getSamples( pass );
    sampler2.getSamples( pass );
    ret = pass;

    // if here then exception was thrown, record failure
    _test( false );    
  } catch( std::exception& e ) {

    // if here exception was thrown, record success
    _test( true );

  } catch(...) {
    // if here, an exception was throw that was not a std::exception
    // hence record failure.
    _test( false );	 
  }
}






void TestDDaceFactorialSampler::testGetSamplesWithNoise()
{
  int i,j;

  DistributionBase::setSeed( seed );
  DDaceFactorialSampler sampler4( 16, 4, true, dists );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = sampler4.getSamples( pass );
  DistributionBase::setSeed( seed );  
  sampler4.getSamples( pass );
  ret = pass;

  if( _test( Arrcmp_ad_est( pass, test_data_wn, 0.001 ) == 0 ) == false )
  {
    for( i = 0; i < 16; i++ )
    {
       cout << "[";
       for( j = 0; j < 2; j++ )
       {
         cout << "(" << pass[i][j] << "|" << test_data_wn[i][j] << "), ";
       }
       cout << "]" << endl;
    }
  }
  _test( Arrcmp_ad_est( ret, test_data_wn, 0.001 ) == 0 );
}

void TestDDaceFactorialSampler::testClone()
{
   DDaceFactorialSampler sampler4( 16, 4, false, dists );
   DDaceFactorialSampler sampler2( 16, 4 );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceFactorialSampler* rtn = (DDaceFactorialSampler*)sampler4.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == sampler4.typeName() );
   _test( rtn->nSamples() == sampler4.nSamples() );
   //_test( rtn->nSymbols() == sampler4.nSymbols() );
   _test( rtn->nInputs() == sampler4.nInputs() );
   _test( rtn->noise() == sampler4.noise() );
   
   if( _test( (rtn->dist()).size() == sampler4.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == sampler4.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == sampler4.dist()[i].upperBound() );
     }
   }
      
 
   a = rtn->lowerBounds();
   b = sampler4.lowerBounds();
   _test( Arrcmp_d( a, b ) == 0 );
   a.resize( 0 );
   b.resize( 0 );

   a = rtn->upperBounds();
   b = sampler4.upperBounds();
   _test( Arrcmp_d( a, b ) == 0 );

   // clean up the dynamic memory
   delete rtn;

   // call clone and test the results
   rtn = (DDaceFactorialSampler*)sampler2.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == sampler2.typeName() );
   _test( rtn->nSamples() == sampler2.nSamples() );
   //_test( rtn->nSymbols() == sampler2.nSymbols() );
   _test( rtn->nInputs() == sampler2.nInputs() );
   _test( rtn->noise() == sampler2.noise() );
   
   if( _test( (rtn->dist()).size() == sampler2.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == sampler2.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == sampler2.dist()[i].upperBound() );
     }
   }
      
 
   a = rtn->lowerBounds();
   b = sampler2.lowerBounds();
   _test( Arrcmp_d( a, b ) == 0 );
   a.resize( 0 );
   b.resize( 0 );

   a = rtn->upperBounds();
   b = sampler2.upperBounds();
   _test( Arrcmp_d( a, b ) == 0 );

   // clean up the dynamic memory
   delete rtn;
}

void TestDDaceFactorialSampler::testPrint()
{  
  DDaceFactorialSampler sampler4( 16, 4, false, dists );
  DDaceFactorialSampler sampler2( 16, 4 );
  char    buf[256];
  char    num[16];
  string  data;
  string  test_str;

  // convert seed from 'long' to 'char*'
  itoa( DistributionBase::seed(), num, 16 );
  
  test_str = "<Factorial samples=\"16\" symbols=\"4\" perturb=\"false\" seed=\"";
  test_str += string( num );
  test_str += "\"/>";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceFactorialSampler4_Log" );
  if( fout )
  {
     sampler4.print( fout );
     fout << endl;
     fout.close();
  }
  else
  {
     _test( false );
  }

  fin.open( "TestDDaceFactorialSampler4_Log" );
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
  
  fout.open( "TestDDaceFactorialSampler2_Log" );
  if( fout )
  {
     sampler2.print( fout );
     fout << endl;
     fout.close();
  }
  else
  { 
     _test( false );
  }

  fin.open( "TestDDaceFactorialSampler2_Log" );
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

void TestDDaceFactorialSampler::testTypeName()
{
  DDaceFactorialSampler sampler4( 16, 4, false, dists );
  DDaceFactorialSampler sampler2( 16, 4 );

  // call type name and check the results
  _test( sampler4.typeName() == "DDaceFactorialSampler" );
  _test( sampler2.typeName() == "DDaceFactorialSampler" );
}

void TestDDaceFactorialSampler::testGetParameter()
{
  DDaceFactorialSampler sampler4( 16, 4, false, dists );
  DDaceFactorialSampler sampler2( 16, 4 );

  // call getParameters and test the output
  try {
     _test( sampler4.getParameter( std::string( "SYMBOLS" ) ) == 4 );
  } catch( std::exception& e ) {
    // if here then exception was thrown record failure
    _test( false );
  }

  try {
    _test( sampler2.getParameter( std::string( "SYMBOLS" ) ) == 4 );
  } catch( std::exception& e ) {
    // if here then exception was thrown record failure
    _test( false );
  }

  try {
     sampler4.getParameter( std::string( "WRONG_INPUT" ) );

     // if here then exception was not thrown record failure
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }

  try {
    sampler2.getParameter( std::string( "WRONG_INPUT" ) );

    // if here then exception was not thrown, record failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }
}

