#include "TestDDaceOASampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

TestDDaceOASampler::TestDDaceOASampler() : seed( 779 )
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
  
  
  
  double   pts[4][3] = { {2.5,2.5,2.5},
                         {2.5,7.5,7.5},
                         {7.5,2.5,7.5},
                         {7.5,7.5,2.5} };
                        
                         
                         
                         

  double pts_wn[4][3] = { {3.895,3.9,3.905},
                          {3.91, 8.915,8.92},
                          {8.925,3.93,8.935},
                          {8.94,8.945,3.95} };

  DistributionBase::setSeed( seed );  
                          
  test_data.resize( 0 );
  test_data_wn.resize( 0 );

  for( i = 0; i < 4; i++ )
  {
     test_data.push_back( tmp );
     test_data_wn.push_back( tmp );
     for( j = 0; j < 3; j++ )
     {
        test_data[i].push_back( pts[i][j] );
        test_data_wn[i].push_back( pts_wn[i][j] );
     }
  }
}

TestDDaceOASampler::~TestDDaceOASampler()
{
   /* Restore the Random Number Generator */	
   DistributionBase::usePseudoRandom(false);		
	
  // reseed the Random Number Generator
  DistributionBase::setSeed( seed );
}

void TestDDaceOASampler::run()
{
   testDDaceOASamplerWithDist();
   testDDaceOASamplerWithoutDist();
   testGetSamplesWithoutNoise();
   testGetSamplesWithNoise();
   testClone();
   testPrint();
   testTypeName();
   testGetParameter();
}

void TestDDaceOASampler::testDDaceOASamplerWithDist()
{  
   // call constructor DDaceOASampler(int nSamples, bool noise, const Array<Distribution>& dist);
   DDaceOASampler samplerwd( 4, false, dists );

   _test( samplerwd.nSamples()      == 4 );
   _test( samplerwd.nInputs()       == 3 );
   //_test( samplerwd.nSymbols()      == 2 );
   _test( samplerwd.noise()         == false );
   _test( samplerwd.dist().size() == 3 );

   for( int i = 0; i < (int) samplerwd.dist().size(); i++ )
   {
     _test( dists[i].lowerBound() == samplerwd.dist()[i].lowerBound() );
     _test( dists[i].upperBound() == samplerwd.dist()[i].upperBound() );
   } 
}

void TestDDaceOASampler::testDDaceOASamplerWithoutDist()
{
   // call constructor DDaceOASampler(int nSamples, int nInputs, bool noise);
   DDaceOASampler samplernd( 4, 3, false );

   _test( samplernd.nSamples()      == 4 );
   _test( samplernd.nInputs()       == 3 );
  // _test( samplernd.nSymbols()      == 2 );
   _test( samplernd.noise()         == false );
   _test( samplernd.dist().size() == 0 );

   for( int i = 0; i < (int) samplernd.dist().size(); i++ )
   {
     _test( dists[i].lowerBound() == samplernd.dist()[i].lowerBound() );
     _test( dists[i].upperBound() == samplernd.dist()[i].upperBound() );
   }
}

void TestDDaceOASampler::testGetSamplesWithoutNoise()
{
  int i,j;
  
  DistributionBase::setSeed( seed );
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) ); 


   DDaceOASampler samplerwd( 4, false, dists );
   DDaceOASampler samplernd( 4, 3, false );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = samplerwd.getSamples( pass );
  samplerwd.getSamples( pass );
  ret = pass;

  if( _test( Arrcmp_ad_est( pass, test_data, 0.001 ) == 0 ) == false )
  {
    for( i = 0; i < 4; i++ )
    {
       cout << "[";
       for( j = 0; j < 3; j++ )
       {
         cout << "(" << pass[i][j] << "|" << test_data[i][j] << "), ";
       }
       cout << "]" << endl;
    }
  }
  _test( Arrcmp_ad( ret, test_data ) == 0 );

//  try {
//    //ret = samplernd.getSamples( pass );
//    samplernd.getSamples( pass );
//    ret = pass;
//
//    // if here then exception was thrown, record failure
//    _test( false );
//  //} catch( ArgumentMisMatchException e ) {
//  } catch (...) {
//    
//    // if here exception was thrown, record success
//    _test( true );
//  }
}

void TestDDaceOASampler::testGetSamplesWithNoise()
{
  int i,j;
  
  DistributionBase::setSeed( seed );    
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 10 ) ) ); 

  DistributionBase::setSeed( seed );                         
   DDaceOASampler samplerwd( 4, true, dists );
   DDaceOASampler samplernd( 4, 3, true );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = samplerwd.getSamples( pass );
  DistributionBase::setSeed( seed );
  samplerwd.getSamples( pass );
  ret = pass;

  if( _test( Arrcmp_ad_est( pass, test_data_wn, 0.001 ) == 0 ) == false )
  {
    for( i = 0; i < 4; i++ )
    {
       cout << "[";
       for( j = 0; j < 3; j++ )
       {
         cout << "(" << pass[i][j] << "|" << test_data_wn[i][j] << "), ";
       }
       cout << "]" << endl;
    }
  }
  _test( Arrcmp_ad_est( ret, test_data_wn, 0.001 ) == 0 );

//  try {
//    //ret = samplernd.getSamples( pass );
//    samplernd.getSamples( pass );
//    ret = pass;
//
//    // if here then exception was thrown, record failure
//    _test( false );
//  //} catch( ArgumentMisMatchException e ) {
//  } catch (...) {
//    
//    // if here exception was thrown, record success
//    _test( true );
//  }
}

void TestDDaceOASampler::testClone()
{
   DDaceOASampler samplerwd( 4, false, dists );
   DDaceOASampler samplernd( 4, 3, false );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceOASampler* rtn = (DDaceOASampler*)samplerwd.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == samplerwd.typeName() );
   _test( rtn->nSamples() == samplerwd.nSamples() );
   _test( rtn->nInputs() == samplerwd.nInputs() );
   //_test( rtn->nSymbols() == samplerwd.nSymbols() );
   _test( rtn->noise() == samplerwd.noise() );
   
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
   rtn = (DDaceOASampler*)samplernd.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == samplernd.typeName() );
   _test( rtn->nSamples() == samplernd.nSamples() );
   _test( rtn->nInputs() == samplernd.nInputs() );
   //_test( rtn->nSymbols() == samplernd.nSymbols() );
   _test( rtn->noise() == samplernd.noise() );
   
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

void TestDDaceOASampler::testPrint()
{  
   DDaceOASampler samplerwd( 4, false, dists );
   DDaceOASampler samplernd( 4, 3, false );

  char    buf[256];
  char    num[16];
  string  data;
  string  test_str;

  // convert seed from 'long' to 'char*'
  itoa( DistributionBase::seed(), num, 16 );
  
  test_str = "<OrthogonalArray samples=\"4\" perturb=\"false\" seed=\"";
  test_str += string( num );
  test_str += "\"/>";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceOASamplerWD_Log" );
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

  fin.open( "TestDDaceOASamplerWD_Log" );
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
  
  fout.open( "TestDDaceOASamplerND_Log" );
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

  fin.open( "TestDDaceOASamplerND_Log" );
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

void TestDDaceOASampler::testTypeName()
{
   DDaceOASampler samplerwd( 4, false, dists );
   DDaceOASampler samplernd( 4, 3, false );

  // call typeName and test the results
  _test( samplerwd.typeName() == std::string( "DDaceOASampler" ) );
  _test( samplernd.typeName() == std::string( "DDaceOASampler" ) );
}

void TestDDaceOASampler::testGetParameter()
{
   DDaceOASampler samplerwd( 3, false, dists );
   DDaceOASampler samplernd( 3, 3, false );

  // call getParameters, should throw exception
  try {
     samplerwd.getParameter( std::string( "PARAM_NAME" ) );

     // if here then no exception was thrown record failure
     _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }

  try {
    samplernd.getParameter( std::string( "PARAM_NAME" ) );

    // if here then no exception was thrown record failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }
}

