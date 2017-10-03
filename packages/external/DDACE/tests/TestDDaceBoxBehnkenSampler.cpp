#include "TestDDaceBoxBehnkenSampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "UniformDistribution.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

TestDDaceBoxBehnkenSampler::TestDDaceBoxBehnkenSampler()
{
  // create distributions need by sampler
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 6, 7 ) ) );
  dists.push_back( Distribution( UniformDistribution( 7, 9 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 5 ) ) );

  // set up the getSamples test data
  int i,j;
  std::vector<double> tmp( 0 );
  double   pts[13][3] = { {6.5,8,2.5},  // these are the correct sample points
                          {7,9,2.5},
                          {7,7,2.5},
                          {6,9,2.5},
                          {6,7,2.5},
                          {7,8,5},
                          {7,8,0},
                          {6,8,5},
                          {6,8,0},
                          {6.5,9,5},
                          {6.5,9,0},
                          {6.5,7,5},
                          {6.5,7,0} };
                          
  test_data.resize( 0 );

  for( i = 0; i < 13; i++ )
  {
     test_data.push_back( tmp );
     for( j = 0; j < 3; j++ )
     {
        test_data[i].push_back( pts[i][j] );
     }
  }
}

TestDDaceBoxBehnkenSampler::~TestDDaceBoxBehnkenSampler()
{
}

void TestDDaceBoxBehnkenSampler::run()
{
   testDDaceBoxBehnkenSampler();
   testGetSamples();
   testClone();
   testPrint();
   testTypeName();
}

void TestDDaceBoxBehnkenSampler::testDDaceBoxBehnkenSampler()
{
    // create the sampler
    try {
      DDaceBoxBehnkenSampler sampler( 13, 3, dists );
      
      // if here then no exception was thrown, record success
      _test( true );
      _test( sampler.nSamples() == 13 );
      _test( sampler.nInputs() == 3 );
      _test( sampler.dist().size() == 3 );

      for( int i = 0; i < (int) sampler.dist().size(); i++ )
      {
         _test( dists[i].lowerBound() == sampler.dist()[i].lowerBound() );
         _test( dists[i].upperBound() == sampler.dist()[i].upperBound() );
      }

    } catch( std::exception& e ) {      
      e.what();
      // if here then an exception was thrown, record failure
      _test( false );
    
    } catch (...) {
    	_test(false);
    }
}

void TestDDaceBoxBehnkenSampler::testGetSamples()
{
  DDaceBoxBehnkenSampler sampler( 13, 3, dists );
  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = sampler.getSamples( pass );
  sampler.getSamples( pass );
  ret = pass;

  _test( Arrcmp_ad( pass, test_data ) == 0 );
  _test( Arrcmp_ad( ret, test_data ) == 0 );
}

void TestDDaceBoxBehnkenSampler::testClone()
{
   DDaceBoxBehnkenSampler sampler( 13, 3, dists );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceBoxBehnkenSampler* rtn = (DDaceBoxBehnkenSampler*)sampler.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == sampler.typeName() );
   _test( rtn->nSamples() == sampler.nSamples() );
   _test( rtn->nInputs() == sampler.nInputs() );
   _test( rtn->noise() == sampler.noise() );
   
   if( _test( (rtn->dist()).size() == sampler.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == sampler.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == sampler.dist()[i].upperBound() );
     }
   }
      
 
   a = rtn->lowerBounds();
   b = sampler.lowerBounds();
   _test( Arrcmp_d( a, b ) == 0 );
   a.resize( 0 );
   b.resize( 0 );

   a = rtn->upperBounds();
   b = sampler.upperBounds();
   _test( Arrcmp_d( a, b ) == 0 );

   // clean up the dynamic memory
   delete rtn;
}

void TestDDaceBoxBehnkenSampler::testPrint()
{
  DDaceBoxBehnkenSampler sampler( 13, 3, dists );
  char    buf[256];
  string  data;
  string  test_str[2];
  
  test_str[0] = "METHOD BoxBehnken";
  test_str[1] = "SAMPLES 13";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceBoxBehnkenSampler_Log" );
  if( fout )
  {
     sampler.print( fout );
     fout << endl;
     fout.close();
  }
  else
  {
     _test( false );
  }

  fin.open( "TestDDaceBoxBehnkenSampler_Log" );
  if( fin )
  {
     for( int i = 0; i < 2; i++ )
     {
       fin.getline( buf, 255 );
       data = buf;

       _test( data == test_str[i] );
     }

     fin.close();
  }
  else
  {
     _test( false );
  }
}

void TestDDaceBoxBehnkenSampler::testTypeName()
{
  DDaceBoxBehnkenSampler sampler( 13, 3, dists );

  // call type name and check the results
  _test( sampler.typeName() == "DDaceBoxBehnkenSampler" );
}

