#include "TestDDaceCentralCompositeSampler.h"
//#include "ArgumentMisMatchException.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

TestDDaceCentralCompositeSampler::TestDDaceCentralCompositeSampler()
{
  // create distributions need by sampler
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 6, 7 ) ) );
  dists.push_back( Distribution( UniformDistribution( 7, 9 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 5 ) ) );

  // set up the getSamples test data
  int i,j;
  std::vector<double> tmp( 0 );
  double   pts[15][3] = { {6.5,8,2.5},
                          {6,8,2.5},
                          {7,8,2.5},
                          {6.5,7,2.5},
                          {6.5,9,2.5},
                          {6.5,8,0},
                          {6.5,8,5},
                          {6,7,0},
                          {7,7,0},
                          {6,9,0},
                          {7,9,0},
                          {6,7,5},
                          {7,7,5},
                          {6,9,5},
                          {7,9,5} };
                          
  test_data.resize( 0 );

  for( i = 0; i < 15; i++ )
  {
     test_data.push_back( tmp );               
     for( j = 0; j < 3; j++ )
     {
        test_data[i].push_back( pts[i][j] );
                                           
     }                                      
  }                                         
}

TestDDaceCentralCompositeSampler::~TestDDaceCentralCompositeSampler()
{
}

void TestDDaceCentralCompositeSampler::run()
{
   testDDaceCentralCompositeSampler();
   testGetSamples();
   testClone();
   testPrint();
   testTypeName();
}

void TestDDaceCentralCompositeSampler::testDDaceCentralCompositeSampler()
{
    // create the sampler
    try {
      DDaceCentralCompositeSampler sampler( 15, 3, dists );
      
      // if here then no exception was thrown, record success
      _test( true );
      _test( sampler.nSamples() == 15 );
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
    }
    catch (...) {
    	_test(false);
    }
}

void TestDDaceCentralCompositeSampler::testGetSamples()
{
  DDaceCentralCompositeSampler sampler( 15, 3, dists );
  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = sampler.getSamples( pass );
  sampler.getSamples( pass );
  ret = pass;

  _test( Arrcmp_ad( pass, test_data ) == 0 );
  _test( Arrcmp_ad( ret, test_data ) == 0 );
}

void TestDDaceCentralCompositeSampler::testClone()
{
   DDaceCentralCompositeSampler sampler( 15, 3, dists );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceCentralCompositeSampler* rtn = (DDaceCentralCompositeSampler*)sampler.clone();

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

void TestDDaceCentralCompositeSampler::testPrint()
{
  DDaceCentralCompositeSampler sampler( 15, 3, dists );
  char    buf[256];
  string  data;
  string  test_str[2];
  
  test_str[0] = "METHOD Central Composite Design";
  test_str[1] = "SAMPLES 15";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceCentralCompositeSampler_Log" );
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

  fin.open( "TestDDaceCentralCompositeSampler_Log" );
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

void TestDDaceCentralCompositeSampler::testTypeName()
{
  DDaceCentralCompositeSampler sampler( 15, 3, dists );

  // call type name and check the results
  _test( sampler.typeName() == "DDaceCentralCompositeSampler" );
}

