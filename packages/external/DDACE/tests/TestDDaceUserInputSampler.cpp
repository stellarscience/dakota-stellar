#include "TestDDaceUserInputSampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>


using namespace std;

TestDDaceUserInputSampler::TestDDaceUserInputSampler()
{
  // set up the getSamples test data
  int i,j;
  std::vector<double> tmp( 0 );
  double   pts[3][3] = { {2, 5, 9},
                         {8, 2, 4},
                         {6, 6, 3} };
                          
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

TestDDaceUserInputSampler::~TestDDaceUserInputSampler()
{
}

void TestDDaceUserInputSampler::run()
{
   testDDaceUserInputSampler();
   testGetSamples();
   testClone();
   testPrint();
   testTypeName();
   testGetParameter();
}

void TestDDaceUserInputSampler::testDDaceUserInputSampler()
{  
   // call constructor DDaceUserInputSampler(const String& ptFilename);
   DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );
   
   _test( sampler.nSamples() == 3 );
   _test( sampler.nInputs() == 3 ); 
}

void TestDDaceUserInputSampler::testGetSamples()
{
  int i,j;

  DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );

  std::vector<DDaceSamplePoint> pass( 0 );
  std::vector<DDaceSamplePoint> ret( 0 );

  // call getSamples and test the results
  //ret = sampler.getSamples( pass );
  sampler.getSamples( pass );
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
}

void TestDDaceUserInputSampler::testClone()
{
   DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );
   std::vector<double> a;
   std::vector<double> b;

   // call clone and test the results
   DDaceUserInputSampler* rtn = (DDaceUserInputSampler*)sampler.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == sampler.typeName() );
   _test( rtn->nSamples() == sampler.nSamples() );
   _test( rtn->nInputs() == sampler.nInputs() );
 
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

void TestDDaceUserInputSampler::testPrint()
{  
  DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );

  char    buf[256];
  char    num[16];
  string  data;
  string  test_str;

  // convert seed from 'long' to 'char*'
  itoa( sampler.nSamples(), num, 16 );
  
  test_str = "<UserInputSampler filename=\"TestDDaceUserInputSamplerData\" samples=\"";
  test_str += string( num );
  test_str += "\"/>";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceUserInputsampler_Log" );
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

  fin.open( "TestDDaceUserInputsampler_Log" );
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

void TestDDaceUserInputSampler::testTypeName()
{
  DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );

  // call typeName and test the results
  _test( sampler.typeName() == std::string( "DDaceUserInputSampler" ) );
}

void TestDDaceUserInputSampler::testGetParameter()
{
  DDaceUserInputSampler sampler( std::string( "TestDDaceUserInputSamplerData" ) );

  // call getParamter, should throw an exception
  try {
     sampler.getParameter( std::string( "INPUT" ) );

     // if here then exception was not thrown record failure
     _test(false);
  } catch( std::exception& e ) {
    // if here then exception was thrown record success
    _test( true );
  }
}

