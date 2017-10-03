#include "TestDDaceSampler.h"
//#include "DDaceArraySampler.h"
#include "UniformDistribution.h"
#include "DDaceFactorialSampler.h"
#include <stdlib.h>
#include <string.h>
#include <string>

using namespace std;

int Arrcmp_sp( const std::vector<DDaceSamplePoint> & a,
		 const std::vector<DDaceSamplePoint> & b )
{
  int i,j;
  int len_a = a.size(); 
  int len_b = b.size();

  if( len_a != len_b ) {
    return 1;    // failed! the arrays don't match
  }

  for( i = 0; i < len_a; i++ ) {
    if( a[i].length() != b[i].length() ) {
      return 1;  // failed! the arrays don't match
    }
    for( j = 0; j < a[i].length(); j++ ) {
      if( a[i][j] != b[i][j] ) {
        return 1;  // failed! the arrays don't match
      }
    }
  }

  return 0;     // success! the arrays match
}

int Arrcmp_d( const std::vector<double> & a, const std::vector<double> & b )
{
  int i;
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b ) {
    return 1;   // failed! the arrays don't match
  }

  for( i = 0; i < len_a; i++ ) {
     if( a[i] != b[i] ) {
       return 1;   // failed! the arrays don't match
     }
  }

  return 0;     // success! the arrays match
}


TestDDaceSampler::TestDDaceSampler() : 
    data_d(0), data_sp(0), lb(0), ub(0)
{
   // int i,j;

   std::vector<double>     tmp( 0 );
/**   
   double num[4][4] = { { 6, 6, 7, 6 },
                        { 2, 0, 2, 5 },
                        { 2, 4, 0, 5 },
                        { 4, 4, 5, 7 } };
   double lowb[4] =     { 2, 0, 0, 5 };
   double upb[4]  =     { 6, 6, 7, 7 };
**/
   // put numbers from above into the bounds arrays
   // and load numbers into the data arrays
//   for( i = 0; i < 4; i++ )
//   {
//      data_d.append( tmp );
//      for( j = 0; j < 4; j++ )
//      {
//         data_d[i].append( num[i][j] );
//      }
//
//      data_sp.append( DDaceSamplePoint( 4, data_d[i] ) );
//
//      lb.append( lowb[i] );
//      ub.append( upb[i] );
//   }
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
	std::vector<double> temp(2);
  for (int i=0; i<16; i++) {
	for(int j=0; j < 2; j++) {
		temp[j] = pts[i][j];
		}
      	data_sp.push_back(DDaceSamplePoint(i, temp));                          
	}
  lb.push_back(0);
  lb.push_back(0);
  ub.push_back(100);
  ub.push_back(100);



}

TestDDaceSampler::~TestDDaceSampler()
{
}

void TestDDaceSampler::run()
{
  testDDaceSampler();
  testDDaceSamplerBase();
  testGetSamples();
  testPrint();
  testTypeName();
  testNSamples();
  testNInputs();
  testGetParameter();
  testDist();
  testLowerBounds();
  testUpperBounds();
  testNoise();
}

void TestDDaceSampler::testDDaceSampler()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;
}

void TestDDaceSampler::testDDaceSamplerBase()
{
  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, true, dists );  
  

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );
}

void TestDDaceSampler::testGetSamples()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  std::vector<DDaceSamplePoint> sp_pass( 0 );
  std::vector<DDaceSamplePoint> sp_ret( 0 );
  

  // call getSamples, should throw an exception
  try {
    //sp_ret = default_data.getSamples( sp_pass );
    default_data.getSamples( sp_pass );
    sp_ret = sp_pass;
    
    // if here then no exception was thrown, record
    // a failure
    _test( false );
  } catch( std::exception& e ) {
    // we're in the catch block so an exception was thrown
    // record a success
    _test( true );
  }

  // call getSamples and compare the returned values
  //sp_ret = base_data.getSamples( sp_pass );
  base_data.getSamples( sp_pass );
  sp_ret = sp_pass;

  _test( Arrcmp_sp( sp_ret, data_sp ) == 0 );
  _test( Arrcmp_sp( sp_pass, data_sp ) == 0 );
}

void TestDDaceSampler::testPrint()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  char buf[1024];
  string data_str;
  string test_str( "<Factorial samples=\"16\" symbols=\"4\" perturb=\"false\" seed=\"0\"/>" );

  ifstream fin;
  ofstream fout;
  fout.open( "TestDDaceSampler_Log" );
  if( !fout )
  {
     cerr << "Couldn't open TestDDaceSampler_Log" << endl;
     return;
  }
  
  try {
    fout << default_data << endl;
    
     // if here then no exception was thrown so record a failure
     _test( false );
  } catch( std::exception& e ) {
    // if here then exceptionw as thrown so record a success
    _test( true );
  }

  // write the data out then close the file
  fout << base_data << endl;
  fout.close();

  // re-open the file
  fin.open( "TestDDaceSampler_Log" );
  if( !fin )
  {
     cerr << "Couldn't open TestDDaceSampler_Log" << endl;
     return;
  }

  // read in the data
  fin.getline( buf,1023 );
  fin.close();
  data_str = buf;

  // check if the data read in is what it is suppose to be
  _test( data_str == test_str );
}

void TestDDaceSampler::testTypeName()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call typeName and test the result
  
  try {  // for this case it should throw an exception
    default_data.typeName();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }
  
  _test( base_data.typeName() == "DDaceFactorialSampler" );
}

void TestDDaceSampler::testNSamples()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call nSamples and test the result

  try {  // for this case it should throw an exception
    default_data.nSamples();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }

  _test( base_data.nSamples() == 16 );
}
  
void TestDDaceSampler::testNInputs()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   
  

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call nInputs and test the result

  try {  // for this case it should throw an exception
    default_data.nInputs();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }

  _test( base_data.nInputs() == 2 );
}

void TestDDaceSampler::testGetParameter()
{
  std::string s;

  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call getParameter, it should throw an exception.
  try {
    default_data.getParameter( s );

    // if we're here then no exception was thrown so record a failure
    _test( false );
  } catch( std::exception& e ) {

    // since we're in the catch block an exception was thrown so 
    // record a success
    _test( true );        
  }

  try {
    base_data.getParameter( s );

    // if we're here then no exception was thrown so record a failure
    _test( false );
  } catch( std::exception& e ) {

    // since we're in the catch block an exception was thrown so
    // record a success
    _test( true );
  }
}

void TestDDaceSampler::testDist()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call dist, it should throw an exception.
  try {
    default_data.dist();

    // if we're here then no exception was thrown so record a failure
    _test( false );
  } catch( std::exception& e ) {

    // since we're in the catch block an exception was thrown so 
    // record a success
    _test( true );        
  }

//  try {
//    base_data.dist();
//
//    // if we're here then no exception was thrown so record a failure
//    _test( false );
//  } catch( ExceptionBase e ) {
//
//    // since we're in the catch block an exception was thrown so
//    // record a success
//    _test( true );
//  }
}

void TestDDaceSampler::testLowerBounds()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

   // call lowerBounds and compare the result

  try {  // for this case it should throw an exception
    default_data.lowerBounds();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }


   _test( Arrcmp_d( base_data.lowerBounds(), lb ) == 0 );
}

void TestDDaceSampler::testUpperBounds()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

   // call lowerBounds and compare the result

  try {  // for this case it should throw an exception
    default_data.upperBounds();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }

   _test( Arrcmp_d( base_data.upperBounds(), ub ) == 0 );
}

void TestDDaceSampler::testNoise()
{
  // create a new DDaceSampler using the default constructor.
  DDaceSampler default_data;

  //DDaceArraySampler ddaceArraySampler( data_d );
  DistributionBase::setSeed(0);
  std::vector<Distribution>       dists;
  dists.resize( 0 );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  dists.push_back( Distribution( UniformDistribution( 0, 100 ) ) );
  DistributionBase::setSeed(0);
  DDaceFactorialSampler sampler( 16, 4, false, dists );   

  // create a new DDaceSampler using the base constructor.
  DDaceSampler base_data( sampler );

  // call noise, it should return false

  try {  // for this case it should throw an exception
    default_data.noise();
 
    // if here then no exception was thrown, record a failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record a success
    _test( true );
  }

  _test( base_data.noise() == false );
}

