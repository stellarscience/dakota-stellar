#include "TestDDaceOALHSampler.h"
//#include "ArgumentMisMatchException.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>

///////////////////////////
// Written By: J Cramp
// Date: July 2004
///////////////////////////

using namespace std;

TestDDaceOALHSampler::TestDDaceOALHSampler() : lb( 1 ), ub( 4 ), seed( 779 )
{
	
   /* For unit tests, we do NOT want a random number generator.              */
   /* We want a number generator that pushes out a known sequence of numbers.*/	
   DistributionBase::usePseudoRandom(true);		
	
  // set the seed to ensure the data matches
  DistributionBase::setSeed( seed );

  // set up the getSamples test data
  int i,j;
  std::vector<double> dtmp( 0 );
  std::vector<int>    itmp( 0 );



  double   pts[4][2] = { {2.335, 2.335},
                         {1.58425,3.8365},
                         {3.08575,1.58425},
                         {3.8365,3.08575} };






  int     ipts[4][2] = { {2,2},
                         {1,4},
                         {3,1},
                         {4,3} };





  int    oapts[4][2] = { {0,0},
                         {0,1},
                         {1,0},
                         {1,1} };
                          
  test_data.resize( 0 );
  test_design.resize( 0 );
  test_oa.resize( 0 );

  for( i = 0; i < 4; i++ )
  {
     test_data.push_back( dtmp );
     test_design.push_back( itmp );
     test_oa.push_back( itmp );
     for( j = 0; j < 2; j++ )
     {
        test_data[i].push_back( pts[i][j] );
        test_design[i].push_back( ipts[i][j] );
        test_oa[i].push_back( oapts[i][j] );
     }
  }
}

TestDDaceOALHSampler::~TestDDaceOALHSampler()
{
   /* Restore the Random Number Generator */	
   DistributionBase::usePseudoRandom(false);		
	
  // reseed the Random Number Generator
  DistributionBase::setSeed( seed );
}

void TestDDaceOALHSampler::run()
{
   testDDaceOALHSampler();
   testGetSamples();
   testClone();
   testPrint();
   testTypeName();
   testGetParameter();
   testGetDesign();
   testGetOA();
}

void TestDDaceOALHSampler::testDDaceOALHSampler()
{  
   // use constructor DDaceOALHSampler( int nSamples, int nInputs,
   //                                   int Strength, int lower, int upper );
   DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );

   // check the lower bounds
   _test( sampler.lowerBound() == lb );
   _test( sampler.upperBound() == ub );
}

void TestDDaceOALHSampler::testGetSamples()
{
   int i,j;

   std::vector<DDaceSamplePoint> pass;
   std::vector<DDaceSamplePoint> ret;
 
   DistributionBase::setSeed( seed );
   DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );
   DistributionBase::setSeed( seed );

   // call getSamples() and test results
   //ret = sampler.getSamples( pass );
   sampler.getSamples( pass );
   DistributionBase::setSeed( seed );
   ret = pass;

   if( _test( Arrcmp_ad_est( pass, test_data, 0.001 ) == 0 ) == false )
     {
       cout << "pass: " << endl;
       for( i = 0; i < (int) pass.size(); i++ )
	 {
           cout << "[";
	   for( j = 0; j < pass[i].length(); j++ )
             {
               cout << "(" << pass[i][j] << "|" << test_data[i][j] << ")";
               if( j != pass[i].length()-1 ) cout << ", ";
             }
           cout << "]" << endl;
         }
     }
  
   _test( Arrcmp_ad_est( ret, test_data, 0.001 ) == 0  );
}

void TestDDaceOALHSampler::testClone()
{
   DDaceOALHSampler  sampler( 4, 2, 2, false, lb, ub );

   // call clone and test the results
   DDaceOALHSampler* rtn = (DDaceOALHSampler*)sampler.clone();

   // check internal values to see if clone() worked
   _test( rtn->typeName() == sampler.typeName() );
   _test( rtn->nSamples() == sampler.nSamples() );
   _test( rtn->nSymbols() == sampler.nSymbols() );
   _test( rtn->nInputs()  == sampler.nInputs() );
   _test( rtn->noise()    == sampler.noise() );
   
   if( _test( (rtn->dist()).size() == sampler.dist().size() ) )
   {
     for( int i = 0; i < (int) (rtn->dist()).size(); i++ )
     {
       _test( (rtn->dist())[i].lowerBound() == sampler.dist()[i].lowerBound() );
       _test( (rtn->dist())[i].upperBound() == sampler.dist()[i].upperBound() );
     }
   }
      
   _test( rtn->lowerBound() == sampler.lowerBound() );
   _test( rtn->upperBound() == sampler.upperBound() );

   // clean up the dynamic memory
   delete rtn;
}

void TestDDaceOALHSampler::testPrint()
{  
  DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );

  char    buf[256];
  char    num[16];
  string  data;
  string  test_str;

  // convert seed from 'long' to 'char*'
  itoa( DistributionBase::seed(), num, 16 );
  
  test_str = "<OrthogonalArrayLatinHypercube samples=\"4\" inputs=\"2\"";
  test_str += " symbols=\"2\" strength=\"2\" frequency=\"1\"";
  test_str += " randomize=\"false\" seed=\"";
  test_str += string( num );
  test_str += "\"/>";
 
  ifstream fin;
  ofstream fout;

  // write data to file then read the data
  // back from the file to verify it out
  // what it was suppose to
  
  fout.open( "TestDDaceOALHSampler_Log" );
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

  fin.open( "TestDDaceOALHSampler_Log" );
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

void TestDDaceOALHSampler::testTypeName()
{
  DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub);

  // call typeName and test the results
  _test( sampler.typeName() == std::string( "DDaceOALHSampler" ) );
}

void TestDDaceOALHSampler::testGetParameter()
{
  DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );

  _test( sampler.getParameter( "SAMPLES" ) == 4 );
  _test( sampler.getParameter( "INPUTS" ) == 2 );
  _test( sampler.getParameter( "SYMBOLS" ) == 2 );
  _test( sampler.getParameter( "STRENGTH" ) == 2 );
  _test( sampler.getParameter( "FREQUENCY" ) == 1 );
  _test( sampler.getParameter( "RANDOMIZED" ) == 0 );

  try {
    sampler.getParameter( "WRONG_INPUT" );
  
    // if here then no exception was thrown, record failure
    _test( false );
  } catch( std::exception& e ) {
    // if here then exception was thrown, record success
    _test( true );
  }
}

void TestDDaceOALHSampler::testGetDesign()
{
   int i,j;
   DistributionBase::setSeed( seed );
   DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );
   DistributionBase::setSeed( seed );
   std::vector<std::vector<int> > ret = sampler.getDesign();
   DistributionBase::setSeed( seed );

   if( _test( Arrcmp_i( ret, test_design ) == 0 ) == false )
     {
       std::vector<std::vector<int> > oa = sampler.getOA();
       
       cout << "Design: " << endl;
       for( i = 0; i < (int) ret.size(); i++ )
	 {
           cout << "[";
	   for( j = 0; j < (int) ret[i].size(); j++ )
             {
               cout << "(" << ret[i][j] << "|" << test_design[i][j] << ")";
               if( j != (int) ret[i].size()-1 ) cout << ", ";
             }
           cout << "]" << endl;
         }
       cout << "OA: " << endl;
       for( i = 0; i < (int) oa.size(); i++ )
	 {
           cout << "[";
	   for( j = 0; j < (int) oa[i].size(); j++ )
             {
               cout << oa[i][j];
               if( j != (int) oa[i].size()-1 ) cout << ", ";
             }
           cout << "]" << endl;
         }
     }
}

void TestDDaceOALHSampler::testGetOA()
{
  int i,j;
 
  // due to the some randomness in the
  // OA generation this particular sampler
  // object uses a different OA than the
  // others created in the tests
  int pts[4][2] = { {0,0},
                    {0,1},
                    {1,0},
                    {1,1} };

  DDaceOALHSampler sampler( 4, 2, 2, false, lb, ub );

  std::vector<std::vector<int> > oa( 4 );
  for( i = 0; i < 4; i++ )
  {
    oa[i].resize( 2 );
    for( j = 0; j < 2; j++ )
    {
      oa[i][j] = pts[i][j];
    }
  }

   std::vector<std::vector<int> > ret = sampler.getOA();
   if( _test( Arrcmp_i( ret, oa ) == 0 ) == false )
     {
       cout << "OA: " << endl;
       for( i = 0; i < (int) ret.size(); i++ )
	 {
           cout << "[";
	   for( j = 0; j < (int) ret[i].size(); j++ )
             {
               cout << "(" << ret[i][j] << "|" << oa[i][j] << ")";
               if( j != (int) ret[i].size()-1 ) cout << ", ";
             }
           cout << "]" << endl;
         }
     }
} 
 

