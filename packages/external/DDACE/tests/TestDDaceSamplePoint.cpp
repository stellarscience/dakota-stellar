#include "TestDDaceSamplePoint.h"

int Arrcmp( const std::vector<double> & a, const std::vector<double> & b )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b ) {
    return 1;
  }

  for( int i = 0; i < len_a; i++ ) {
    if( a[i] != b[i] ) {
      return 1;
    }
  }

  return 0;
}

int Arrcmp( const std::vector<std::string> & a, 
		const std::vector<std::string> & b )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b ) {
    return 1;
  }

  for( int i = 0; i < len_a; i++ ) {
    if( a[i] != b[i] ) {
      return 1;
    }
  }

  return 0;
}

TestDDaceSamplePoint::TestDDaceSamplePoint()
{
}

TestDDaceSamplePoint::~TestDDaceSamplePoint()
{
}


void TestDDaceSamplePoint::run()
{
  
  testDDaceSamplePoint();       
  testDDaceSamplePointDouble();  
  //testDDaceSamplePointString(); 
  //testGetDataType();              
  //testSetDataType();            
  testIndex();               
  testLength();             
  testOperator();           
  //testGetDoubleValue();     
  //testGetStringValue();     
  testParameters();               
  //testGetDoubleValues();          
  //testGetStringValues();          
  //testPrint();                    
}

void TestDDaceSamplePoint::testDDaceSamplePoint()
{
  // create a sample point using the default constructor
  // test to see if the constructor set the values correctly
  DDaceSamplePoint default_data;
  _test( default_data.index() == 0 );
  _test( default_data.length() == 0 );
  //_test( default_data.getDataType() == DDaceSamplePoint::DOUBLE );

}

void TestDDaceSamplePoint::testDDaceSamplePointDouble()
{
  // create an array of doubles
  std::vector<double> numbers;
  for( int i = 0; i < 10; i++ ) {
    numbers.push_back( i/10.0 );
  }

  // create a sample point using the double constructor
  // test to see if the constructor set the values correctly
  DDaceSamplePoint double_data( 0, numbers );
  _test( double_data.index() == 0 );
  _test( double_data.length() == 10 );
  //_test( double_data.getDataType() == DDaceSamplePoint::DOUBLE );
}

//void TestDDaceSamplePoint::testDDaceSamplePointString()
//{
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  // test to see if the constructor set the values correctly
//  DDaceSamplePoint string_data( 0, strings );
//  _test( string_data.index() == 0 );
//  _test( string_data.length() == 10 );
//  _test( string_data.getDataType() == DDaceSamplePoint::STRING );
//}

//void TestDDaceSamplePoint::testGetDataType()
//{
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // test default constructor
//  _test( default_data.getDataType() == DDaceSamplePoint::DOUBLE );
//
//  // test double constructor
//  _test( double_data.getDataType() ==  DDaceSamplePoint::DOUBLE );
//
//  // test string constructor
//  _test( string_data.getDataType() == DDaceSamplePoint::STRING );
//}

//void TestDDaceSamplePoint::testSetDataType()
//{
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // change the data type and check to see if it worked
//  default_data.setDataType( DDaceSamplePoint::STRING );
//  _test( default_data.getDataType() == DDaceSamplePoint::STRING );
//
//  double_data.setDataType( DDaceSamplePoint::STRING );
//  _test( double_data.getDataType() == DDaceSamplePoint::STRING );
//
//  string_data.setDataType( DDaceSamplePoint::DOUBLE );
//  _test( string_data.getDataType() == DDaceSamplePoint::DOUBLE );
//
//  // restore the values and check to see if it worked
//  default_data.setDataType( DDaceSamplePoint::DOUBLE );
//  _test( default_data.getDataType() == DDaceSamplePoint::DOUBLE );
//
//  double_data.setDataType( DDaceSamplePoint::DOUBLE );
//  _test( double_data.getDataType() == DDaceSamplePoint::DOUBLE );
//
//  string_data.setDataType( DDaceSamplePoint::STRING );
//  _test( string_data.getDataType() == DDaceSamplePoint::STRING );
//}

void TestDDaceSamplePoint::testIndex()
{
  // create a sample point using the default constructor
  DDaceSamplePoint default_data;

  // create an array of doubles
  std::vector<double> numbers;
  for( int i = 0; i < 10; i++ ) {
    numbers.push_back( i/10.0 );
  }

  // create a sample point using the double constructor
  DDaceSamplePoint double_data( 0, numbers );

//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );

  // check the returned index values
  _test( default_data.index() == 0 );
  _test( double_data.index() == 0 );
  //_test( string_data.index() == 0 );
}

void TestDDaceSamplePoint::testLength()
{
  // create a sample point using the default constructor
  DDaceSamplePoint default_data;

  // create an array of doubles
  std::vector<double> numbers;
  for( int i = 0; i < 10; i++ ) {
    numbers.push_back( i/10.0 );
  }

  // create a sample point using the double constructor
  DDaceSamplePoint double_data( 0, numbers );

//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );

  // check the returned length values
  _test( default_data.length() == 0 );
  _test( double_data.length() == 10 );
  //_test( string_data.length() == 10 );
}

void TestDDaceSamplePoint::testOperator()
{
  // create a sample point using the default constructor
  DDaceSamplePoint default_data;

  // create an array of doubles
  std::vector<double> numbers;
  for( int i = 0; i < 10; i++ ) {
    numbers.push_back( i/10.0 );
  }

  // create a sample point using the double constructor
  DDaceSamplePoint double_data( 0, numbers );

//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );

//  // check the [] operator
//  try {
//     _test( default_data[0] == 0 );
//
//  } catch( IndexRangeException e ) {
//  }

  try {
     //_test( double_data[3] == 0.3 );
     _test( fabs(double_data[3]-0.3) < 1.0e-10);
  } catch( std::exception& e ) {
  }

//  try {
//     _test( string_data[3] == 0 ); // returns 0 because dataType=STRING
//  } catch( IndexRangeException e ) {
//  }
}

//void TestDDaceSamplePoint::testGetDoubleValue()
//{
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // check the return value of the function
//  try {
//     _test( default_data.getDoubleValue( 0 ) == 0 );
//  } catch( IndexRangeException e ) {
//  }
//
//  try {
//     _test( double_data.getDoubleValue( 3 ) == 0.3 );
//  } catch( IndexRangeException e ) {
//  }
//
//  try {
//     _test( string_data.getDoubleValue( 3 ) == 0 ); // returns 0 because dataType=STRING
//  } catch( IndexRangeException e ) {
//  }
//}

//void TestDDaceSamplePoint::testGetStringValue()
//{
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // create an empty strin to test against
//  String empty;
//
//  // check the return value of the function
//  try {
//     _test( default_data.getStringValue( 0 ) == empty );  // return "" because dataType=DOUBLE
//  } catch( IndexRangeException e ) {
//  }
//
//  try {
//     _test( double_data.getStringValue( 3 ) == empty );   // same
//  } catch( IndexRangeException e ) {
//  }
//
//  try {
//     _test( string_data.getStringValue( 3 ) == String( "100" ) );
//  } catch( IndexRangeException e ) {
//  }
//}

void TestDDaceSamplePoint::testParameters()
{
  int i;

  // create a sample point using the default constructor
  DDaceSamplePoint default_data;

  // create an array of doubles
  std::vector<double> numbers;
  for( int i = 0; i < 10; i++ ) {
    numbers.push_back( i/10.0 );
  }

  // create a sample point using the double constructor
  DDaceSamplePoint double_data( 0, numbers );

//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // create an empty array to test against
  std::vector<double> arr( 0 );
  _test( Arrcmp( default_data.parameters(), arr ) == 0 );
//
//  // populate an array of 0's to test against
//  for( i = 0; i < 10; i++ ) {
//    arr.append( 0 );
//  } 
//  _test( Arrcmp( string_data.parameters(), arr ) == 0 );
 
  // populate an array of doubles to test against
  arr.resize( 0 );  // clear out the array
  for( i = 0; i < 10; i++ ) {
    arr.push_back( i/10.0 );
  }
  _test( Arrcmp( double_data.parameters(),arr ) == 0 );
}

//void TestDDaceSamplePoint::testGetDoubleValues()
//{
//  int i;
//
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // create an empty array to test against
//  Array<double> arr( 0 );
//  _test( Arrcmp( default_data.getDoubleValues(), arr ) == 0 );
//
//  // populate the array with 0's to test against
//  for( i = 0; i < 10; i++ ) {
//    arr.append( 0 );
//  }
//  _test( Arrcmp( string_data.getDoubleValues(), arr ) == 0 );
//
//  // populate the array with doubles to test against
//  arr.resize( 0 );   // clear the array
//  for( i = 0; i < 10; i++ ) {
//    arr.append( i/10.0 );
//  }
//  _test( Arrcmp( double_data.getDoubleValues(), arr ) == 0 );
//}

//void TestDDaceSamplePoint::testGetStringValues()
//{
//  int i;
//
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  // create an array of empty strings to test against
//  Array<String> arr( 0, "" );
//  _test( Arrcmp( default_data.getStringValues(), arr ) == 0 );
//
//  // resize the array to 10 empty strings
//  arr.resize( 0 );
//  for( i = 0; i < 10; i++ ) {
//    arr.append( String( "" ) );
//  }
//  _test( Arrcmp( double_data.getStringValues(), arr ) == 0 );
//
//  // populate the array with "100"'s
//  arr.resize( 0 );
//  for( i = 0; i < 10; i++ ) {
//    arr.append( String( "100" ) );
//  }
//  _test( Arrcmp( string_data.getStringValues(), arr ) == 0 );
//}

//void TestDDaceSamplePoint::testPrint()
//{
//  int i;
//  char buf[1024];
//
//  // create a sample point using the default constructor
//  DDaceSamplePoint default_data;
//
//  // create an array of doubles
//  Array<double> numbers;
//  for( int i = 0; i < 10; i++ ) {
//    numbers.append( i/10.0 );
//  }
//
//  // create a sample point using the double constructor
//  DDaceSamplePoint double_data( 0, numbers );
//
//  // create an array of strings
//  String str( "100" );
//  Array<String> strings;
//  for( int i = 0; i < 10; i++ ) {
//    strings.append( str );
//  }
//  
//  // create a sample point using the string constructor
//  DDaceSamplePoint string_data( 0, strings );
//
//  string data[6];
//  string test_str[] =  { "[ 0 () ]",
//                         "[ 0 () ]",
//                         "[ 0 (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) ]",
//                         "[ 0 (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) ]",
//                         "[ 0 (100, 100, 100, 100, 100, 100, 100, 100, 100, 100) ]",
//                         "[ 0 (100, 100, 100, 100, 100, 100, 100, 100, 100, 100) ]" };
//
//  ifstream fin;
//  ofstream fout;
//  fout.open( "TestDDaceSamplePoint_Log" );
//
//  if (fout)
//    {
//
//      fout << default_data << endl;
//      default_data.print( fout );
//      fout << endl;
//
//      fout << double_data << endl;
//      double_data.print( fout );
//      fout << endl;
//
//      fout << string_data << endl;
//      string_data.print( fout );
//      fout << endl;
//    }
//
//  fout.close();
//
//  // reopen the file
//  fin.open( "TestDDaceSamplePoint_Log" );
//  
//  if( fin )
//  {
//      for( i = 0; i < 6; i++ )
//      {
//         fin.getline( buf, 1023 );
//         data[i] = buf;
//         _test( data[i] == test_str[i] );
//      }
//  }
//  fin.close();
//}
