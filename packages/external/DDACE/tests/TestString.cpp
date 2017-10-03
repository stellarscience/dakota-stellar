#include "TestString.h"

using namespace std;

TestString::TestString()
{
}

TestString::~TestString()
{
}

void TestString::run()
{
  testStringConstructors();
  testNewStringConstructors();
  testAccessors();
  testOperators();
  testSearch();
  testMiscFuncs();
}

void TestString::testStringConstructors()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the String constructor
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Default Constructor
  String defaultStr;
  _test( !strcmp(defaultStr.c_str(), "" ));

  // Single char constructor
  const char* cstr = "c";
  String charStr( 'c' );
  _test( !strcmp(charStr.c_str(), cstr));

  // C string constructor
  String cstrStr( "Hello World" );
  const char* c2str = "Hello World";
  _test( !strcmp(cstrStr.c_str(), c2str) );

  // Selected Element Constructor
  String someElems( "Hello World", 5);
  _test( !strcmp(someElems.c_str(), "Hello"));
  String noElems( "Hello World" , 0);
  _test( !strcmp(noElems.c_str(), "") ); 

  // Copy Constructor
  String copyStr( "Hello World");
  String other(copyStr);
  _test( !strcmp(other.c_str(), copyStr.c_str()));
  // 5 tests in this function
}

void TestString::testNewStringConstructors()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using the new STL String Constructors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Default constructor
  String defaultStr;
//  string tmp = defaultStr.stlString();
  _test( !defaultStr.compare("") );

  // Single char constructor
  String charStr('c');
//  tmp = charStr.stlString();
  _test( !charStr.compare("c"));

  String nullcharStr('\0');
//  tmp = nullcharStr.stlString();
  _test( !nullcharStr.compare(""));

  // C string constructor
  String cstrStr( "Hello World" );
  _test( !cstrStr.compare("Hello World"));

  String nullcStr( "" );
  _test( !nullcStr.compare(""));

  // Selected Element constructor
  String someElems( "Hello World", 5);
  _test( !someElems.compare("Hello"));
  String noElems( "Hello World" , 0);
  _test( !noElems.compare("") ); 

  // Copy constructor
  String copyStr( "Hello World");
  String other(copyStr);
  _test( !copyStr.compare(other.c_str()));
  // note that the null string copy constructor is not possible, at least
  // using gcc
}

void TestString::testAccessors()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests using old Accessors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~
  String accessorStr("Hello World");
  String nullStr;
  const char* tmp2 = accessorStr.c_str();

  //Test length
  _test( accessorStr.length() == 11); 
  _test( nullStr.length() == 0);


  // Read-only Element Access
  _test( accessorStr[1] == 'e');
}

void TestString::testOperators()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Tests for Operators. Note that there can only be one set of
  // operators at a time in a class, therefore there will not be
  // a testNewOperators() method.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  String operatorStr("Hello World");
  String nullStr;
  String otherStr("Hello World");
  String other2Str("Hi Planet");

  _test( operatorStr == "Hello World" );
  _test( nullStr == "" );
  _test( operatorStr == otherStr );
  _test( operatorStr != other2Str );
  _test( operatorStr != nullStr );

//  _test( operatorStr < other2Str );

  operatorStr += other2Str;
  _test( operatorStr == "Hello WorldHi Planet");
  operatorStr += nullStr;
  _test( operatorStr == "Hello WorldHi Planet");

  String yetAnotherStr = otherStr + other2Str;
  _test( yetAnotherStr == "Hello WorldHi Planet");
  yetAnotherStr = otherStr + nullStr;
  _test( yetAnotherStr == "Hello World" );

  yetAnotherStr = otherStr + 2;
  _test( yetAnotherStr == "Hello World2" );
  yetAnotherStr = nullStr + 2;
  _test( yetAnotherStr == "2" );

  yetAnotherStr = otherStr + 3.432;
  _test( yetAnotherStr == "Hello World  3.4319999999999999e+00" );
  yetAnotherStr = nullStr + 3.432;
  _test( yetAnotherStr == "  3.4319999999999999e+00" );


 // 14 tests in this function, 29 tests up to this point
}

void TestString::testSearch()
{
  String search1 = "Hello World";
  String search2 = "lo";
  String search3 = "qzt";
  String search4 = "Hel";
  String search5 = " World";
  String search6, search7;
  String nullStr;

// find method. taking both a c string and a String as input
  _test( search1.find("lo") == 3 );
  _test( search1.find(search2) == 3 );
  _test( search1.find("qzt") == -1 );
  _test( search1.find(search3) == -1);

  _test( nullStr.find("lo") == -1 );
  _test( nullStr.find(search2) == -1 );
/**
// before method. again, taking both a c string and a String as input
  _test( search1.before("lo") == search4 );
  _test( search1.before(search2) == search4 );
  _test( search1.before("qzt") == search1 );
  _test( search1.before(search3) == search1 );
  _test( search1.before(nullStr) == nullStr );

  _test( nullStr.before("lo") == nullStr );
  _test( nullStr.before(search2) == nullStr );
  _test( nullStr.before(nullStr) == nullStr );

// after method. taking both a c strign and a String as input
  _test( search1.after("lo") == search5 );
  _test( search1.after(search2) == search5 );
  _test( search1.after("qzt") == nullStr );
  _test( search1.after(search3) == nullStr );
  _test( search1.after(nullStr) == nullStr );

  _test( search1.after("rld") == nullStr);
  _test( search1.after("d") == nullStr);

  _test( nullStr.after("qzt") == nullStr );
  _test( nullStr.after(search3) == nullStr );
  _test( nullStr.after(nullStr) == nullStr );
*/
/**
// between method. only one version
  _test( search1.between("H","r",search6, search7) == "ello Wo");
  _test( search6 == "" &&  search7 == "ld" );
  cout << "s1: " << search1.between("H","r",search6, search7) << "***" << endl;
  cout << "s6: " << search6 << "***" << "s7: " << search7 << "***" << endl;
  _test( search1.between("llo","or", search6, search7) == " W" );
  _test( search6 == "He" && search7 == "ld");
  _test( search1.between("qzt","rld",search6, search7) == "");
  _test( search6 == "Hello Wo" &&  search7 == "");
  _test( search1.between("H","qzt", search6, search7) == "");
  _test( search6 == "ello World" && search7 == "");
  _test( search1.between("ell","qzt", search6, search7) == "");
  _test( search6 == "H" && search7 == "o World");
  _test( search1.between(nullStr,"Wor", search6, search7) == "");
  _test( search6 == "Hello " && search7 == "ld" );
  _test( search1.between("el", nullStr, search6, search7) == "");
  _test( search6 == "H" && search7 == "lo World");
  _test( search1.between(nullStr, nullStr, search6, search7) == "");
  _test( search6 == "Hello World" && search7 == "");

  _test( nullStr.between("H","qzt",search6, search7) == "");
  _test( search6 == "" && search7 == "");
  _test( nullStr.between(nullStr,"qzt", search6, search7) == "");
  _test( search6 == "" && search7 == "");
  _test( nullStr.between("H",nullStr, search6, search7) == "" );
  _test( search6 == "" && search7 == "");
  _test( nullStr.between(nullStr,nullStr, search6, search7) == "");
  _test( search6 == "" && search7 == "");


// subString method
  _test( search1.subString(0,0) == nullStr );
  _test( search1.subString(0,5) == "Hello" );
  _test( search1.subString(11,11) == nullStr );
  _test( search1.subString(6,11) == "World" );
  
  String search8("12345");
  _test( search8.subString(1,2) == "2");

  _test( nullStr.subString(0,0) == "");
*/
}

void TestString::testMiscFuncs()
{
  String misc1( "Hello World" );
  String misc2( " \n\t\r" );
  String misc3( "    Hi Planet" );
  String misc4( "GOOD DAY");
  String nullStr;
/**
  _test( misc1.allCaps() == "HELLO WORLD" );
  _test( misc2.allCaps() == misc2 );
  _test( misc4.allCaps() == misc4 );
  _test( nullStr.allCaps() == nullStr );

// trimInitialWhitespace.
  _test( misc1.trimInitialWhitespace() == misc1 );
  _test( misc2.trimInitialWhitespace() == nullStr );
  _test( misc3.trimInitialWhitespace() == "Hi Planet" );
  _test( nullStr.trimInitialWhitespace() == nullStr );
*/
// conversion tools
/**
  String num1("1");
  String num2("3.029384203");
  String bool1("True");
  String bool2("TrUE");
 
  _test( num1.atoi() == 1 );
  _test( num2.atof() - 3.029384203 <= .000000001 || num2.atof() - 3.029384203 >= .000000001 );
  _test( bool1.atob() );
  _test( bool2.atob() );
  _test( !misc1.atob() );
*/
}
