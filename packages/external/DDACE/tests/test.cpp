// test.cpp

/*
 * C/C++ Users Journal Sept 2000 <br>
 * The Simplest Automated Unit Test Framework That Could Possibly Work <br>
 * Chuck Allison <br>
 */

#include "test.h"
#include <iostream>
#include <string>
#include <typeinfo>     // Visual Studio requires /GR""

#ifdef _MSC_VER
//Allow return-less mains:
#pragma warning(disable: 4541)
#endif

using namespace std;

bool Test::do_test(bool cond, const std::string& lbl,
                   const char* fname, long lineno)
{
    if (!cond)
    {
        do_fail(lbl, fname, lineno);
        return false;
    }
    else
        _succeed();

    return true;
}

void Test::do_fail(const std::string& lbl,
                   const char* fname, long lineno)
{
    ++m_nFail;
    if (m_osptr)
    {
        *m_osptr << typeid(*this).name()
                             << "failure: (" << lbl << ") , "
                                 << fname
                 << " (line " << lineno << ")\n";
    }
}

long Test::report() const
{
    // Use RTTI to get class name then
    // parse out the junk at the start
    // of the string.
    // Requires that all class names 
    // begin with "Test".
    string name( typeid(*this).name() );
    name.erase( 0, name.find( "Test" ) );

    if (m_osptr)
        {
            *m_osptr << "Test \""  << name << "\":" << endl 
                     //<< this->getClassName().c_str() << "\":" << endl
                     << "\tPassed: " << m_nPass
                     << "\tFailed: " << m_nFail
                     << endl;
        }
    return m_nFail;
}


