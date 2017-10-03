#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include "TestMainEffectsExcelOutput.h"
#include "TestMainEffectsConverter.h"
#ifdef HAVE_STANDALONE
#include "TestOneWayANOVA.h"
#endif // HAVE_STANDALONE
#include "TestFactor.h"
#include "TestResponse.h"
#include "TestDDaceArraySampler.h"
//#include "TestDDaceStringArraySampler.h"
#include "TestMean.h"
#include "TestStdDeviation.h"
#include "TestDDaceSamplePoint.h"
#include "TestDDaceSampler.h"
#include "TestDDaceBoxBehnkenSampler.h"
#include "TestDDaceCentralCompositeSampler.h"
#include "TestDDaceFactorialSampler.h"
#include "TestDDaceLHSampler.h"
#include "TestDDaceOASampler.h"
#include "TestDDaceRandomSampler.h"
#include "TestDDaceUserInputSampler.h"
#include "TestUniformDistribution.h"
#include "TestNormalDistribution.h"
#include "TestDistribution.h"
#include "TestDDaceOALHSampler.h"
#ifdef HAVE_STANDALONE
#include "TestMarsAnalyzer.h"
#endif // HAVE_STANDALONE
#include "TestMainEffectsAnalyzer.h"
#include "TestPseudoRandom.h"
#ifdef _MSC_VER
#include <direct.h>
#define getcwd _getcwd
#define chdir _chdir
#else
#include <unistd.h>
#endif

using namespace std;

int main(int argc, char **argv) {

try{
	
    /* When we are using the command line,                                 */
    /* the Makefile cd's to this directory.                                */
    /* When we are using an IDE,                                           */
    /* the current working directory is wherever we launched the IDE from. */
    /* If the current working directory does not end in "tests"            */
    /* then cd to tests.                                                   */
   char nameOfCurrentWorkingDirectory[1024];
   char nameOfThisFolder[] = "tests";
   getcwd(nameOfCurrentWorkingDirectory, 1024);
   char *ptr = strstr(nameOfCurrentWorkingDirectory, nameOfThisFolder);
   if (ptr==NULL) {
   	    chdir(nameOfThisFolder);
   }
   
//    PMachine::init(argc, (void**)argv);


    Suite s("DDace Test Suite", &cout);
    s.addTest(new TestMainEffectsExcelOutput);
    s.addTest(new TestMainEffectsConverter);
#ifdef HAVE_STANDALONE
    s.addTest(new TestOneWayANOVA);
#endif // HAVE_STANDALONE
    s.addTest(new TestFactor);
    s.addTest(new TestResponse);
    s.addTest(new TestPseudoRandom);
    s.addTest(new TestDDaceArraySampler);
    //s.addTest(new TestDDaceStringArraySampler);
    s.addTest(new TestDDaceSamplePoint);
    s.addTest(new TestDDaceSampler);
    s.addTest(new TestDDaceBoxBehnkenSampler);
    s.addTest(new TestDDaceCentralCompositeSampler);
    s.addTest(new TestDDaceFactorialSampler);
    s.addTest(new TestDDaceLHSampler);
    s.addTest(new TestDDaceOASampler);
    s.addTest(new TestDDaceRandomSampler);
    s.addTest(new TestDDaceUserInputSampler);
    s.addTest(new TestUniformDistribution);
    s.addTest(new TestNormalDistribution);
    s.addTest(new TestDDaceOALHSampler);
#ifdef HAVE_STANDALONE
    s.addTest(new TestMarsAnalyzer);
#endif // HAVE_STANDALONE
    s.addTest(new TestMainEffectsAnalyzer);
    // if you are adding a test which involves
    // random numbers in anyway put it after this
    // line or else you are going to break
    // some of the previous tests
    s.addTest(new TestDistribution);
    s.addTest(new TestMean);
    s.addTest(new TestStdDeviation);
    s.run();
 

    /* if I do NOT have any command line arguments */
    /* then print the entire report.               */
    long nFail = 0;
    if (argc==1) { 
        nFail = s.report();
        s.free();
        cout << "\nTotal Failures: " << nFail << endl;

    /* if I have one or more command line arguments */
    /*     if there are no failures then print 0    */
    /*     else print the entire report             */
    } else {
        nFail = s.getNumFailed();
        cout << nFail << endl;
        if (nFail != 0) {
            s.report();
            cout << "\nTotal Failures: " << nFail << endl;
        }
	s.free();
    }
    
/**
    PMachine::synchronize(PWorld());

    PMachine::finalize();
**/
return(nFail);

}
catch(std::exception& e)
{
	cerr << "The following exception occured: " << e.what() << endl;
}
/**catch(...)
{
	cerr << "An exception occurred, but I don't know what it was, sorry." << endl;
}
*/
}

