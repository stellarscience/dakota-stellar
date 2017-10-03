#include "TestMarsAnalyzer.h"
#include <cstdlib>
#include "DDaceSampler.h"
#include "DDaceUserInputSampler.h"
#include "FuncApprox.h"
#include "Mars.h"
#include "NormalDistribution.h"
#include "UniformDistribution.h"
#include "Distribution.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void myFunc(const std::vector<double>& in, std::vector<double>& out);

TestMarsAnalyzer::TestMarsAnalyzer() : seed( 889 )
{

	std::vector< std::vector < double > > funcOut;
	ofstream outfile("marsOutput.txt");
	int nSamples = 100;

  	// set random number generator seed
  	DistributionBase::setSeed( seed );

 	 // create distributions need by sampler
  	dists.resize( 0 );
  	dists.push_back( Distribution( UniformDistribution( -1.0, 1.0 ) ) );
  	dists.push_back( Distribution( NormalDistribution( Mean(0.0),
						  StdDeviation(1.0),
						  3.0)));

	// Resize, seed, and initialize.
	DistributionBase::setSeed( seed );
	DDaceSampler sampler = DDaceUserInputSampler("samplerData.dat");
	DDaceSamplePoint sample;

	std::vector<std::string> varNames(2);
	varNames[0] = "x1";
	varNames[1] = "x2";

	std::vector<std::string> outputNames(1, "y(x1,x2)");

	// Use the sampler to populate data
  	DistributionBase::setSeed( seed );
	sampler.getSamples( data );

	int len = data.size();
	funcOut.resize(len);
	for(int i = 0; i < len; i++)
	{
	  sample = data[i];
	  std::vector<double> value(1);
	  myFunc(sample.parameters(), value);
	  funcOut[i] = value;
	}

  	DistributionBase::setSeed( seed );

	//cout << "x = " << data << endl;

	//cout << "f = " << funcOut << endl;

	// y = sliceArray(funcOut, 0);

	// this loop replaces the sliceArray function from the Array class
	y = std::vector<double>(funcOut.size());
	for(int i = 0; i < (int) funcOut.size(); i++) y[i] = funcOut[i][0];

	//cout << "y = " << y << endl;

	FuncApprox fa = Mars(2, nSamples);

	fa.setBounds(sampler.lowerBounds(), sampler.upperBounds());

	fa.setNPtsPerDim(11); // compute 11 MARS values for each var

	settings = std::vector<double>(2, 0.0);

	try
	{
  		DistributionBase::setSeed( seed );
		fa.write2DGridData("funcApprox.dat", data, y, 0, 1, settings);
	}
	catch(std::exception& e)
	{ cerr << "in write2DGrid " << endl; }


}

TestMarsAnalyzer::~TestMarsAnalyzer()
{
}

void TestMarsAnalyzer::run()
{
   testMarsOutput();
}

bool TestMarsAnalyzer::close(double a, double b)
{
   return (a-b < .001 || b-a < 0.001);
}

void TestMarsAnalyzer::testMarsOutput()
{
	fstream testData("funcApprox.dat");
	fstream correctData("funcApproxCorrect.dat");
	string test, correct;
	while(!testData.eof() || !correctData.eof())
	{
	  testData >> test ;
	  correctData >> correct;
	  _test( close(atof(test.c_str()), atof(correct.c_str()))); 
	}
	_test( testData.eof() && correctData.eof() );
	correctData.close();
	testData.close();

}


void myFunc(const std::vector<double>& in, std::vector<double>& out)
    {
//    int wait = (int) floor(1 + 3*DistributionBase::uniformUnitDeviate());
    //sleep(wait);
    out.resize(1);
    out[0] = in[0]*in[0] + 9.0*in[1]*in[1];
}
