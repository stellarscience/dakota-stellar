#include "TestOneWayANOVA.h"
#include "arrcmp.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;


TestOneWayANOVA::TestOneWayANOVA()
{
}

TestOneWayANOVA::~TestOneWayANOVA()
{
}

void TestOneWayANOVA::run()
{
   testConstructor();
   testGetANOVATables();
}

void TestOneWayANOVA::testConstructor(){

    std::string reference = "";
    reference.append("ANOVA Table for Factor 0\n");
    reference.append("Source of              Sum of    Mean Sum of     \n");
    reference.append("Variation        DoF   Squares   Squares        Fdata       p\n");
    reference.append("Between Groups    1    139.599   139.599      40.1621    0.0000\n");
    reference.append("Within Groups     8     27.807     3.47587\n");
    reference.append("Total             9    167.406\n");

    std::vector<int> dataInput(10);
    std::vector<double> dataOutput(10);
    dataInput[0] = 0;  dataOutput[0]=11.06;
    dataInput[1] = 0;  dataOutput[1]=15.4;
    dataInput[2] = 0;  dataOutput[2]=11.81;
    dataInput[3] = 1;  dataOutput[3]=20.48;
    dataInput[4] = 1;  dataOutput[4]=20.32;    
    dataInput[5] = 0;  dataOutput[5]=10.96;
    dataInput[6] = 0;  dataOutput[6]=15.86;
    dataInput[7] = 0;  dataOutput[7]=10.8;
    dataInput[8] = 1;  dataOutput[8]=19.73;
    dataInput[9] = 1;  dataOutput[9]=20.57;
    
    DDaceMainEffects::Response response(dataOutput);
    DDaceMainEffects::Factor factor = 
         DDaceMainEffects::Factor(dataInput, 2, response);
    
    std::vector<DDaceMainEffects::Factor> factors;
    factors.push_back(factor);
    
    DDaceMainEffects::OneWayANOVA oneWayANOVA(factors);	
    
    oneWayANOVA.printANOVATables();
    std::string output = oneWayANOVA.getANOVATables();
    _test(output==reference);
	
}

void TestOneWayANOVA::testGetANOVATables()
{
	
    std::string reference = "";
    reference.append("Source of              Sum of      Mean Sum of\n");
    reference.append("Variation        DoF   Squares     Squares     Fdata      p \n");
    reference.append("Between Groups     1    139.599   139.599      40.1621    0.0000\n");
    reference.append("Within Groups      8     27.807     3.47587\n");
    reference.append("Total              9    167.406\n");

    std::vector<int> dataInput(10);
    std::vector<double> dataOutput(10);
    dataInput[0] = 0;  dataOutput[0]=11.06;
    dataInput[1] = 0;  dataOutput[1]=15.4;
    dataInput[2] = 0;  dataOutput[2]=11.81;
    dataInput[3] = 1;  dataOutput[3]=20.48;
    dataInput[4] = 1;  dataOutput[4]=20.32;    
    dataInput[5] = 0;  dataOutput[5]=10.96;
    dataInput[6] = 0;  dataOutput[6]=15.86;
    dataInput[7] = 0;  dataOutput[7]=10.8;
    dataInput[8] = 1;  dataOutput[8]=19.73;
    dataInput[9] = 1;  dataOutput[9]=20.57;
    
    DDaceMainEffects::Response response(dataOutput);
    DDaceMainEffects::Factor factor = DDaceMainEffects::Factor(dataInput, 2, response);
    
    std::vector<DDaceMainEffects::Factor> factors;
    factors.push_back(factor);
    
    DDaceMainEffects::OneWayANOVA oneWayANOVA(factors);	
    
    std::string output = oneWayANOVA.getANOVATables();
    _test(output==reference);
    
}


