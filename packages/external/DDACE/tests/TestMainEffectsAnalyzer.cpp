#include "TestMainEffectsAnalyzer.h"
#include "DDaceSamplePoint.h"
#include "Distribution.h"
#include "arrcmp.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;



TestMainEffectsAnalyzer::TestMainEffectsAnalyzer() 
    : data(10,vector<DataValue>(5)), columnHeaders(5), seed( 889 )
{
        columnHeaders[0] = 
            ColumnHeader("Screw Length", "in", ColumnHeader::FACTOR);
        columnHeaders[1] = 
            ColumnHeader("Speed", "in/sec", ColumnHeader::FACTOR);
        columnHeaders[2] = 
            ColumnHeader("L/D", "", ColumnHeader::FACTOR);
        columnHeaders[3] = 
            ColumnHeader("Elastic Modulus", 
                         "x10^6 psi", 
                         ColumnHeader::RESPONSE);
        columnHeaders[4] = 
            ColumnHeader("Yield Stress", 
                         "x10^3 psi", 
                         ColumnHeader::RESPONSE);
        

       /* trial 1, run 1 */
       data[0][0] = DataValue(1);      // Screw Length
       data[0][1] = DataValue(1);      // Speed
       data[0][2] = DataValue(1);      // D/L
       data[0][3] = DataValue(11.06);  //Elastic Modulus
       data[0][4] = DataValue(158.9);  //Yield Stress

       /* trial 1, run 2 */
       data[1][0] = DataValue(1);      // Screw Length
       data[1][1] = DataValue(1);      // Speed
       data[1][2] = DataValue(1);      // D/L
       data[1][3] = DataValue(10.96);  //Elastic Modulus
       data[1][4] = DataValue(157.02);  //Yield Stress

       /* trial 2, run 1 */
       data[2][0] = DataValue(1);      // Screw Length
       data[2][1] = DataValue(2);      // Speed
       data[2][2] = DataValue(1);      // D/L
       data[2][3] = DataValue(15.4);   //Elastic Modulus
       data[2][4] = DataValue(191.58); //Yield Stress

       /* trial 2, run 2 */
       data[3][0] = DataValue(1);      // Screw Length
       data[3][1] = DataValue(1);      // Speed
       data[3][2] = DataValue(1);      // D/L
       data[3][3] = DataValue(15.86);  //Elastic Modulus
       data[3][4] = DataValue(191.24); //Yield Stress

       /* trial 3, run 1 */
       data[4][0] = DataValue(1);      // Screw Length
       data[4][1] = DataValue(2);      // Speed
       data[4][2] = DataValue(1);      // D/L
       data[4][3] = DataValue(11.81);  //Elastic Modulus
       data[4][4] = DataValue(158.7);  //Yield Stress

       /* trial 3, run 2 */
       data[5][0] = DataValue(1);      // Screw Length
       data[5][1] = DataValue(2);      // Speed
       data[5][2] = DataValue(1);      // D/L
       data[5][3] = DataValue(10.8);   //Elastic Modulus
       data[5][4] = DataValue(152.09); //Yield Stress

       /* trial 4, run 1 */
       data[6][0] = DataValue(2);      // Screw Length
       data[6][1] = DataValue(1);      // Speed
       data[6][2] = DataValue(2);      // D/L
       data[6][3] = DataValue(20.48);  //Elastic Modulus
       data[6][4] = DataValue(166.48); //Yield Stress

       /* trial 4, run 2 */
       data[7][0] = DataValue(2);      // Screw Length
       data[7][1] = DataValue(1);      // Speed
       data[7][2] = DataValue(2);      // D/L
       data[7][3] = DataValue(19.73);  //Elastic Modulus
       data[7][4] = DataValue(165.5);  //Yield Stress

       /* trial 5, run 1 */
       data[8][0] = DataValue(2);      // Screw Length
       data[8][1] = DataValue(1);      // Speed
       data[8][2] = DataValue(2);      // D/L
       data[8][3] = DataValue(20.32);  //Elastic Modulus
       data[8][4] = DataValue(165.28); //Yield Stress

       /* trial 5, run 2 */
       data[9][0] = DataValue(2);      // Screw Length
       data[9][1] = DataValue(1);      // Speed
       data[9][2] = DataValue(2);      // D/L
       data[9][3] = DataValue(20.57);  //Elastic Modulus
       data[9][4] = DataValue(167.43); //Yield Stress


   x = MainEffectsAnalyzer3(columnHeaders, data);
}

TestMainEffectsAnalyzer::~TestMainEffectsAnalyzer()
{
}

void TestMainEffectsAnalyzer::run()
{
   testAllMainEffects();
}

void TestMainEffectsAnalyzer::testAllMainEffects()
{
   _test(x.toIndexInputColumn("Screw Length") == 0);
   _test(x.toIndexInputColumn("SL") == 0);
   _test(x.toIndexInputColumn("a") == 0);
   _test(x.toIndexInputColumn("A") == 0);

   _test(x.toIndexInputColumn("Speed") == 1);
   _test(x.toIndexInputColumn("S") == 1);
   _test(x.toIndexInputColumn("b") == 1);
   _test(x.toIndexInputColumn("B") == 1);

   _test(x.toIndexInputColumn("L/D") == 2);
   _test(x.toIndexInputColumn("LD") == 2);
   _test(x.toIndexInputColumn("c") == 2);
   _test(x.toIndexInputColumn("C") == 2);

   _test(x.toIndexInputColumn("Elastic Modulus") == 3);
   _test(x.toIndexInputColumn("EM") == 3);
   _test(x.toIndexInputColumn("d") == 3);
   _test(x.toIndexInputColumn("D") == 3);

   _test(x.toIndexInputColumn("Yield Stress") == 4);
   _test(x.toIndexInputColumn("YS") == 4);
   _test(x.toIndexInputColumn("e") == 4);
   _test(x.toIndexInputColumn("E") == 4);  

   _test(x.getNumberOfObservations(0,3) == 10);
   _test(x.getNumberOfObservations("Screw Length",3) == 10);
   _test(x.getNumberOfObservations("SL",3) == 10);
   _test(x.getNumberOfObservations("a",3) == 10);
   _test(x.getNumberOfObservations("A",3) == 10);
   _test(x.getNumberOfObservations(0,"Elastic Modulus") == 10);
   _test(x.getNumberOfObservations(0,"EM") == 10);
   _test(x.getNumberOfObservations(0,"d") == 10);
   _test(x.getNumberOfObservations(0,"D") == 10);
   _test(x.getNumberOfObservations(0,4) == 10);

   DataValue factorValue = DataValue(1);
   _test(x.getNumberOfObservations(0,factorValue,3) == 6);
   _test(x.getNumberOfObservations("Screw Length",factorValue,3) == 6);
   _test(x.getNumberOfObservations("SL",factorValue,3) == 6);
   _test(x.getNumberOfObservations("a",factorValue,3) == 6);
   _test(x.getNumberOfObservations("A",factorValue,3) == 6);
   _test(x.getNumberOfObservations(0,factorValue,"Elastic Modulus") == 6);
   _test(x.getNumberOfObservations(0,factorValue,"EM") == 6);
   _test(x.getNumberOfObservations(0,factorValue,"d") == 6);
   _test(x.getNumberOfObservations(0,factorValue,"D") == 6);
   _test(x.getNumberOfObservations(0,1,3) == 6);
   _test(x.getNumberOfObservations("Screw Length",1,3) == 6);
   _test(x.getNumberOfObservations("SL",1,3) == 6);
   _test(x.getNumberOfObservations("a",1,3) == 6);
   _test(x.getNumberOfObservations("A",1,3) == 6);
   _test(x.getNumberOfObservations(0,1,"Elastic Modulus") == 6);
   _test(x.getNumberOfObservations(0,1,"EM") == 6);
   _test(x.getNumberOfObservations(0,1,"d") == 6);
   _test(x.getNumberOfObservations(0,1,"D") == 6);
   _test(x.getNumberOfObservations(0,factorValue,4) == 6);
   

    factorValue = DataValue(2);
   _test(x.getNumberOfObservations(0,factorValue,3) == 4);
   _test(x.getNumberOfObservations(0,factorValue,4) == 4);

   _test(closeEnough(x.getSumOfObservations(0,3), 156.99, 0.0001));
   _test(closeEnough(x.getSumOfObservations(0,4), 1674.22, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getSumOfObservations(0,factorValue,3),75.89, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfObservations(0,factorValue,3), 81.10, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getSumOfObservations(0,factorValue,4), 1009.53, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfObservations(0,factorValue,4), 664.69, 0.0001));

   _test(closeEnough(x.getAverageObservation(0,3), 15.699, 0.0001));
   _test(closeEnough(x.getAverageObservation(0,4), 167.422, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getAverageObservation(0,factorValue,3), 12.648333, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getAverageObservation(0,factorValue,3), 20.275, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getAverageObservation(0,factorValue,4), 168.255, 0.0001));

    factorValue = DataValue(2);
   _test(closeEnough(x.getAverageObservation(0,factorValue,4), 166.1725, 0.0001));  

   _test(closeEnough(x.getSumOfSquares(0,3), 167.405490, 0.0001));
   _test(closeEnough(x.getSumOfSquares(0,4), 1652.045360, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getSumOfSquares(0,factorValue,3), 27.378883, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfSquares(0,factorValue,3), 0.428100, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getSumOfSquares(0,factorValue,4), 1638.712350, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfSquares(0,factorValue,4), 2.924675, 0.0001));

   _test(closeEnough(x.getVariance(0,3), 18.600610, 0.0001));
   _test(closeEnough(x.getVariance(0,4), 183.560596, 0.0001));

   factorValue = DataValue(1);
   _test(closeEnough(x.getVariance(0,factorValue,3), 5.475777, 0.0001));

   factorValue = DataValue(2);
   _test(closeEnough(x.getVariance(0,factorValue,3), 0.142700, 0.0001));
   factorValue = DataValue(1);
   _test(closeEnough(x.getVariance(0,factorValue,4), 327.742470, 0.0001));
   factorValue = DataValue(2);
   _test(closeEnough(x.getVariance(0,factorValue,4), 0.974892, 0.0001)); 

   _test(x.getD(0,3) == 9);
   _test(x.getD(0,4) == 9);

   factorValue = DataValue(1);
   _test(x.getD(0,factorValue,3) == 5);
   factorValue = DataValue(2);
   _test(x.getD(0,factorValue,3) == 3);
   factorValue = DataValue(1);
   _test(x.getD(0,factorValue,4) == 5);
   factorValue = DataValue(2);
   _test(x.getD(0,factorValue,4) == 3); 

   std::vector< DataValue > uniqueFactors = 
            x.getNonEmptyUniqueFactors(1);
   std::vector< DataValue >::iterator iterator(uniqueFactors.begin());
   string tmp(((DataValue)(*iterator)).toString());
   string tmp2("DataValue:dataType=integer value=1");
   _test(!tmp.compare(tmp2));
   iterator++;
   tmp = ((DataValue)(*iterator)).toString();
   tmp2 = "DataValue:dataType=integer value=2";
   _test(!tmp.compare(tmp2));

   _test(closeEnough(x.getSumOfSquaresBetweenGroups(0,3), 139.598507, 0.0001));
   _test(closeEnough(x.getSumOfSquaresBetweenGroups(0,4), 10.408335, 0.0001));

   _test(closeEnough(x.getSumOfSquaresWithinGroups(0,3), 27.806983, 0.0001));
   _test(closeEnough(x.getSumOfSquaresWithinGroups(0,4), 1641.637025, 0.0001));

   _test(x.getDBetweenGroups(0,3) == 1);
   _test(x.getDBetweenGroups(0,4) == 1);

   _test(x.getDWithinGroups(0,3) == 8);
   _test(x.getDWithinGroups(0,4) == 8);

   _test(closeEnough(x.getVarianceBetweenGroups(0,3), 139.598507, 0.0001));
   _test(closeEnough(x.getVarianceBetweenGroups(0,4), 10.408335, 0.0001));

   _test(closeEnough(x.getVarianceWithinGroups(0,3), 3.475873, 0.0001));
   _test(closeEnough(x.getVarianceWithinGroups(0,4), 205.204628, 0.0001));


   _test(closeEnough(x.getFdata(0,3), 40.162143, 0.0001));
   _test(closeEnough(x.getFdata(0,4), 0.050722, 0.0001));

   _test(x.getNumberOfObservations
                   ("Screw Length",
                    "Elastic Modulus") == 10);
   _test(x.getNumberOfObservations
                   ("Screw Length",
                    1,
                    "Elastic Modulus") == 6);           
       factorValue = DataValue(2);
   _test(x.getNumberOfObservations
                   ("Screw Length",
                    factorValue,
                    "Elastic Modulus") == 4);  

   _test(closeEnough(x.getSumOfObservations
                   ("SL",
                    "EM"), 156.99, 0.0001));
   _test(closeEnough(x.getSumOfObservations
                   ("SL",
                    1,
                    "EM"), 75.89, 0.0001));           
       factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfObservations
                   ("SL",
                    factorValue,
                    "EM"), 81.1, 0.0001));  

   _test(closeEnough(x.getAverageObservation
                   ("A",
                    "D"), 15.699, 0.0001));
   _test(closeEnough(x.getAverageObservation
                   ("A",
                    1,
                    "D"), 12.648333, 0.0001));           
       factorValue = DataValue(2);
   _test(closeEnough(x.getAverageObservation
                   ("A",
                    factorValue,
                    "D"), 20.275, 0.0001));  

   _test(closeEnough(x.getSumOfSquares
                   ("a",
                    "d"), 167.405490, 0.0001));

   _test(closeEnough(x.getSumOfSquares
                   ("a",
                    1,
                    "d"), 27.378883, 0.0001));           
       factorValue = DataValue(2);
   _test(closeEnough(x.getSumOfSquares
                   ("a",
                    factorValue,
                    "d"), 0.4281, 0.0001));  

   _test(closeEnough(x.getVariance
                   (0,
                    3), 18.600610, 0.0001));
   _test(closeEnough(x.getVariance
                   (0,
                    1,
                    3), 5.475777, 0.0001));
       factorValue = DataValue(2);
   _test(closeEnough(x.getVariance
                   (0,
                    factorValue,
                    3), 0.1427, 0.0001));  

   _test(closeEnough(x.getSumOfSquaresBetweenGroups
                  ("Screw Length",
                   "EM"), 139.598507, 0.0001));

   _test(closeEnough(x.getSumOfSquaresWithinGroups
                  ("Screw Length",
                   "d"), 27.806983, 0.0001));

   _test(x.getDBetweenGroups
                  ("Screw Length",
                   "D") == 1);

   _test(x.getDWithinGroups
                  ("Screw Length",
                   3) == 8);

   _test(closeEnough(x.getVarianceBetweenGroups
                  ("SL",
                   "d"), 139.598507, 0.0001));

   _test(closeEnough(x.getVarianceWithinGroups
                  ("SL",
                   "D"), 3.475873, 0.0001));

   _test(closeEnough(x.getFdata
                  ("SL",
                   3), 40.162143, 0.0001));

}

