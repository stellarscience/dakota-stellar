#include "TestDDaceArraySampler.h"

TestDDaceArraySampler::TestDDaceArraySampler() {
}

TestDDaceArraySampler::~TestDDaceArraySampler(){
}


void TestDDaceArraySampler::run() {


	//PMachine::init(argc, (void**)argv);

	testDDaceSamplePoint();

        //PMachine::finalize();
}


void TestDDaceArraySampler::testDDaceSamplePoint() {


	
	    /* create a sampler                        */
            /* In this example, our sampler is storing */
            /* our data samples or measurments.        */
            std::vector < std::vector < double > > array(3,std::vector<double>(2));
            array[0][0] = 1.0; array[0][1] = 2.0;
            array[1][0] = 11.0; array[1][1] = 12.0;
            array[2][0] = 21.0; array[2][1] = 22.0;	
            DDaceSampler sampler = DDaceArraySampler(array);

            /* The input values should be 1, 2, 11, 12, 21, 22 */
            std::vector< DDaceSamplePoint > data;
            sampler.getSamples(data); 
            _test(data[0][0]==1.0);
            _test(data[0][1]==2.0);
            _test(data[1][0]==11.0);
            _test(data[1][1]==12.0);
            _test(data[2][0]==21.0);
            _test(data[2][1]==22.0);


            /* get the lower and upper bounds of each input variable */
            std::vector < double > lowerBounds;
            lowerBounds = sampler.lowerBounds();
            std::vector < double > upperBounds;
            upperBounds = sampler.upperBounds();

            /* The lower bound for the first variable should be 1.0 */
            /* The upper bound for the first variable should be 21.0 */
            _test(lowerBounds[0]==1.0);
            _test(upperBounds[0]==21.0);            

            /* The lower bound for the 2nd variable should be 2.0 */
            /* The upper bound for the 2nd variable should be 22.0 */
            _test(lowerBounds[1]==2.0);
            _test(upperBounds[1]==22.0);            
	    
            /* Create a name for every input variable in our samples */
	    std::vector<std::string> variableNames(2);
	    variableNames[0] = "x1";
	    variableNames[1] = "x2";

            /* Create a name for every output variable in our samples */
	    std::vector<std::string> outputNames(1, "y(x1,x2)");

//            /* encapsulate everything into a DDace object */
//	    DDace ddace
//	        (sampler,      //contains (x,y) pair of input data
//	         variableNames, //names of all variables
//	         outputNames,  //name of all user-defined functions
//	         "userInputArchive.xml");  //name of output file
//
//            PMachine::synchronize(PWorld());
//
//            /* extract the input data */            
//            Array< DDaceSamplePoint > dataInput;
//            Array< Array < double > > dataOutput;
//            ddace.getResults(dataInput, dataOutput);
//
//            PMachine::synchronize(PWorld());
//
//            /* The input data should have 3 samples or runs */
//            if (PMachine::getRank()==0) {
//                int numberOfSamples = dataInput.length();
//                _test(numberOfSamples==3);
//            }
//
//            PMachine::synchronize(PWorld());
//
//            /* Each sample or run should contain values for 2 input variables */
//            if (PMachine::getRank()==0) {
//                int numberOfSamples = dataInput.length();
//                int numberOfInputVariables = 0;
//                if (numberOfSamples > 0)
//                    numberOfInputVariables = dataInput[0].length();
//
//                _test(numberOfInputVariables==2);
//                _test(dataInput[0][0]==1.0);
//                _test(dataInput[0][1]==2.0);
//                _test(dataInput[1][0]==11.0);
//                _test(dataInput[1][1]==12.0);
//                _test(dataInput[2][0]==21.0);
//                _test(dataInput[2][1]==22.0);
//            } 
  
            /* The name of the sampler should be a DDaceArraySampler */
//            PMachine::synchronize(PWorld());
            _test(!sampler.typeName().compare("DDaceArraySampler"));        
             
           
//            /* Get the sampler from the ddace object */ 
//	    PMachine::synchronize(PWorld());
//            if (PMachine::getRank()==0) {
//                DDaceSampler ddaceSampler;
//                ddace.getSampler(ddaceSampler);
//
//
//                /* The input values should be 1, 2, 11, 12, 21, 22 */
//                ddaceSampler.getSamples(data); 
//                _test(data[0][0]==1.0);
//                _test(data[0][1]==2.0);
//                _test(data[1][0]==11.0);
//                _test(data[1][1]==12.0);
//                _test(data[2][0]==21.0);
//                _test(data[2][1]==22.0);
//            }

 //           PMachine::synchronize(PWorld());

}


;

