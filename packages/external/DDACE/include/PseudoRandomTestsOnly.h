#ifndef PSEUDORANDOMTESTSONLY_
#define PSEUDORANDOMTESTSONLY_

#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
#include "ddace_config.h"
#endif

#include <cstdlib>
#include <cstdio>


class PseudoRandomTestsOnly {


    protected: 
        int index;
        
        
    public:
    
    
        /**
         * Construct a PseudoRandomTestsOnly 
         */
        PseudoRandomTestsOnly();
    
       /**
        * Get a fake number that is between 0 and 1.0, 
        * excluding 1.0. 
        * @return fake number that is between 0 and 1.0, excluding 1.0.
        */
       double getPseudoRandom();
       
       
       /**
        * set the seed
        * @param seed Set the seed for this fake random number generator.
        */
        void setSeed(int seed);
         
    
};



#endif 
