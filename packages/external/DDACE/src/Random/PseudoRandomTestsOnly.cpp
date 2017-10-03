#include "PseudoRandomTestsOnly.h"

PseudoRandomTestsOnly::PseudoRandomTestsOnly() {
	this->index = 0;
}


void PseudoRandomTestsOnly::setSeed(int seed) {
	this->index = seed % 1000;
}

double PseudoRandomTestsOnly::getPseudoRandom() {
	
	/* generate a pseudo random number that is between 0 and 1 */
	double number = ((double)this->index) / 1000.0;
	
	/* increment the index so that we will generate a different number */
	/* the next time getPseudoRandom() is invoked.                 */
	this->index++;
	if (this->index >= 1000) 
          this->index=0;
	
	/* return the generated number */
	return(number);
}
