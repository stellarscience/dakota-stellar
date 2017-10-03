#include "DDaceArraySampler.h"

using namespace std;

string DDaceArraySampler::typeName_ = "DDaceArraySampler";

DDaceArraySampler::DDaceArraySampler(std::vector< std::vector < double > > array)
	: DDaceSamplerBase(0, 0, false, std::vector<Distribution>(0)),
	 pts_(0), lowerBounds_(), upperBounds_()
{
    this->setInputData(array);
}




void DDaceArraySampler::setInputData(std::vector < std::vector < double > >&array) {



	/* number of data samples or experiments      */
        /* A data sample contains one value for every */
        /* input variable (x1, x2, x3, ...)           */
	nSamples_ = array.size();

        /* The number of input variables (x1,x2,x3,...) */
	nInputs_ = array[0].size();
	
        /* We will copy the input and output values into this array */
	pts_.resize(nSamples_);

        /*
         * For each input variable (x1 or x2 or x3 or ...)
         * we want to find the range of values (i.e. the smallest
         * value for the variable and the largest value for the 
         * variable).
         * NOTE:  lowerBound_[n] contains the lowest value
         * of the nth input variable.
         * We are going to initialize lowerBound_[n] to the
         * the FIRST value of the nth input variable.  Later,
         * we will adjust this assumption.
         */
	upperBounds_.resize(nInputs_);
	lowerBounds_.resize(nInputs_);
	for (int j=0; j<nInputs_; j++) {
	    lowerBounds_[j] = upperBounds_[j] = 0.0;
            if (nSamples_>0)
               lowerBounds_[j] = upperBounds_[j] = array[0][j];
	}

        /*
         * For each experiment or run, we measured ONE value
         * for each of the input variables (x1, x2, x3, ...).
         * Collect these values into a DDaceSamplePoint.
         * That is, DDaceSamplePoint[n] contains the values
         * we measured for x1, x2, x3,... during the nth run.
         *
         * ALSO, if a value for a variable is lower than lowerBound_
         * then reset the lowerBound.  Ditto for upperBound_.
         */
         for (int i=0; i<nSamples_; i++) {
             if ((int) array[i].size() != nInputs_){
                     throw std::runtime_error("DDaceArraySampler(): mismatched input line lengths");
             }
             std::vector<double> p(nInputs_);
             for (int j=0; j<nInputs_; j++) {
                 p[j] = array[i][j];
                 if (p[j] < lowerBounds_[j]) lowerBounds_[j] = p[j];
                 if (p[j] > upperBounds_[j]) upperBounds_[j] = p[j];
             }
             pts_[i] = DDaceSamplePoint(i, p);
         }

}


const std::vector<Distribution>& DDaceArraySampler::dist() const
{
	throw runtime_error("DDaceArraySampler::dist() should not be called");
	return dist_;
}

std::vector<double> DDaceArraySampler::lowerBounds() const
{
	return lowerBounds_;
}

std::vector<double> DDaceArraySampler::upperBounds() const
{
	return upperBounds_;
}


DDaceSamplerBase* DDaceArraySampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceArraySampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}

void DDaceArraySampler::print(ostream& os) const
{
	os << "<ArraySampler " 
	   << "\" samples=\"" << nSamples_ << "\"/>" ;
}



int DDaceArraySampler::getParameter(const string& /*parameterName*/) const 
{
	throw runtime_error("DDaceArraySampler::getParameter class has no parameters");
	return 0;
}
std::vector<std::vector<int> > DDaceArraySampler::getP() const
{
        throw runtime_error("DDaceArraySampler::getP not defined for this class");
        std::vector<std::vector<int> > tmp;
        return tmp;
}

