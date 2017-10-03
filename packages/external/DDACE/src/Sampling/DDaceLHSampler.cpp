#include "DDaceLHSampler.h"

using namespace std;

std::string DDaceLHSampler::typeName_ = "DDaceLHSampler";

DDaceLHSampler::DDaceLHSampler(int nSamples, int nReplications, bool noise,
		const std::vector<Distribution>& dist)
		: DDaceSamplerBase(nSamples, dist.size(), noise, dist), 
		permutationMatrix_(nSamples, vector<int>((int) dist.size())),
		nSymbols_(0), 
		nReplications_(nReplications)
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols

	nSymbols_ = nSamples_ / nReplications_;

    /**
     * Check that nInputs equals dist.size(), this must be true or else
     * there is the possibility of a seg fault in the method getSamples().
     */
    if (nInputs_ != (int) dist.size())
        throw runtime_error("DDaceLHSampler: nInputs not equal to dist.length()");


	// set up sample pattern
	initPattern();
}

DDaceLHSampler::DDaceLHSampler(int nSamples, int nInputs, int nReplications, bool noise)
	: DDaceSamplerBase(nSamples, nInputs, noise,
	 std::vector<Distribution>(nInputs, UniformDistribution(0,10))), 
		permutationMatrix_(nSamples, vector<int>(nInputs)), nSymbols_(0),
		nReplications_(nReplications)
{
	// The following code is adapted from Charles Tong's original
	// implementation

	// compute number of symbols

	nSymbols_ = nSamples_ / nReplications_;

	// set up sample pattern
	initPattern();
}

void DDaceLHSampler::initPattern()
{
	int i, j, k;		
	// randomize the pattern

	for (i=0; i<nSamples_; i+=nSymbols_)
		{
			for (j=0; j<nSymbols_; j++)
				{
					for (k=0; k<nInputs_; k++)
						{
							permutationMatrix_[i+j][k] = j;
						}
				}
		}


	vector<int> ivec2(nSamples_);
	vector<int> ivec1(nSymbols_);
	
	for (i=0; i<nSamples_; i+=nSymbols_)
		{
			for (j=0; j<nInputs_; j++)
				{
                                        //cout << "setting up ivec1 " << endl;
					ivec1 = randomIVector(nSymbols_);
				        //cout << "we have ivec1 " << endl;
					for (k=0; k<nSymbols_; k++)
						{
							ivec2[k] = permutationMatrix_[i+ivec1[k]][j];
						}
					for (int k=0; k<nSymbols_; k++)
						{
							permutationMatrix_[i+k][j] = ivec2[k];
						}
				}
		}
}

/**
 * If a constructor which doesn't set the distribution is used to create the
 * DDaceLHSampler object the getSamples() method will throw an
 * ArgumentMisMatchException unless the setDist() method is used first to
 * pass in an array of distributions of length nInputs_.
 */

vector<DDaceSamplePoint>& DDaceLHSampler::getSamples(vector<DDaceSamplePoint>& samplePoints) const
{
	int i, s;
	double denom = (double) nSymbols_;
	samplePoints.resize(nSamples_);
	vector<double> x(nInputs_);

	for (s=0; s<nSamples_; s++)
		{
			for (i=0; i<nInputs_; i++)
				{
					double p = (double) permutationMatrix_[s][i];
					if (noise_)
						{
							p = (p + DistributionBase::uniformUnitDeviate())/denom;
						}
					else
						{
							p = (p + 0.5)/denom;
							//cout << "p = " << p << endl;
						}
					// compute x as CDF^-1(p).
					x[i] = dist_[i].getDeviate(p);					
				}
			samplePoints[s] = DDaceSamplePoint(s, x);
		}

	return samplePoints;
}
	
DDaceSamplerBase* DDaceLHSampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceLHSampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}


void DDaceLHSampler::print(ostream& os) const
{
	os << "<LatinHypercube "
		 << "samples=\"" << nSamples_ << "\" ";
	os << "replications=\"" << nReplications_ << "\" ";
	os << "perturb=\"";
	if (noise_)
		{
			os << "true\" ";
		}
	else
		{
			os << "false\" ";
		}
	os << "seed=\"" << DistributionBase::seed() << "\"/>";
}

int DDaceLHSampler::getParameter(const std::string& parameterName) const 
{
	std::string tmp(parameterName);
	transform(tmp.begin(),tmp.end(),tmp.begin(),(int(*) (int)) toupper);
	if (tmp == "REPLICATIONS")
		{
			return nReplications_;
		}
	throw runtime_error("DDaceLHSampler::getParameter invalid parameter name");
}

