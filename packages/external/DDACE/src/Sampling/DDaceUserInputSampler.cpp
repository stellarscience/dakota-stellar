#include "DDaceUserInputSampler.h"

using namespace std;

std::string DDaceUserInputSampler::typeName_ = "DDaceUserInputSampler";


DDaceUserInputSampler::DDaceUserInputSampler(const string& ptFilename)
	: DDaceSamplerBase(0, 0, false, std::vector<Distribution>(0)),
		ptFilename_(ptFilename), pts_(0), lowerBounds_(), upperBounds_()
{
	ifstream in(ptFilename_.c_str());
	if (!in) 
		{
			throw runtime_error("DDaceUserInputSampler ctor: could not open input file");
		}
	
	std::vector<std::vector<string> > lines = tokenizeFile(in, '#');
	
	nSamples_ = lines.size();
	nInputs_ = lines[0].size();
	
	pts_.resize(nSamples_);
	upperBounds_.resize(nInputs_);
	lowerBounds_.resize(nInputs_);

	for (int j=0; j<nInputs_; j++)
		{
			lowerBounds_[j] = 1.0e50;
			upperBounds_[j] = -1.0e50;
		}

	for (int i=0; i<nSamples_; i++)
		{
			if ((int) lines[i].size() != nInputs_)
				{
					throw runtime_error("DDaceUserInputSampler(): mismatched input line lengths");
				}
			std::vector<double> p(nInputs_);
			for (int j=0; j<nInputs_; j++)
				{
					p[j] = atof(lines[i][j].c_str());
					if (p[j] < lowerBounds_[j]) lowerBounds_[j] = p[j];
					if (p[j] > upperBounds_[j]) upperBounds_[j] = p[j];
				}
			pts_[i] = DDaceSamplePoint(i, p);
		}
}

const std::vector<Distribution>& DDaceUserInputSampler::dist() const
{
	throw runtime_error("DDaceUserInputSampler::dist() should not be called");
	return dist_;
}

std::vector<double> DDaceUserInputSampler::lowerBounds() const
{
	return lowerBounds_;
}

std::vector<double> DDaceUserInputSampler::upperBounds() const
{
	return upperBounds_;
}


DDaceSamplerBase* DDaceUserInputSampler::clone() const
{
	DDaceSamplerBase* rtn = new DDaceUserInputSampler(*this);
	if (rtn==0) throw std::bad_alloc();
	return rtn;
}

void DDaceUserInputSampler::print(ostream& os) const
{
	os << "<UserInputSampler filename=\"" << ptFilename_ 
	   << "\" samples=\"" << nSamples_ << "\"/>" ;
}



int DDaceUserInputSampler::getParameter(const string& /*parameterName*/) const 
{
	throw runtime_error("DDaceUserInputSampler::getParameter class has no parameters");
	return 0;
}

std::vector<std::vector<std::string> > DDaceUserInputSampler::tokenizeFile(
							std::istream& is,
                                                        char comment)
{
        std::string line;
        char cstr[500];
        std::vector<std::vector<std::string> > rtn(0);
        std::vector<std::string> lines;
        while(!is.eof()) {
                is.getline(cstr,499);
                lines.push_back(std::string(cstr));
        }
                                                                                
        rtn.reserve(lines.size());
        int count = 0;
        for (int i=0; i< (int) lines.size(); i++)
                {
                        if (lines[i].length() == 0) continue;
                        std::vector<std::string> tokens = stringTokenizer(lines[i]);
                        if (tokens.size() == 0) continue;
                        rtn.push_back(tokens);
                        count++;
                }
                                                                                
        return rtn;
}


/** The StringTokenizer takes a String in and returns an Array of Strings
 *  that are composed of the individual words of the original String. It
 *  uses the findNextNonWhitespace and findNextWhitespace functions to
 *  determine where the words begin and end.
 *  In the days of yore, this was in StrUtils, a class in the trecherous
 *  CPPUtilities library.
 */
std::vector<std::string> DDaceUserInputSampler::stringTokenizer(const std::string& str)
{
        std::vector<std::string> rtn(0);
        int start = 0;
        while(start < (int) str.length())
                {
                        start =  findNextNonWhitespace(str, start);
                        int stop = findNextWhitespace(str, start);
                        if (start-stop == 0) return rtn;
                        std::string sub = str.substr(start, stop);
                        rtn.push_back(sub);
                        start =  findNextNonWhitespace(str, stop);
                }
        return rtn;
}


/** This function takes a String and an integer value as its input
  * The integer value is used as the starting position of the search.
  * It then begins looking for spaces, tabs, carriage returns,
  * or line breaks.
  * Once upon a time, this was in StrUtils, before the liberation from
  * CPPUtilities.
  */

int DDaceUserInputSampler::findNextWhitespace(const std::string& str, 
						int offset)
{
        for (int i=0; i< (int) (str.length()-offset); i++)
                {
                        if (str[i+offset]==' ' ||
                            str[i+offset]=='\t' ||
                            str[i+offset]=='\n' ||
                            str[i+offset]=='\r')
                                                                                
                                {
                                        return i+offset;
                                }
                }
        return str.length();
}

int DDaceUserInputSampler::findNextNonWhitespace(const std::string& str,
						int offset)
{
        for (int i=0; i<( (int) str.length()-offset); i++)
                {
                        if (!(str[i+offset]==' ' ||
                                str[i+offset]=='\t' ||
                                str[i+offset]=='\n' ||
                                str[i+offset]=='\r'))
                                {
                                        return i+offset;
                                }
                }
        return str.length();
}

std::vector<std::vector<int> > DDaceUserInputSampler::getP() const 
{
        throw runtime_error("DDaceSamplerBase::getP not defined for base class");
        std::vector<std::vector<int> > tmp;
        return tmp;
}


