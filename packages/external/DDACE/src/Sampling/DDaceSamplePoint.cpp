#include "DDaceSamplePoint.h"


void DDaceSamplePoint::print(std::ostream& os) const 
{
	os << "[ " << index_ << " ";

	os << "(";
	int ultimate = x_.size()-1;

	for (int i = 0; i < ultimate; ++i)
		os << x_[i] << ", ";

	if(ultimate >= 0) os << x_[ultimate] << ") ]";
}

