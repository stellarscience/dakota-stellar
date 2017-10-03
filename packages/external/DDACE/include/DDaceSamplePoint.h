#ifndef DDACESAMPLEPOINT_H
#define DDACESAMPLEPOINT_H

/**
 * DDaceSamplePoint stores a vector of doubles representing a point
 * in parameter space. It also keeps its index in the table
 * of sample points, so that later we can reconstruct where to put the 
 * sample values
 */

#include <vector>
#include <iostream>


class DDaceSamplePoint
{
 public:
	DDaceSamplePoint() : index_(0), x_(0) {;}
	DDaceSamplePoint(int index, const std::vector<double>& x)
		: index_(index), x_(x)
		{;}

	int index() const {return index_;}
	int length() const {return x_.size();}
	const double& operator[](int i) const {return x_[i];}
	const std::vector<double>& parameters() const {return x_;}

//	void bcast(int sender, const PComm& comm);

	void print(std::ostream& os) const;
 private:
	int index_;
	std::vector<double> x_;
};

inline std::ostream& operator<<(std::ostream& os, const DDaceSamplePoint& pt)
{
	pt.print(os);
	return os;
}
		



#endif
