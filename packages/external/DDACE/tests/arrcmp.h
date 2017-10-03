#ifndef ARRCMP_H
#define ARRCMP_H

#include "SmartPtr.h"
#include <vector>
#include "DDaceSamplePoint.h"


int Arrcmp_d( std::vector<double>& a, std::vector<double>& b );

int Arrcmp_d_est( std::vector<double>& a, 
		std::vector<double>& b,
		const float errlim );

int Arrcmp_i( std::vector<std::vector<int> >& a, 
		std::vector<std::vector<int> >& b );

int Arrcmp_ad( std::vector<DDaceSamplePoint>& a, 
		std::vector< std::vector<double> >& b );

int Arrcmp_ad_est( std::vector<DDaceSamplePoint>& a, 
			std::vector< std::vector<double> >& b, 
			const float errlim );

int DDaceSamplePoint_cmp( const DDaceSamplePoint& a, 
				const DDaceSamplePoint& b );

// itoa is a native windows function
#ifndef _MSC_VER
void itoa( int n, char buf[], const int size );
#endif

bool closeEnough( double a, double b, const float errlim );


#endif
