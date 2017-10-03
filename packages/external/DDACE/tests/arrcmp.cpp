#include "arrcmp.h"
//#include "DDaceSamplePoint.h"
#include <math.h>

///////////////////////////
// Written By: J Cramp
// Date: May 2004
///////////////////////////

int Arrcmp_d( std::vector<double>& a, std::vector<double>& b )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b )
  {
    return 1;  // failed!
  }

  for( int i = 0; i < len_a; i++ )
  {
      if( a[i] != b[i] )
      { 
        return 1;  // failed!
      }
  }
 
  return 0; // success!
}

int Arrcmp_d_est( std::vector<double>& a, 
			std::vector<double>& b, 
			const float errlim )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b )
  {
    return 1;  // failed!
  }

  for( int i = 0; i < len_a; i++ )
  {
      //if( abs(a[i] - b[i]) >= errlim )
      //{ 
      //  return 1;  // failed!
      //}
      if (a[i] >= b[i]) {
          if ((a[i] - b[i]) >= errlim)
              return 1; // failed!
      } else {
          if ((b[i] - a[i]) >= errlim)     
              return 1; // failed!
      }
      
      
  }
 
  return 0; // success!
}

int Arrcmp_i( std::vector<std::vector<int> >& a,
		 std::vector<std::vector<int> >& b )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b )
  {
    return 1;  // failed!
  }

  for( int i = 0; i < len_a; i++ )
  {
    if( a[i].size() != b[i].size() )
    {
      return 1;  // failed!
    }
    for( int j = 0; j < (int) a[i].size(); j++ )
    {
      if( a[i][j] != b[i][j] )
      {
        return 1;  // failed!
      }
    }
  }

  return 0;  // success!
}

int Arrcmp_ad( std::vector<DDaceSamplePoint>& a, 
		std::vector< std::vector<double> >& b )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b )
  {
    return 1;  // failed!
  }

  for( int i = 0; i < len_a; i++ )
  {
     if( a[i].length() != (int) b[i].size() )
     {
        return 1;  // failed!
     }
     for( int j = 0; j < a[i].length(); j++ )
     {
        if( a[i][j] != b[i][j] )
        { 
          return 1;  // failed!
        }
     }
  }
 
  return 0; // success!
}

int Arrcmp_ad_est( std::vector<DDaceSamplePoint>& a,
			std::vector< std::vector<double> >& b, 
			const float errlim )
{
  int len_a = a.size();
  int len_b = b.size();

  if( len_a != len_b )
  {
    return 1;  // failed!
  }

  for( int i = 0; i < len_a; i++ )
  {
     if( a[i].length() != (int) b[i].size() )
     {
        return 1;  // failed!
     }
     for( int j = 0; j < a[i].length(); j++ )
     {
        if( !closeEnough( a[i][j], b[i][j], errlim ) )
        { 
          return 1;  // failed!
        }
     }
  }
 
  return 0; // success!
}

int DDaceSamplePoint_cmp( const DDaceSamplePoint& a, const DDaceSamplePoint& b )
{
  int len_a = a.length();
  int len_b = b.length();

  if( len_a != len_b )
  {
    return 1;  // Failed! sample points not the same
  }

  for( int i = 0; i < len_a; i++ )
  {
     if( a[i] != b[i] )
     {
        return 1;  // Failed! sample points not the same
     }
  }

  return 0;  // Success! sample points are the same
}

// itoa is a native windows function
#ifndef _MSC_VER 
void itoa( int n, char buf[], const int size )
{
  int    c     = 0;
  int    tmp;
  double i;
  const char*  alpha = "0123456789";

  // determine highest power of 10 needed
  i = -1;
  while( (n/pow(10.0,++i)) > 1 )
      ;

  // i is one too big
  i--;

  // map each digit to its character
  for( ; i >= 0 && c < size; i-- )
  {
     tmp = (n/((int)pow(10.0,i)));
     buf[c++] = alpha[tmp];
     n -= (tmp * ((int)pow(10.0,i)));
  }

  // insert null character
  buf[c] = '\0';
}
#endif

bool closeEnough( double a, double b, const float errlim )
{
  if( a > b )
  {
    return (a-b < errlim);
  }
  return (b-a < errlim );
}

  

