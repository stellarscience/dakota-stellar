#ifdef HAVE_CONFIG
#include "fsudace_config.h"
#endif
# ifdef HAVE_STD
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

int main ( int argc, char *argv[] );
void quality_handle ( char *input_filename );

//****************************************************************************

int main ( int argc, char *argv[] )

//****************************************************************************
//
//  Purpose:
//
//    FSU_TABLE_QUALITY determines quality measures for a given set of points.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Modified:
//
//    06 October 2004
//
//  Author:
//
//    Max Gunzburger
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Local parameters:
//
//    Local, int NDIM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[NDIM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
{ 
  int i;
  char input_filename[81];

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "FSU_TABLE_QUALITY:\n";
  cout << "  Compute measures of uniform dispersion for a pointset.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "FSU_TABLE_QUALITY:\n";
    cout << "  Please enter the name of a file to be analyzed.\n";

    cin.getline ( input_filename, sizeof ( input_filename ) );

    quality_handle ( input_filename );

  }
  else 
  {
    for ( i = 1; i < argc; i++ ) 
    {
      quality_handle ( argv[i] );
    }
  } 

  cout << "\n";
  cout << "FSU_TABLE_QUALITY:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

void quality_handle ( char *input_filename )

//****************************************************************************
//
//  Purpose:
//
//    QUALITY_HANDLE handles a single file.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Modified:
//
//    04 October 2004
//
//  Author:
//
//    Max Gunzburger
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//  Local parameters:
//
//    Local, int NDIM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[NDIM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
{
  int n;
  int ndim;
  int ns = 100000;
  int seed_init = 123456789;
  double *z; 

  dtable_header_read ( input_filename, &ndim, &n );
// 
//  Read the point set.
//
  z = dtable_data_read ( input_filename, ndim, n );

  cout << "\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\".\n";
  cout << "\n";
  cout << "  The spatial dimension NDIM =     " << ndim      << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n" << flush;
  cout << "\n";

  cout << "  The regularity measure         Chi = "
    << chi_measure ( ndim, n, z, ns, seed_init ) << "\n";
  cout << "  2nd moment determinant measure   D = "
    << d_measure ( ndim, n, z, ns, seed_init ) << "\n";
  cout << "  The point distribution norm      H = "
    << h_measure ( ndim, n, z, ns, seed_init ) << "\n";
  cout << "  2nd moment trace measure       Tau = "
    << tau_measure ( ndim, n, z, ns, seed_init ) << "\n";

  delete [] z;

  return;
}
