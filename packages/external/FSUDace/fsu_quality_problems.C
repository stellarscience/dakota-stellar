#ifdef HAVE_CONFIG
#include "fsudace_config.h"
#endif
# ifdef HAVE_STD
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
#   include <cstring>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
#   include <string.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

int main ( void );
void quality_test01 ( int ndim, int n, double z[], int ns, int seed_init );
void quality_test02 ( int ndim, int n, double z[], int ns, int seed_init );
void quality_test03 ( int ndim, int n, double z[], int ns, int seed_init );
void quality_test04 ( int ndim, int n, double z[], int ns, int seed_init );

//******************************************************************************

int main ( void )

//******************************************************************************
//
//  Purpose:
//
//    FSU_QUALITY_PROBLEMS calls the FSU_QUALITY routines.
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
{
  char input_filename[80];
  int n;
  int ndim;
  int ns;
  int seed_init;
  double *z;

  timestamp ( );

  cout << "\n";
  cout << "FSU_QUALITY_PROBLEMS\n";
  cout << "  Test the C++ FSU_QUALITY library.\n";

  ns = 100000;
  seed_init = 123456789;
  strcpy ( input_filename, "halton_02_00100.txt" );

  dtable_header_read ( input_filename, &ndim, &n );

  cout << "\n";
  cout << "FSU_QUALITY_PROBLEMS:\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\"\n";
  cout << "\n";
  cout << "  The spatial dimension NDIM =     " << ndim      << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n";
  cout << "\n";

  z = dtable_data_read ( input_filename, ndim, n );

  dmat_transpose_print_some ( ndim, n, z, 1, 1, 5, 5, 
    "  5x5 portion of data read from file:" );

  quality_test01 ( ndim, n, z, ns, seed_init );
  quality_test02 ( ndim, n, z, ns, seed_init );
  quality_test03 ( ndim, n, z, ns, seed_init );
  quality_test04 ( ndim, n, z, ns, seed_init );

  delete [] z;

  cout << "\n";
  cout << "FSU_QUALITY_PROBLEMS\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//******************************************************************************

void quality_test01 ( int ndim, int n, double z[], int ns, int seed_init )

//******************************************************************************
//
//  Purpose:
//
//    QUALITY_TEST01 tests CHI_MEASURE.
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
{
  cout << "\n";
  cout << "QUALITY_TEST01\n";
  cout << "  CHI_MEASURE computes the CHI measure of quality.\n";
  cout << "  The regularity measure         Chi = "
       << chi_measure ( ndim, n, z, ns, seed_init ) << "\n";

  return;
}
//******************************************************************************

void quality_test02 ( int ndim, int n, double z[], int ns, int seed_init )

//******************************************************************************
//
//  Purpose:
//
//    QUALITY_TEST02 tests D_MEASURE.
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
{
  cout << "\n";
  cout << "QUALITY_TEST02\n";
  cout << "  D_MEASURE computes the D measure of quality.\n";
  cout << "  2nd moment determinant measure   D = "
       << d_measure ( ndim, n, z, ns, seed_init ) << "\n";

  return;
}
//******************************************************************************

void quality_test03 ( int ndim, int n, double z[], int ns, int seed_init )

//******************************************************************************
//
//  Purpose:
//
//    QUALITY_TEST03 tests H_MEASURE.
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
{
  cout << "\n";
  cout << "QUALITY_TEST03\n";
  cout << "  H_MEASURE computes the H measure of quality.\n";
  cout << "  The point distribution norm      H = "
       << h_measure ( ndim, n, z, ns, seed_init ) << "\n";

  return;
}
//******************************************************************************

void quality_test04 ( int ndim, int n, double z[], int ns, int seed_init )

//******************************************************************************
//
//  Purpose:
//
//    QUALITY_TEST04 tests TAU_MEASURE.
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
{
  cout << "\n";
  cout << "QUALITY_TEST04\n";
  cout << "  TAU_MEASURE computes the Tau measure of quality.\n";
  cout << "  2nd moment trace measure       Tau = "
       << tau_measure ( ndim, n, z, ns, seed_init ) << "\n";

  return;
}
