//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <limits>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/LinearLagrangeInterpolationSurrogate.h>
#include <queso/InterpolationSurrogateData.h>

double two_d_fn( double x, double y );

int main(int argc, char ** argv)
{
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc, &argv);
  QUESO::FullEnvironment env(MPI_COMM_WORLD, "test_InterpolationSurrogate/queso_input.txt", "", NULL);
#else
  QUESO::FullEnvironment env("test_InterpolationSurrogate/queso_input.txt", "", NULL);
#endif

  int return_flag = 0;

  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix>
    paramSpace(env,"param_", 2, NULL);

  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins[0] = -2.5;
  paramMins[1] = 3.0;

  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs[0] = 1.4;
  paramMaxs[1] = 4.1;

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMins, paramMaxs);

  std::vector<unsigned int> n_points(2);
  n_points[0] = 101;
  n_points[1] = 51;

  QUESO::InterpolationSurrogateData<QUESO::GslVector, QUESO::GslMatrix>
    data(paramDomain,n_points);

  std::vector<double> values(n_points[0]*n_points[1]);

  double spacing_x = (paramMaxs[0] - paramMins[0])/(n_points[0]-1);
  double spacing_y = (paramMaxs[1] - paramMins[1])/(n_points[1]-1);

  for( unsigned int i = 0; i < n_points[0]; i++ )
    {
      for( unsigned int j = 0; j < n_points[1]; j++ )
        {
          unsigned int n = i + j*n_points[0];

          double x = paramMins[0] + i*spacing_x;
          double y = paramMins[1] + j*spacing_y;

          values[n] = two_d_fn(x,y);
        }
    }

  data.set_values( values );

  QUESO::LinearLagrangeInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>
    two_d_surrogate( data );

  QUESO::GslVector domainVector(paramSpace.zeroVector());
  domainVector[0] = -0.4;
  domainVector[1] = 3.764;

  double test_val = two_d_surrogate.evaluate(domainVector);

  double exact_val = two_d_fn(domainVector[0],domainVector[1]);

  double tol = 2.0*std::numeric_limits<double>::epsilon();

  double rel_error = (test_val - exact_val)/exact_val;

  if( std::fabs(rel_error) > tol )
    {
      std::cerr << "ERROR: Tolerance exceeded for 2D Lagrange interpolation test."
                << std::endl
                << " test_val  = " << test_val << std::endl
                << " exact_val = " << exact_val << std::endl
                << " rel_error = " << rel_error << std::endl
                << " tol       = " << tol << std::endl;

      return_flag = 1;
    }

#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return return_flag;
}

double two_d_fn( double x, double y )
{
  return 3.0 + 2.5*x - 3.1*y +0.1*x*y;
}
