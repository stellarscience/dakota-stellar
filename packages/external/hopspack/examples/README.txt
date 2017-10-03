$Id: README.txt 217 2013-11-25 21:59:49Z tplante $
$URL: https://software.sandia.gov/svn/hopspack/trunk/examples/README.txt $
************************************************************************
        HOPSPACK: Hybrid Optimization Parallel Search Package
                Copyright 2009-2013 Sandia Corporation

   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
   the U.S. Government retains certain rights in this software.
************************************************************************


This file provides a brief description of the example problems.
Other subdirectories (e.g., "linked-evaluator-example") are described
by README files in the directory.

To run an example, change to its directory and execute HOPSPACK with
the *_params.txt file as input; for example:
  > cd 1-var-bnds-only
  > ../../HOPSPACK_main_serial example1_params.txt

Each example contains an executable that evaluates the objective function.
The executable must be in the execution PATH; otherwise, you may see the
error message
  ERROR: Call failed: 'var_bnds_only ...' <SystemCall>
On a Linux or Mac OSX machine this is easily fixed by adding the current
directory to the PATH environment variable:
  > export PATH=$PATH:.
For more information, please read the HOPSPACK User Manual.


Examples supplied with HOPSPACK:
  1 - minimize subject to variable bound constraints
  2 - minimize subject to variable bounds and linear constraints
  3 - same as example 2 but with an additional degenerate inequality
  4 - maximize subject to nonlinear inequality constraint
  5 - minimize a nonconvex problem using a multi-start algorithm


-- example 1 --
  The optimization problem is

    minimize    f(x1,x2) = x1^2 + (2 * x2^2)

    subject to  -1 <= x1 <= 1
                -1 <= x2 <= 1

  The exact solution point is x1=0, x2=0, with objective f(0,0) = 0.
  The serial HOPSPACK executable with given configuration file should converge
  after about 30 function evaluations.  The solution should be a point
  where each variable is within 0.01 of the exact solution, because this is
  the default GSS "Step Tolerance".  A typical solution displayed by HOPSPACK is
    x=[ 1.250e-02 -1.250e-02 ], F= 4.687e-04.

  Try changing the Step Tolerance parameter in the "Citizen 1" sublist to be:
    "Step Tolerance" double 0.002
  HOPSPACK should take more evaluations and reach a more accurate solution point.
  If you make the Step Tolerance very small, then HOPSPACK may stop because
  the Objective Target is reached.  Make this parameter zero to force a highly
  accurate search.


-- example 2 --
  The optimization problem is

    minimize    f(x[1],x[2],x[3],x[4]) = \sum (x[i] - 10)^2

    subject to  - x[1] - x[2] - x[3] - x[4] >= -10
                  x[1] - x[2] + x[3] - x[4] <=  -1
                 2x[1]        +2x[3] -7x[4] =    3
                -10 <= x[i] <= 10

  The exact solution point is (2.25, 4.6429, 2.25, 0.8571), where f = 232.416.
  Both linear inequalities are active at the solution, but none of the bounds.
  Some additional GSS configuration parameters are shown in this example;
  refer to the HOPSPACK User Manual for the complete list.


-- example 3 --
  The optimization problem is

    minimize    f(x[1],x[2],x[3],x[4]) = \sum (x[i] - 10)^2

    subject to  - x[1] - x[2] - x[3] - x[4] >= -10
                  x[1] - x[2] + x[3] - x[4] <=  -1
                 2x[1]        +2x[3]        <=   9
                 2x[1]        +2x[3] -7x[4] =    3
                -10 <= x[i] <= 10

  The exact solution point is (2.25, 4.6429, 2.25, 0.8571), where f = 232.416.
  All linear inequalities are active at the solution, but the 3rd inequality
  is weakly active.

  Example 3 is the same as example 2, except for an additional linear constraint
  that is degenerate (redundant) at the solution point.  This example shows
  how GSS is able to handle the degeneracy by automatically invoking the CDD
  library (to see details, set the "Display" to 3 in the "Citizen 1" sublist).


-- example 4 --
  The optimization problem is

    maximize    f(x,y) = x + 2y

    subject to  x + y >= 1
                1 - (x - 1)^2 - (y - 1)^2 >= 0

  The exact solution point is x = 1.44721, y = 1.89443, where f = 5.23607.
  The nonlinear inequality constraint is active at the solution, while the
  linear equality is not.

  The serial HOPSPACK executable with given configuration file should converge
  after about 400 function evaluations.

  Accuracy can be improved when there is an active nonlinear constraint,
  but multiple configuration parameters are involved.  Example 4 accuracy is
  improved by modifying these parameters:
      "Final Step Tolerance" double 1.0e-5        # (default = 0.001)
      "Penalty Parameter Increase" double 1.5     # (default = 2.0)
  The solution with default parameters  (1.403, 1.915) F = 5.233
                      then improves to  (1.444, 1.896) F = 5.236


-- example 5 --
  The optimization problem is the well-known six-humped camel back problem
  ("Towards Global Optimization", LCW Dixon and GP Szego, Eds.,
   North Holland, 1975):

    minimize    f(x,y) = x^2(4 - 2.1x^2 + x^4/3) + xy + 4y^2(y^2 - 1)

    subject to  -3 <= x <= 3
                -2 <= y <= 3

  The problem has two global minima and four local minima:
    x =  0.089842     x =  1.70361      x =  1.60711
    y = -0.712656     y = -0.796084     y =  0.568651
    f = -1.03163      f = -0.215464     f =  2.10425

    x = -0.089842     x = -1.70361      x = -1.60711
    y =  0.712656     y =  0.796084     y = -0.568651
    f = -1.03163      f = -0.215464     f =  2.10425

  Example 5 uses the HOPSPACK multi-start citizen GSS-MS to solve from 20
  different start points.  Points are chosen at random from within the bound
  constraints.  The HOPSPACK GSS solver is applied to each start point, and
  the result of each subproblem displayed.  In most cases GSS converges to
  one of the global minima.
