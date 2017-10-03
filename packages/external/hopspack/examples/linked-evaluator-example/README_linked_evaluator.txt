$Id: README_linked_evaluator.txt 220 2014-01-02 21:24:59Z tplante $
$URL: https://software.sandia.gov/svn/hopspack/trunk/examples/linked-evaluator-example/README_linked_evaluator.txt $
************************************************************************
        HOPSPACK: Hybrid Optimization Parallel Search Package
                Copyright 2009-2014 Sandia Corporation

   Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
   the U.S. Government retains certain rights in this software.
************************************************************************


Files in this directory provide an example that evaluates the optimization
objective by direct C++ call instead of the default "system call".
Compiling these files and linking with HOPSPACK allows the user to eliminate
external system calls connected by I/O files.

The problem solved is identical with example 4:

    maximize    f(x,y) = x + 2y

    subject to  x + y >= 1
                1 - (x - 1)^2 - (y - 1)^2 >= 0

  The exact solution point is x = 1.44721, y = 1.89443, where f = 5.23607.
  The nonlinear inequality constraint is active at the solution, while the
  linear equality is not.


Here is a description of the files in this directory:

  ExampleLinkedEvaluator.hpp
      - Header file for the custom evaluator class
        that derives from HOPSPACK::Evaluator.

  ExampleLinkedEvaluator.cpp
      - Source code for the custom evaluator class
        that evaluates f(x,y) and nonlinear constraints.

  CMakeLists.txt 
      - File that instructs CMake to build this class instead
        of the usual HOPSPACK evaluator.  Any dependencies on
        C++ libraries should be expressed in this file.

  HOPSPACK_EvaluatorFactory.cpp
      - Simple modification of the factory that ignores the
        Evaluator Type parameter to use ExampleLinkedEvaluator instead.


To build the example using CMake you must replace files in the HOPSPACK
source tree with files in this directory, and rebuild the desired version
of HOPSPACK.  The following steps will build the multithreaded version
of HOPSPACK on a Linux machine.  Consult the "Building HOPSPACK" section
of the User Manual for examples that build on other operating systems.

  1. Start with a clean source directory (back up the old one).
     Note that any build directories created using CMake will detect source
     file modifications; hence, it is best to also start with clean
     target build directories.

  2. Copy source files, replacing any existing files.
       CMakeLists.txt                 to  src/src-evaluator/.
       ExampleLinkedEvaluator.hpp     to  src/src-evaluator/.
       ExampleLinkedEvaluator.cpp     to  src/src-evaluator/.
       HOPSPACK_EvaluatorFactory.cpp  to  src/src-evaluator/.

  3. Create a clean target directory for the multithreaded build.
     Here it is assumed to be at the same level as hopspack-2.0.
     Change to this directory and run CMake:
     > mkdir build_link_example
     > cd build_link_example
     > cmake ../hopspack-2.0 -Dmt=yes
     > make
     > make install

  4. Run the new version of HOPSPACK using the appropriate configuration file.
     Note that configuration parameters in the "Evaluator" sublist are ignored.
     > cd examples/4-nonlinear-constraints
     > ../../HOPSPACK_main_threaded example4_params.txt
     Accuracy of the solution point can be improved by setting these input file
     parameters (see examples/README.txt):
        "Final Step Tolerance" double 1.0e-5
        "Penalty Parameter Increase" double 1.5

     HOPSPACK will run with incorrect configuration parameters (for instance,
     example 1), but the results are meaningless.
