HOPSPACK 2.0.2

Thank you for using HOPSPACK, a software framework for solving
derivative-free optimization problems using parallel resources.
HOPSPACK is developed and maintained by Sandia National Laboratories,
distributed with permission of the United States Department of Energy.

This distribution contains two compiled binary executables, some
example problems, and a user manual.  The "serial" executable runs
HOPSPACK as a single thread on your machine.  The "multithread"
executable runs HOPSPACK using as many threads as specified in the
configuration parameters file (2 threads in the examples).

Binary executables should work on most 32-bit x86 Linux distributions.
They were built under Red Had Enterprise Linux WS 4 using g++ 3.4.6.
If the executable fails to run the examples provided, please email
Todd Plantenga (tplante@sandia.gov).

To run the first example with the "serial" executable, type the following
(first make sure that "." is in your $PATH environment variable):
> cd examples/1-var-bnds-only
> ../../HOPSPACK_main_serial example1_params.txt

This brief README file is only intended for getting you started.
Further information is available in
- the User Manual (HopspackUserManual_2_0_2.pdf)
- examples/README.txt
- https://software.sandia.gov/trac/hopspack


Thanks,
  Todd Plantenga, Sandia National Laboratories  (tplante@sandia.gov)
  Tamara Kolda,   Sandia National Laboratories
