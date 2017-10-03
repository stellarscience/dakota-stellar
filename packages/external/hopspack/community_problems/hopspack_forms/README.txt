Todd Plantenga notes, especially on converting from APPSPACK to HOPSPACK.
June 2009.


----- Directories ----------------------------------------------------

"case_" directories are the original Community paper formulation of each
test problem.  Linear and nonlinear constraints are checked at evaluation,
and the objective returned as invalid if infeasible.

TBD "lincon_" directories incorporate the linear constraints into the HOPSPACK
description.

TBD "nonlin_" directories incorporate TBD linear constraints into the
HOPSPACK description, and return nonlinear constraints separate from the
objective so they can be used in a penalty method.


----- Building -------------------------------------------------------

1. Run make in modflw96/src to create FORTRAN executable ../bin/modflw96.
   Put it in $PATH (or symbolically link from all evaluation directories).
   - The HC problem uses modflw96lenx instead, but I don't know how to build it
     (Genetha thinks you have to modify the resolution parameter internally).
2. cd to one of the problem evaluation directories.
3. Run "make" to create evalXXX.
4. Run "make test" and compare test_output.1 with test_expected_output.1.
   The test point was usually chosen as the first APPSPACK evaluation
   point that did not have a constraint violation.


----- Running --------------------------------------------------------

The usual HOPSPACK calls will work.  For example:
  HOPSPACK_main_serial  hopspack_parameters.txt


----- Notes ----------------------------------------------------------

The only change I had to make to the evaluation code in the "case_"
directories was to read the evaluation request type, and write the length
of the output vector (both in file EvalFn.c).

About the problems:
- Case UNC5 is in the Community paper.
  The goal is to find positions for 5 wells (10 variables total).
- Case CON5 is not in the paper.
  Seems to be a variant of UNC5 hard-coded parameters, with an otherwise
  identical structure.  Genetha Gray says UNC means "unconfined flow" and
  CON means "confined flow".
- Case UNC6 is in the Community paper.
  The goal is to find positions and flow rates for 6 wells (18 variables total).
  Optimization should find that one of the wells has nearly zero flow.
- Case CON6 is not in the paper.
  Seems to be a variant of UNC6 hard-coded parameters, with an otherwise
  identical structure.
- Case HC is in the Community paper.
  The goal is to find positions and flow rates for 4 wells (12 variables total).
  An additional "gradient" constraint is included with the C3 check made
  by the FORTRAN evaluation code.

Solutions tend to look different when run asynchronously.  Tammy says the
objective functions are shallow, so locations can move without affecting the
objective value.  Also, there is no constraint to fix a particular well at a
particular location, and sometimes two identical answers appear different
because a couple of wells swapped location.

The APPSPACK Community code worked as follows:
- APPSPACK Evaluator makes a system call to the shell script APPScomUNC5
  - Creates and populates files for MODFLOW data (eg, A2-nnn.mfn)
  - Copies the command line arguments (file names) to input.nnn
  - Calls FORTRAN program compiled from EvalFn
    - Reads the point's location
    - Checks linear constraint feasibility (C1 and C2), aborting if infeasible
    - Makes a system call to modflw96
      - Performs PDE computation, outputs are in a file
    - Reads the result
    - Checks nonlinear constraint feasibility (C3), aborting if infeasible
    - Writes the result in a file for APPSPACK
  - Deletes temp files

Genetha Gray still works with the Community problems, but uses MODFLOW 2000.
Her team uses GA+GSS and TGP.  They do not incorporate constraints, but modify
the "infinite objective" to be something more like a barrier value.  They use
GA to formulate HC with integer variables, which makes it easier and removes
the ad-hoc threshold.
