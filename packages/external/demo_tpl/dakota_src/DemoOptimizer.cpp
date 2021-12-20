/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       DemoTPLOptimizer
//- Description: Wrapper class for Demo_Opt
//- Owner:       Russell Hooper
//- Checked by:  

// Dakota headers
// This header file provides the database that holds options specified
// in the Dakota input file.

#include "ProblemDescDB.hpp"

// Demo_Opt headers
// There are two sets of source files needed for integrating a TPL
// into Dakota: the source code for the TPL and the source files that
// provide the interface between Dakota and the TPL.  This source file
// and the first header comprise the interface.  The second header is
// associated the the TPL source.  Replace the second header with any
// necessary header files from the TPL that is being integrated.

#include "DemoOptimizer.hpp"
#include "demo_opt.hpp"

//
// - DemoTPLOptimizer implementation
//

// -----------------------------------------------------------------
// Some calls to the Optimizer (Demo_Opt) are unique to that particular
// TPL, whereas other calls, eg to exchange data, employ general
// utilities such as adapters contained within a helper class.
//
// Calls specific to the TPL are indicated via 
//
//                      TPL_SPECIFIC
//
// within the associated comments.
// -----------------------------------------------------------------



// All of the interface source should be included in the Dakota
// namespace.

namespace Dakota {


// -----------------------------------------------------------------
/** Implementation of DemoTPLOptimizer class. */

// Standard constructor for DemoTPLOptimizer.  Sets up Demo_Opt solver
// based on information from the database.  problem_db and model are
// Dakota objects from which various information can be accessed.  No
// additional implementation is needed in the wrapper.  DemoOptTraits
// is also a Dakot object, but it requires some additional
// implementation; instructions are in DemoOptimizer.hpp.

// Demo_Opt and demoOpt are specific to the TPL being integrated.  In
// this example, the assumption is that the TPL is object oriented and
// that its main class (Demo_Opt) has a function evaluation method
// (ObjectiveFn()) that needs to be implemented.  The TPL solver
// (demoOpt) is instantiated.  TPLs do not have to be in C++, nor to
// the necessary objects have to be instantiated/initialized in the
// constructor as long as they are created before they are accessed.

DemoTPLOptimizer::DemoTPLOptimizer(ProblemDescDB& problem_db, Model& model):
  Optimizer(problem_db, model, std::shared_ptr<TraitsBase>(new DemoOptTraits())),
  Demo_Opt::ObjectiveFn(),
  Demo_Opt::NonlinearEqFn(),
  Demo_Opt::NonlinearIneqFn(),
  demoOpt(std::make_shared<Demo_Opt>())
{
  // Call a helper function to set method parameters.  It is
  // implemented later in this source file.

  set_demo_parameters();

  // Register ourself as the callback interface for objective function
  // evaluations and nonlinear equality constraints.  This assumes that
  // the TPL makes a function call to do objective function evaluations,
  // and a pointer to the function must be provided to the TPL.  This code
  // should be replaced with whatever mechanism the TPL being integrated
  // uses for setting that function pointer.  There are other ways that
  // objective functions can be implemented that will be added to future
  // versions of this example.

  // ------------------  TPL_SPECIFIC  ------------------
  demoOpt->register_obj_fn(this);

  // Conditionally register ourself as a constraint callback depending
  // on if the problem uses them
  if( get_num_nln_eq(true) > 0 )
    demoOpt->register_nln_eq_fn(this);   // TPL_SPECIFIC
  if( get_num_nln_ineq(true) > 0 )
    demoOpt->register_nln_ineq_fn(this); // TPL_SPECIFIC
}


// -----------------------------------------------------------------

// core_run redefines the Optimizer virtual function to perform the
// optimization using Demo_Opt and catalogue the results.  core_run
// will be called by Dakota, and core_run will call the TPL optimizer.

void DemoTPLOptimizer::core_run()
{
  // Call a helper function to set the initial values of the
  // variables.  It is implemented later in this source file.

  initialize_variables_and_constraints();

  // Invoke the TPL's method to actually perform the optimization.
  // This code should be replaced by whatever the TPL's mechanism is
  // for running its solver.

  // ------------------  TPL_SPECIFIC  ------------------
  demoOpt->execute(true);

  // The TPL should provide the optimal value of the objective
  // function and the associated variable values.  For the purposes of
  // this example, the values should be returned to standard C++ data
  // types, double for the function value and std::vector<double> for
  // the variable values.

  if (!localObjectiveRecast) {
    // Replace this line with however the TPL being integrated returns
    // the optimal function value.  To use this demo with minimal
    // changes, the returned value needs to be (converted to) a
    // double.
    double best_f = demoOpt->get_best_f(); // TPL_SPECIFIC

    // If the TPL defaults to doing minimization, no need to do
    // anything with this code.  It manages needed sign changes
    // depending on whether minimize or maximize has been specified in
    // the Dakota input file.
    const BoolDeque& max_sense = iteratedModel.primary_response_fn_sense();
    RealVector best_fns(iteratedModel.response_size());

    // Get best (single) objcetive value respecting max/min expectations
    best_fns[0] = (!max_sense.empty() && max_sense[0]) ?  -best_f : best_f;

    // Get best Nonlinear Equality Constraints from TPL
    if( numNonlinearEqConstraints > 0 )
    {
      auto best_nln_eqs = demoOpt->get_best_nln_eqs();
      //std::copy( best_nln_eqs.begin(), best_nln_eqs.end(), &best_fns(0)+1);
      dataTransferHandler->get_best_nonlinear_eq_constraints_from_tpl(
                                          best_nln_eqs,
                                          best_fns);
    }

    // Get best Nonlinear Inequality Constraints from TPL
    if( numNonlinearIneqConstraints > 0 )
    {
      auto best_nln_ineqs = demoOpt->get_best_nln_ineqs(); // TPL_SPECIFIC

      dataTransferHandler->get_best_nonlinear_ineq_constraints_from_tpl(
                                          best_nln_ineqs,
                                          best_fns);
    }

    bestResponseArray.front().function_values(best_fns);
  }

  std::vector<double> best_x = demoOpt->get_best_x(); // TPL_SPECIFIC

  // Set Dakota optimal value data.
  set_variables<>(best_x, iteratedModel, bestVariablesArray.front());

} // core_run


// -----------------------------------------------------------------

// Dakota will call initialize_run() for any one-time setup.  If the
// TPL being integrated requires such, it should be implemented here.
// Replace the demoOpt method call with any initialization (or call to
// initialization)needed.

void DemoTPLOptimizer::initialize_run()
{
  Optimizer::initialize_run();
  demoOpt->initialize(true); // TPL_SPECIFIC
}


// -----------------------------------------------------------------

// This helper function sets the TPL algorithmic parameters.  For this
// initial example, the only mechanism to do this is for the TPL to
// directly read in a file with the information.  The formatting and
// reading of the file are the TPL's responsibility.  The file name is
// specified in the Dakota input file.  This code extracts the
// filename and passes it to the TPL.

void DemoTPLOptimizer::set_demo_parameters()
{
  int    max_fn_evals;
  int    max_iters;
  double conv_tol;
  double min_var_chg;
  double obj_target;

  // Values for common stopping criteria can be obtained from Dakota.
  // If user has provided values in the input file, those values will
  // be returned.  Otherwise, Dakota defaults will be returned.
  get_common_stopping_criteria(max_fn_evals, max_iters, conv_tol, min_var_chg, obj_target );

  // ------------------  TPL_SPECIFIC  ------------------
  demoOpt->set_param("Maximum Evaluations", max_fn_evals);
  demoOpt->set_param("Maximum Iterations",  max_iters);
  demoOpt->set_param("Function Tolerance",  obj_target);

  // Check for native Demo_Opt input file.  The file name needs to be
  // included in the Dakota input file.
  String adv_opts_file = probDescDB.get_string("method.advanced_options_file");
  if (!adv_opts_file.empty())
  {
    if (!boost::filesystem::exists(adv_opts_file))
    {
      Cerr << "\nError: Demo_Opt options_file '" << adv_opts_file
	   << "' specified, but file not found.\n";
      abort_handler(METHOD_ERROR);
    }
  }

  // Replace this line by whatever the TPL being integrated uses to
  // set its input file name.

  demoOpt->set_solver_options(adv_opts_file, true); // TPL_SPECIFIC

} // set_demo_parameters

// -----------------------------------------------------------------

// This helper function gets the initial values of the variables and
// the values of the bound constraints.  They are returned in standard
// C++ data types.  The values are passed on to the TPL.

void DemoTPLOptimizer::initialize_variables_and_constraints()
{

  // Get the number of variables, the initial values, and the values
  // of bound constraints.  They are returned to standard C++ data
  // types.  This example considers only continuous variables.  Other
  // types of variables and constraints will be added at a later time.
  // Note that double is aliased to Real in Dakota.
  int num_total_vars = numContinuousVars;
  std::vector<Real> init_point(num_total_vars);
  std::vector<Real> lower(num_total_vars),
                    upper(num_total_vars);

  // More on DemoOptTraits can be found in DemoOptimizer.hpp.
  get_variables(iteratedModel, init_point);
  get_variable_bounds_from_dakota<DemoOptTraits>( lower, upper );

  // Replace this line by whatever the TPL being integrated uses to
  // ingest variable values and bounds, including any data type
  // conversion needed.

  // ------------------  TPL_SPECIFIC  ------------------
  demoOpt->set_problem_data(init_point,   //  "Initial Guess"
                            lower     ,   //  "Lower Bounds"
                            upper      ); //  "Upper Bounds"

} // initialize_variables_and_constraints

// -----------------------------------------------------------------

// This is the implementation of the objective function evaluation.
// This assumes a function callback approach, i.e., the TPL optimizer
// calls this function whenever it needs an evaluation done.  Other
// ways to interface to function will be added in the future.  This
// interface should be replaced with what ever interface the TPL uses.

Real
DemoTPLOptimizer::compute_obj(const std::vector<double> & x, bool verbose)
{
  // Tell Dakota what variable values to use for the function
  // valuation.  x must be (converted to) a std::vector<double> to use
  // this demo with minimal changes.
  set_variables<>(x, iteratedModel, iteratedModel.current_variables());

  // Evaluate the function at the specified x.
  iteratedModel.evaluate();

  // Retrieve the the function value and sign it appropriately based
  // on whether minimize or maximize has been specified in the Dakota
  // input file.
  double f = dataTransferHandler->get_response_value_from_dakota(iteratedModel.current_response());

  return f;
}


// -----------------------------------------------------------------

// This is the implementation of the nonlinear equality constraint evaluation.
// This assumes a function callback approach, i.e., the TPL optimizer
// calls this function whenever it needs an evaluation done.  Other
// ways to interface to function will be added in the future.  This
// interface should be replaced with what ever interface the TPL uses.

int
DemoTPLOptimizer::get_num_nln_eq(bool verbose)
{
  return dataTransferHandler->num_dakota_nonlin_eq_constraints();
}


// -----------------------------------------------------------------

void
DemoTPLOptimizer::compute_nln_eq(std::vector<Real> &c, const std::vector<Real> &x, bool verbose)
{
  // Tell Dakota what variable values to use for the nonlinear constraint
  // evaluations.  x must be (converted to) a std::vector<double> to use
  // this demo with minimal changes.
  set_variables<>(x, iteratedModel, iteratedModel.current_variables());

  // Evaluate the function at the specified x.
  iteratedModel.evaluate();

  // Use an adapter to copy data
  dataTransferHandler->get_nonlinear_eq_constraints_from_dakota(iteratedModel.current_response(), c);

} // nonlinear eq constraints value


// -----------------------------------------------------------------

// This is the implementation of the nonlinear inequality constraint evaluation.
// This assumes a function callback approach, i.e., the TPL optimizer
// calls this function whenever it needs an evaluation done.  Other
// ways to interface to function will be added in the future.  This
// interface should be replaced with what ever interface the TPL uses.

int
DemoTPLOptimizer::get_num_nln_ineq(bool verbose)
{
  return dataTransferHandler->num_tpl_nonlin_ineq_constraints();
}

// -----------------------------------------------------------------

void
DemoTPLOptimizer::compute_nln_ineq(std::vector<Real> &c, const std::vector<Real> &x, bool verbose)
{
  // Tell Dakota what variable values to use for the nonlinear constraint
  // evaluations.  x must be (converted to) a std::vector<double> to use
  // this demo with minimal changes.
  set_variables<>(x, iteratedModel, iteratedModel.current_variables());

  // Evaluate the function at the specified x.
  iteratedModel.evaluate();

  // Use an adapter to copy data from Dakota into Demo_Opt
  dataTransferHandler->get_nonlinear_ineq_constraints_from_dakota(iteratedModel.current_response(), c);

} // nonlinear ineq constraints value

} // namespace Dakota
