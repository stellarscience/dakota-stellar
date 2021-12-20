/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2020
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#include "PEBBLMinimizer.hpp"
#include "PRPMultiIndex.hpp"
#include "dakota_data_io.hpp"

namespace Dakota 
{
extern PRPCache data_pairs; // global container

PebbldMinimizer::PebbldMinimizer(ProblemDescDB& problem_db, Model& model):
  Minimizer(problem_db, model, std::shared_ptr<TraitsBase>(new PebbldTraits()))
{
  // While this copy will be replaced in best update, initialize here
  // since relied on in Minimizer::initialize_run when a sub-iterator
  bestVariablesArray.push_back(
    iteratedModel.current_variables().copy());

  // Instantiate the approximate sub-problem minimizer
  const String& subprob_method_ptr
    = probDescDB.get_string("method.sub_method_pointer");
  const String& subprob_method_name
    = probDescDB.get_string("method.sub_method_name");

  if (!subprob_method_ptr.empty()) {
    // Approach 1: method spec support for subProbMinimizer
    const String& model_ptr = probDescDB.get_string("method.model_pointer");
    size_t method_index = probDescDB.get_db_method_node(); // for restoration
    probDescDB.set_db_method_node(subprob_method_ptr); // method only
    // sub-problem minimizer will use shallow copy of iteratedModel
    // (from problem_db.get_model())
    subProbMinimizer = probDescDB.get_iterator();//(iteratedModel);
    // suppress DB ctor default and don't output summary info
    subProbMinimizer.summary_output(false);
    // verify method's modelPointer is empty or consistent
    const String& am_model_ptr = probDescDB.get_string("method.model_pointer");
    if (!am_model_ptr.empty() && am_model_ptr != model_ptr)
      Cerr << "Warning: BandB method_pointer specification includes an\n"
	   << "         inconsistent model_pointer that will be ignored."
	   << std::endl;
    probDescDB.set_db_method_node(method_index); // restore method only
  }
  else if (!subprob_method_name.empty())
    // Approach 2: instantiate on-the-fly w/o method spec support
    subProbMinimizer
      = probDescDB.get_iterator(subprob_method_name, iteratedModel);

  branchAndBound = new PebbldBranching();
  branchAndBound->setModel(model);
  branchAndBound->setIterator(subProbMinimizer);
}

PebbldMinimizer::PebbldMinimizer(Model &model)
	: Minimizer(BRANCH_AND_BOUND, model, std::shared_ptr<TraitsBase>(new PebbldTraits())) 
{//branchAndBound(model)
};

PebbldMinimizer::PebbldMinimizer(Model &model, int random_seed, size_t max_iter, size_t max_eval) 
       : Minimizer(BRANCH_AND_BOUND, model, std::shared_ptr<TraitsBase>(new PebbldTraits()))
{//branchAndBound(model,random_seed, max_iter, max_eval)
};

/** redefines default iterator logic to execute branch and bound and
    extract optimization results. **/
void PebbldMinimizer::core_run()
{
	InitializeTiming();
	branchAndBound->search();

	pebbl::arraySolution<double>* newSolution = dynamic_cast<pebbl::arraySolution<double>*>(branchAndBound->getSolution());
	  
	RealVector variables(numContinuousVars);
	RealVector best_fns(1);
	for (size_t i=0; i<numContinuousVars; i++)
	  variables[i] = newSolution->array[i];
	bestVariablesArray.front().continuous_variables(variables);
	best_fns[0] = newSolution->value;
	bestResponseArray.front().function_values(best_fns);
}

/** Redefines default iterator results printing to include
    optimization results (objective functions and constraints). */
void PebbldMinimizer::print_results(std::ostream& s, short results_state)
{
  size_t i, num_best = bestVariablesArray.size();
  if (num_best != bestResponseArray.size()) {
    Cerr << "\nError: mismatch in lengths of bestVariables and bestResponses."
         << std::endl;
    abort_handler(-1); 
  } 

  const String& interface_id = iteratedModel.interface_id();
  int eval_id;
  activeSet.request_values(1);

  // -------------------------------------
  // Single and Multipoint results summary
  // -------------------------------------
  for (i=0; i<num_best; ++i) {
    // output best variables
    s << "<<<<< Best parameters          ";
    if (num_best > 1) s << "(set " << i+1 << ") ";
    s << "=\n" << bestVariablesArray[i];
    // output best response
    const RealVector& best_fns = bestResponseArray[i].function_values();
    if (optimizationFlag) {
      if (numUserPrimaryFns > 1) s << "<<<<< Best objective functions ";
      else                       s << "<<<<< Best objective function  ";
      if (num_best > 1) s << "(set " << i+1 << ") "; s << "=\n";
      write_data_partial(s, (size_t)0, numUserPrimaryFns, best_fns);
    }
    else 
      print_residuals(numUserPrimaryFns, best_fns,
                              RealVector(),
                              num_best, i,
                              s);
    size_t num_cons = numFunctions - numUserPrimaryFns;
    //        if (num_cons) {
    //          s << "<<<<< Best constraint values   ";
    //          if (num_best > 1) s << "(set " << i+1 << ") "; s << "=\n";
    //         write_data_partial(s, numUserPrimaryFns, num_cons, best_fns);
    //        }
    // lookup evaluation id where best occurred.  This cannot be catalogued 
    // directly because the optimizers track the best iterate internally and 
    // return the best results after iteration completion.  Therfore, perform a
    // search in data_pairs to extract the evalId for the best fn eval.
    PRPCacheHIter cache_it = lookup_by_val(data_pairs, interface_id,
					   bestVariablesArray[i], activeSet);
    if (cache_it == data_pairs.get<hashed>().end()) {
      s << "<<<<< Best data not found in evaluation cache\n\n";
      eval_id = 0;
    }
    else {
      eval_id = cache_it->eval_id();
      if (eval_id > 0)
	s << "<<<<< Best data captured at function evaluation " << eval_id
	  << "\n\n";
      else // should not occur
	s << "<<<<< Best data not found in evaluations from current execution,"
	  << "\n      but retrieved from restart archive with evaluation id "
	  << -eval_id << "\n\n";
    }
  }
}

} // namespace Dakota
