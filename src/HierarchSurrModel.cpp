/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       HierarchSurrModel
//- Description: Implementation code for the HierarchSurrModel class
//- Owner:       Mike Eldred
//- Checked by:

#include "HierarchSurrModel.hpp"
#include "ProblemDescDB.hpp"

static const char rcsId[]=
  "@(#) $Id: HierarchSurrModel.cpp 6656 2010-02-26 05:20:48Z mseldre $";

namespace Dakota {


HierarchSurrModel::HierarchSurrModel(ProblemDescDB& problem_db):
  SurrogateModel(problem_db),
  corrOrder(problem_db.get_short("model.surrogate.correction_order")),
  correctionMode(SINGLE_CORRECTION)//, componentParallelKey(_NPOS,_NPOS)
{
  // Hierarchical surrogate models pass through numerical derivatives
  supportsEstimDerivs = false;
  // initialize ignoreBounds even though it's irrelevant for pass through
  ignoreBounds = problem_db.get_bool("responses.ignore_bounds");
  // initialize centralHess even though it's irrelevant for pass through
  centralHess = problem_db.get_bool("responses.central_hess");

  const StringArray& ordered_model_ptrs
    = problem_db.get_sa("model.surrogate.ordered_model_pointers");

  size_t i, num_models = ordered_model_ptrs.size(),
           model_index = problem_db.get_db_model_node(); // for restoration

  const std::pair<short,short>& cv_view = currentVariables.view();
  orderedModels.resize(num_models);
  for (i=0; i<num_models; ++i) {
    problem_db.set_db_model_nodes(ordered_model_ptrs[i]);
    orderedModels[i] = problem_db.get_model();
    check_submodel_compatibility(orderedModels[i]);
    //if (cv_view != orderedModels[i].current_variables().view()) {
    //  Cerr << "Error: variable views in hierarchical models must be "
    //       << "identical." << std::endl;
    //  abort_handler(-1);
    //}
  }

  problem_db.set_db_model_nodes(model_index); // restore

  // default index values, to be overridden at run time
  if (num_models == 1) { // first and last solution level (1 model)
    surrModelKey.resize(2); truthModelKey.resize(2);
    surrModelKey[0] = surrModelKey[1] = truthModelKey[0] = 0;
    truthModelKey[1] = orderedModels[0].solution_levels() - 1;
  }
  else { // first and last model form (solution levels ignored)
    surrModelKey.resize(1); truthModelKey.resize(1);
    surrModelKey[0] = 0; truthModelKey[0] = num_models - 1;
    //surrModelKey[1] = truthModelKey[1] = _NPOS;
  }
  check_model_interface_instance();

  // Correction is required in HierarchSurrModel for some responseModes.
  // Enforcement of a correction type for these modes occurs in
  // surrogate_response_mode(short).
  if (corrType) // initialize DiscrepancyCorrection using initial LF indices
    deltaCorr[fidelity_keys()].initialize(surrogate_model(), surrogateFnIndices,
					  corrType, corrOrder);

  //truthResponseRef[truthModelKey] = currentResponse.copy();
}


void HierarchSurrModel::check_submodel_compatibility(const Model& sub_model)
{
  SurrogateModel::check_submodel_compatibility(sub_model);
  
  bool error_flag = false;
  // Check for compatible array sizing between sub_model and currentResponse.
  // HierarchSurrModel creates aggregations and DataFitSurrModel consumes them.
  // For now, allow either a factor of 2 or 1 from aggregation or not.  In the
  // future, aggregations may span a broader model hierarchy (e.g., factor =
  // orderedModels.size()).  In general, the fn count check needs to be
  // specialized in the derived classes.
  size_t sm_qoi = sub_model.qoi();//, aggregation = numFns / sm_qoi;
  if ( numFns % sm_qoi ) { //|| aggregation < 1 || aggregation > 2 ) {
    Cerr << "Error: incompatibility between approximate and actual model "
	 << "response function sets\n       within HierarchSurrModel: "<< numFns
	 << " approximate and " << sm_qoi << " actual functions.\n       "
	 << "Check consistency of responses specifications." << std::endl;
    error_flag = true;
  }

  // TO DO: Bayes exp design (hi2lo) introduces new requirements on a
  // hierarchical model, and MF active subspaces will as well.
  // > For (simulation-based) OED, one option is to enforce consistency in
  //   inactive state (config vars) and allow active parameterization to vary.
  // > For hi2lo, this implies that the active variable subset could be null
  //   for HF, as the active calibration variables only exist for LF.
  size_t sm_icv = sub_model.icv(),  sm_idiv = sub_model.idiv(),
    sm_idsv = sub_model.idsv(),     sm_idrv = sub_model.idrv(),
    icv  = currentVariables.icv(),  idiv = currentVariables.idiv(),
    idsv = currentVariables.idsv(), idrv = currentVariables.idrv();
  if (sm_icv != icv || sm_idiv != idiv || sm_idsv != idsv || sm_idrv != idrv) {
    Cerr << "Error: incompatibility between approximate and actual model "
	 << "variable sets within\n       HierarchSurrModel: inactive "
	 << "approximate = " << icv << " continuous, " << idiv
	 << " discrete int, " << idsv << " discrete string, and " << idrv
	 << " discrete real and\n       inactive actual = " << sm_icv
	 << " continuous, " << sm_idiv << " discrete int, " << sm_idsv
	 << " discrete string, and " << sm_idrv << " discrete real.  Check "
	 << "consistency of variables specifications." << std::endl;
    error_flag = true;
  }

  if (error_flag)
    abort_handler(-1);
}


void HierarchSurrModel::
derived_init_communicators(ParLevLIter pl_iter, int max_eval_concurrency,
                           bool recurse_flag)
{
  // responseMode is a run-time setting (in SBLMinimizer, it is switched among
  // AUTO_CORRECTED_SURROGATE, BYPASS_SURROGATE, and UNCORRECTED_SURROGATE;
  // in NonDExpansion, it is switching between MODEL_DISCREPANCY and
  // UNCORRECTED_SURROGATE).  Since it is neither static nor generally
  // available at construct/init time, take a conservative approach with init
  // and free and a more aggressive approach with set.

  if (recurse_flag) {
    size_t i, model_index = probDescDB.get_db_model_node(), // for restoration
              num_models = orderedModels.size();

    // init and free must cover possible subset of active responseModes and
    // ordered model fidelities, but only 2 models at mpst will be active at
    // runtime.  In order to reduce the number of parallel configurations, we
    // group the responseModes into two sets: (1) the correction-based set
    // commonly used in surrogate-based optimization et al., and (2) the
    // aggregation-based set commonly used in multilevel/multifidelity UQ.

    // TO DO: would like a better detection option, but passing the mode from
    // the Iterator does not work in parallel w/o an additional bcast (Iterator
    // only instantiated on iteratorComm rank 0).  For now, we will infer it
    // from an associated method spec at init time.
    // Note: responseMode gets bcast at run time in component_parallel_mode()
    bool extra_deriv_config
      = (probDescDB.get_ushort("method.algorithm") & MINIMIZER_BIT);
    //(responseMode == UNCORRECTED_SURROGATE ||
    // responseMode == BYPASS_SURROGATE ||
    // responseMode == AUTO_CORRECTED_SURROGATE);

    for (i=0; i<num_models; ++i) {
      Model& model_i = orderedModels[i];
      // superset of possible init calls (two configurations for i > 0)
      probDescDB.set_db_model_nodes(model_i.model_id());
      model_i.init_communicators(pl_iter, max_eval_concurrency);
      if (extra_deriv_config) // && i) // mid and high fidelity only?
        model_i.init_communicators(pl_iter, model_i.derivative_concurrency());
    }


    /* This version inits only two models
    Model& lf_model = surrogate_model();
    Model& hf_model = truth_model();

    // superset of possible init calls (two configurations for HF)
    probDescDB.set_db_model_nodes(lf_model.model_id());
    lf_model.init_communicators(pl_iter, max_eval_concurrency);

    probDescDB.set_db_model_nodes(hf_model.model_id());
    hf_model.init_communicators(pl_iter, hf_model.derivative_concurrency());
    hf_model.init_communicators(pl_iter, max_eval_concurrency);
    */


    /* This version does not support runtime updating of responseMode
    switch (responseMode) {
    case UNCORRECTED_SURROGATE:
      // LF are used in iterator evals
      lf_model.init_communicators(pl_iter, max_eval_concurrency);
      break;
    case AUTO_CORRECTED_SURROGATE:
      // LF are used in iterator evals
      lf_model.init_communicators(pl_iter, max_eval_concurrency);
      // HF evals are for correction and validation:
      // concurrency = one eval at a time * derivative concurrency per eval
      hf_model.init_communicators(pl_iter, hf_model.derivative_concurrency());
      break;
    case BYPASS_SURROGATE:
      // HF are used in iterator evals
      hf_model.init_communicators(pl_iter, max_eval_concurrency);
      break;
    case MODEL_DISCREPANCY: case AGGREGATED_MODELS:
      // LF and HF are used in iterator evals
      lf_model.init_communicators(pl_iter, max_eval_concurrency);
      hf_model.init_communicators(pl_iter, max_eval_concurrency);
      break;
    }
    */

    probDescDB.set_db_model_nodes(model_index); // restore all model nodes
  }
}


void HierarchSurrModel::
derived_set_communicators(ParLevLIter pl_iter, int max_eval_concurrency,
                          bool recurse_flag)
{
  miPLIndex = modelPCIter->mi_parallel_level_index(pl_iter);// run time setting

  // HierarchSurrModels do not utilize default set_ie_asynchronous_mode() as
  // they do not define the ie_parallel_level

  // This aggressive logic is appropriate for invocations of the Model via
  // Iterator::run(), but is fragile w.r.t. invocations of the Model outside
  // this scope (e.g., Model::evaluate() within SBLMinimizer).  The default
  // responseMode value is {AUTO_,UN}CORRECTED_SURROGATE, which mitigates
  // the specific case of SBLMinimizer, but the general fragility remains.
  if (recurse_flag) {

    // bcast not needed for recurse_flag=false in serve_run call to set_comms
    //if (pl_iter->server_communicator_size() > 1)
    //  parallelLib.bcast(responseMode, *pl_iter);

    switch (responseMode) {

    // CASES WITH A SINGLE ACTIVE MODEL:

    case UNCORRECTED_SURROGATE: {
      Model& lf_model = surrogate_model();
      lf_model.set_communicators(pl_iter, max_eval_concurrency);
      asynchEvalFlag     = lf_model.asynch_flag();
      evaluationCapacity = lf_model.evaluation_capacity();
      break;
    }
    case BYPASS_SURROGATE: {
      Model& hf_model = truth_model();
      hf_model.set_communicators(pl_iter, max_eval_concurrency);
      asynchEvalFlag     = hf_model.asynch_flag();
      evaluationCapacity = hf_model.evaluation_capacity();
      break;
    }

    // CASES WHERE ANY/ALL MODELS COULD BE ACTIVE:

    case AUTO_CORRECTED_SURROGATE: {
      // Lowest fidelity model is interfaced with minimizer:
      Model& model_0 = orderedModels[0];
      model_0.set_communicators(pl_iter, max_eval_concurrency);
      asynchEvalFlag     = model_0.asynch_flag();
      evaluationCapacity = model_0.evaluation_capacity();

      // TO DO: this will not be true for multigrid optimization:
      bool use_deriv_conc = true; // only verifications/corrections
      // Either need detection logic, a passed option, or to abandon the
      // specialization and just generalize init/set/free to use the max
      // of the two values...

      // Loop over all higher fidelity models:
      size_t i, num_models = orderedModels.size(); int cap_i;
      for (i=1; i<num_models; ++i) {
	Model& model_i = orderedModels[i];
	if (use_deriv_conc) {
	  int deriv_conc_i = model_i.derivative_concurrency();
	  model_i.set_communicators(pl_iter, deriv_conc_i);
	  if (deriv_conc_i > 1 && model_i.asynch_flag()) asynchEvalFlag = true;
	}
	else {
	  model_i.set_communicators(pl_iter, max_eval_concurrency);
	  if (model_i.asynch_flag()) asynchEvalFlag = true;
	}
	cap_i = model_i.evaluation_capacity();
	if (cap_i > evaluationCapacity) evaluationCapacity = cap_i;
      }
      break;
    }
    case MODEL_DISCREPANCY: case AGGREGATED_MODELS: {
      size_t i, num_models = orderedModels.size(); int cap_i;
      asynchEvalFlag = false; evaluationCapacity = 1;
      for (i=0; i<num_models; ++i) {
	Model& model_i = orderedModels[i];
	model_i.set_communicators(pl_iter, max_eval_concurrency);
	if (model_i.asynch_flag()) asynchEvalFlag = true;
	cap_i = model_i.evaluation_capacity();
	if (cap_i > evaluationCapacity) evaluationCapacity = cap_i;
      }
      break;
    }
    }
  }
}


void HierarchSurrModel::
derived_free_communicators(ParLevLIter pl_iter, int max_eval_concurrency,
                           bool recurse_flag)
{
  if (recurse_flag) {

    size_t i, num_models = orderedModels.size();
    bool extra_deriv_config = true;//(responseMode == UNCORRECTED_SURROGATE ||
    // responseMode == BYPASS_SURROGATE ||
    // responseMode == AUTO_CORRECTED_SURROGATE);
    for (i=0; i<num_models; ++i) {
      Model& model_i = orderedModels[i];
      // superset of possible init calls (two configurations for i > 0)
      model_i.free_communicators(pl_iter, max_eval_concurrency);
      if (extra_deriv_config) // && i) // mid and high fidelity only?
        model_i.free_communicators(pl_iter, model_i.derivative_concurrency());
    }


    /* This version frees only two models:
    // superset of possible free calls (two configurations for HF)
    surrogate_model().free_communicators(pl_iter, max_eval_concurrency);
    Model& hf_model = truth_model();
    hf_model.free_communicators(pl_iter, hf_model.derivative_concurrency());
    hf_model.free_communicators(pl_iter, max_eval_concurrency);
    */


    /* This version does not support runtime updating of responseMode:
    switch (responseMode) {
    case UNCORRECTED_SURROGATE:
      lf_model.free_communicators(pl_iter, max_eval_concurrency);
      break;
    case AUTO_CORRECTED_SURROGATE:
      lf_model.free_communicators(pl_iter, max_eval_concurrency);
      hf_model.free_communicators(pl_iter, hf_model.derivative_concurrency());
      break;
    case BYPASS_SURROGATE:
      hf_model.free_communicators(pl_iter, max_eval_concurrency);
      break;
    case MODEL_DISCREPANCY: case AGGREGATED_MODELS:
      lf_model.free_communicators(pl_iter, max_eval_concurrency);
      hf_model.free_communicators(pl_iter, max_eval_concurrency);
      break;
    }
    */
  }
}


/** Inactive variables must be propagated when a HierarchSurrModel
    is employed by a sub-iterator (e.g., OUU with MLMC or MLPCE).
    In current use cases, this can occur once per sub-iterator
    execution within Model::initialize_mapping(). */
bool HierarchSurrModel::initialize_mapping(ParLevLIter pl_iter)
{
  Model::initialize_mapping(pl_iter);

  // push inactive variable values/bounds from currentVariables and
  // userDefinedConstraints into orderedModels
  size_t i, num_models = orderedModels.size();
  for (i=0; i<num_models; ++i) {
    orderedModels[i].initialize_mapping(pl_iter);
    init_model(orderedModels[i]);
  }

  return false; // no change to problem size
}


/** Inactive variables must be propagated when a HierarchSurrModel
    is employed by a sub-iterator (e.g., OUU with MLMC or MLPCE).
    In current use cases, this can occur once per sub-iterator
    execution within Model::initialize_mapping(). */
bool HierarchSurrModel::finalize_mapping()
{
  size_t i, num_models = orderedModels.size();
  for (i=0; i<num_models; ++i)
    orderedModels[i].finalize_mapping();

  Model::finalize_mapping();

  return false; // no change to problem size
}


void HierarchSurrModel::build_approximation()
{
  Cout << "\n>>>>> Building hierarchical approximation.\n";

  // perform the eval for the low fidelity model
  // NOTE: For SBO, the low fidelity eval is performed externally and its
  // response is passed into compute_correction.
  // -->> move LF model out and restructure if(!approxBuilds)
  //ActiveSet temp_set = lf_model.current_response().active_set();
  //temp_set.request_values(1);
  //if (sameModelInstance)
  //  lf_model.solution_level_index(surrogate_level_index());
  //lf_model.evaluate(temp_set);
  //const Response& lo_fi_response = lf_model.current_response();

  Model& hf_model = truth_model();
  if (hierarchicalTagging) {
    String eval_tag = evalTagPrefix + '.' +
                      boost::lexical_cast<String>(surrModelEvalCntr+1);
    hf_model.eval_tag_prefix(eval_tag);
  }

  // set HierarchSurrModel parallelism mode to HF model
  component_parallel_mode(TRUTH_MODEL);

  // update HF model with current variable values/bounds/labels
  update_model(hf_model);

  // store inactive variable values for use in determining whether an
  // automatic rebuild of an approximation is required
  // (reference{C,D}{L,U}Bnds are not needed in the hierarchical case)
  const Variables& hf_vars = hf_model.current_variables();
  copy_data(hf_vars.inactive_continuous_variables(),    referenceICVars);
  copy_data(hf_vars.inactive_discrete_int_variables(),  referenceIDIVars);
  referenceIDSVars = hf_vars.inactive_discrete_string_variables();
  copy_data(hf_vars.inactive_discrete_real_variables(), referenceIDRVars);

  // compute the response for the high fidelity model
  ShortArray total_asv, hf_asv, lf_asv;
  DiscrepCorrMap::iterator dc_it = deltaCorr.find(fidelity_keys());
  if (dc_it!=deltaCorr.end() && dc_it->second.initialized())
    total_asv.assign(numFns, dc_it->second.data_order());
  else
    total_asv.assign(numFns, 1); // default: values only if no deriv correction
  asv_split(total_asv, hf_asv, lf_asv, true);

  if ( truthResponseRef.find(truthModelKey) == truthResponseRef.end() )
    truthResponseRef[truthModelKey] = currentResponse.copy();

  ActiveSet hf_set = currentResponse.active_set(); // copy
  hf_set.request_vector(hf_asv);
  if (sameModelInstance)
    hf_model.solution_level_index(truth_level_index());
  hf_model.evaluate(hf_set);
  truthResponseRef[truthModelKey].update(hf_model.current_response());

  // could compute the correction to LF model here, but rely on an external
  // call for consistency with DataFitSurr and to facilitate SBO logic.  In
  // particular, lo_fi_response involves find_center(), hard conv check, etc.
  //deltaCorr[indices].compute(..., truthResponseRef, lo_fi_response);

  Cout << "\n<<<<< Hierarchical approximation build completed.\n";
  ++approxBuilds;
}


/*
bool HierarchSurrModel::
build_approximation(const RealVector& c_vars, const Response& response)
{
  // NOTE: this fn not currently used by SBO, but it could be.

  // Verify data content in incoming response
  const ShortArray& asrv = response.active_set_request_vector();
  bool data_complete = true; short corr_order = dataCorr.correction_order();
  for (size_t i=0; i<numFns; i++)
    if ( ( corr_order == 2 &&  (asrv[i] & 7) != 7 ) ||
	 ( corr_order == 1 &&  (asrv[i] & 3) != 3 ) ||
	 ( corr_order == 0 && !(asrv[i] & 1) ) )
      data_complete = false;
  if (data_complete) {
    Cout << "\n>>>>> Updating hierarchical approximation.\n";

    // are these updates necessary?
    Model& hf_model = truth_model();
    currentVariables.continuous_variables(c_vars);
    update_model(hf_model);
    const Variables& hf_vars = hf_model.current_variables();
    copy_data(hf_vars.inactive_continuous_variables(), referenceICVars);
    copy_data(hf_vars.inactive_discrete_variables(),   referenceIDVars);

    truthResponseRef.update(response);

    Cout << "\n<<<<< Hierarchical approximation update completed.\n";
  }
  else {
    Cerr << "Warning: cannot use anchor point in HierarchSurrModel::"
	 << "build_approximation(RealVector&, Response&).\n";
    currentVariables.continuous_variables(c_vars);
    build_approximation();
  }
  return false; // correction is not embedded and must be computed (by SBO)
}
*/


/** Compute the response synchronously using LF model, HF model, or
    both (mixed case).  For the LF model portion, compute the high
    fidelity response if needed with build_approximation(), and, if
    correction is active, correct the low fidelity results. */
void HierarchSurrModel::derived_evaluate(const ActiveSet& set)
{
  ++surrModelEvalCntr;

  // define LF/HF evaluation requirements
  ShortArray hi_fi_asv, lo_fi_asv;
  bool hi_fi_eval, lo_fi_eval, mixed_eval;
  Response lo_fi_response, hi_fi_response; // don't use truthResponseRef
  switch (responseMode) {
  case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
  case AGGREGATED_MODELS:
    asv_split(set.request_vector(), hi_fi_asv, lo_fi_asv, false);
    hi_fi_eval = !hi_fi_asv.empty(); lo_fi_eval = !lo_fi_asv.empty();
    mixed_eval = (hi_fi_eval && lo_fi_eval);            break;
  case BYPASS_SURROGATE:
    hi_fi_eval = true; lo_fi_eval = mixed_eval = false; break;
  case MODEL_DISCREPANCY:
    hi_fi_eval = lo_fi_eval = mixed_eval = true;        break;
  }

  Model& lf_model = surrogate_model();
  Model& hf_model = truth_model();
  if (hierarchicalTagging) {
    String eval_tag = evalTagPrefix + '.' +
                      boost::lexical_cast<String>(surrModelEvalCntr+1);
    if (sameModelInstance) lf_model.eval_tag_prefix(eval_tag);
    else {
      if (hi_fi_eval) hf_model.eval_tag_prefix(eval_tag);
      if (lo_fi_eval) lf_model.eval_tag_prefix(eval_tag);
    }
  }

  if (sameModelInstance) update_model(lf_model);

  // Notes on repetitive setting of model.solution_level_index():
  // > when LF & HF are the same model, then setting the index for low or high
  //   invalidates the other fidelity definition.
  // > within a single derived_evaluate(), could protect these updates with
  //   "if (sameModelInstance && mixed_eval)", but this does not guard against
  //   changes in eval requirements from the previous evaluation.  Detecting
  //   the current solution index state is currently as expensive as resetting
  //   it, so just reset each time.

  // ------------------------------
  // Compute high fidelity response
  // ------------------------------
  if (hi_fi_eval) {
    component_parallel_mode(TRUTH_MODEL); // TO DO: sameModelInstance
    if (sameModelInstance)
      hf_model.solution_level_index(truth_level_index());
    else
      update_model(hf_model);
    switch (responseMode) {
    case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
    case AGGREGATED_MODELS: {
      ActiveSet hi_fi_set(hi_fi_asv, set.derivative_vector());
      hf_model.evaluate(hi_fi_set);
      if (mixed_eval)
        hi_fi_response = (sameModelInstance) ? // deep copy or shared rep
	  hf_model.current_response().copy() : hf_model.current_response();
      else {
        currentResponse.active_set(hi_fi_set);
        currentResponse.update(hf_model.current_response());
      }
      break;
    }
    case BYPASS_SURROGATE:
      hf_model.evaluate(set);
      currentResponse.active_set(set);
      currentResponse.update(hf_model.current_response());
      break;
    case MODEL_DISCREPANCY:
      hf_model.evaluate(set);
      hi_fi_response = (sameModelInstance) ? hf_model.current_response().copy()
                       : hf_model.current_response(); // shared rep
      break;
    }
  }

  // -----------------------------
  // Compute low fidelity response
  // -----------------------------
  if (lo_fi_eval) {
    // pre-process
    switch (responseMode) {
    case AUTO_CORRECTED_SURROGATE:
      // if build_approximation has not yet been called, call it now
      if (!approxBuilds || force_rebuild())
        build_approximation();
      break;
    }
    // compute the LF response
    component_parallel_mode(SURROGATE_MODEL); // TO DO: sameModelInstance
    if (sameModelInstance)
      lf_model.solution_level_index(surrogate_level_index());
    else
      update_model(lf_model);
    ActiveSet lo_fi_set;
    switch (responseMode) {
    case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
    case AGGREGATED_MODELS:
      lo_fi_set.request_vector(lo_fi_asv);
      lo_fi_set.derivative_vector(set.derivative_vector());
      lf_model.evaluate(lo_fi_set);
      break;
    case MODEL_DISCREPANCY:
      lf_model.evaluate(set);
      break;
    }

    // post-process
    switch (responseMode) {
    case AUTO_CORRECTED_SURROGATE: {
      // LF resp should not be corrected directly (see derived_synchronize())
      lo_fi_response = lf_model.current_response().copy();
      recursive_apply(currentVariables, lo_fi_response);
      if (!mixed_eval) {
        currentResponse.active_set(lo_fi_set);
        currentResponse.update(lo_fi_response);
      }
      break;
    }
    case UNCORRECTED_SURROGATE:
      if (mixed_eval)
        lo_fi_response = lf_model.current_response(); // shared rep
      else {
        currentResponse.active_set(lo_fi_set);
        currentResponse.update(lf_model.current_response());
      }
      break;
    }
  }

  // ------------------------------
  // perform any LF/HF aggregations
  // ------------------------------
  switch (responseMode) {
  case MODEL_DISCREPANCY: {
    // don't update surrogate data within deltaCorr[indices]'s Approximations;
    // just update currentResponse (managed as surrogate data at a higher level)
    bool quiet_flag = (outputLevel < NORMAL_OUTPUT);
    currentResponse.active_set(set);
    deltaCorr[fidelity_keys()].compute(hi_fi_response,
				       lf_model.current_response(),
				       currentResponse, quiet_flag);
    break;
  }
  case AGGREGATED_MODELS:
    aggregate_response(hi_fi_response, lf_model.current_response(),
                       currentResponse);
    break;
  case UNCORRECTED_SURROGATE:   case AUTO_CORRECTED_SURROGATE:
    if (mixed_eval) {
      currentResponse.active_set(set);
      response_combine(hi_fi_response, lo_fi_response, currentResponse);
    }
    break;
  }
}


/** Compute the response asynchronously using LF model, HF model, or
    both (mixed case).  For the LF model portion, compute the high
    fidelity response with build_approximation() (for correcting the
    low fidelity results in derived_synchronize() and
    derived_synchronize_nowait()) if not performed previously. */
void HierarchSurrModel::derived_evaluate_nowait(const ActiveSet& set)
{
  ++surrModelEvalCntr;

  Model& lf_model = surrogate_model();
  Model& hf_model = truth_model();

  ShortArray hi_fi_asv, lo_fi_asv;
  bool hi_fi_eval, lo_fi_eval, asynch_lo_fi = lf_model.asynch_flag(),
                               asynch_hi_fi = hf_model.asynch_flag();
  switch (responseMode) {
  case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
  case AGGREGATED_MODELS:
    asv_split(set.request_vector(), hi_fi_asv, lo_fi_asv, false);
    hi_fi_eval = !hi_fi_asv.empty();  lo_fi_eval = !lo_fi_asv.empty();  break;
  case BYPASS_SURROGATE:
    hi_fi_eval = true; lo_fi_eval = false;                              break;
  case MODEL_DISCREPANCY:
    hi_fi_eval = lo_fi_eval = true;                                     break;
  }

  if (hierarchicalTagging) {
    String eval_tag = evalTagPrefix + '.' +
                      boost::lexical_cast<String>(surrModelEvalCntr+1);
    if (sameModelInstance) lf_model.eval_tag_prefix(eval_tag);
    else {
      if (hi_fi_eval) hf_model.eval_tag_prefix(eval_tag);
      if (lo_fi_eval) lf_model.eval_tag_prefix(eval_tag);
    }
  }

  if (sameModelInstance) update_model(lf_model);

  // perform Model updates and define active sets for LF and HF evaluations
  ActiveSet hi_fi_set, lo_fi_set;
  if (hi_fi_eval) {
    // update HF model
    if (!sameModelInstance) update_model(hf_model);
    // update hi_fi_set
    hi_fi_set.derivative_vector(set.derivative_vector());
    switch (responseMode) {
    case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
    case AGGREGATED_MODELS:
      hi_fi_set.request_vector(hi_fi_asv);             break;
    case BYPASS_SURROGATE: case MODEL_DISCREPANCY:
      hi_fi_set.request_vector(set.request_vector());  break;
    }
  }
  if (lo_fi_eval) {
    // if build_approximation has not yet been called, call it now
    if ( responseMode == AUTO_CORRECTED_SURROGATE &&
         ( !approxBuilds || force_rebuild() ) )
      build_approximation();
    // update LF model
    if (!sameModelInstance) update_model(lf_model);
    // update lo_fi_set
    lo_fi_set.derivative_vector(set.derivative_vector());
    switch (responseMode) {
    case UNCORRECTED_SURROGATE: case AUTO_CORRECTED_SURROGATE:
    case AGGREGATED_MODELS:
      lo_fi_set.request_vector(lo_fi_asv);             break;
    case MODEL_DISCREPANCY:
      lo_fi_set.request_vector(set.request_vector());  break;
    }
  }

  // HierarchSurrModel's asynchEvalFlag is set if _either_ LF or HF is
  // asynchronous, resulting in use of derived_evaluate_nowait().
  // To manage general case of mixed asynch, launch nonblocking evals first,
  // followed by blocking evals.

  // For notes on repetitive setting of model.solution_level_index(), see
  // derived_evaluate() above.

  // launch nonblocking evals before any blocking ones
  if (hi_fi_eval && asynch_hi_fi) { // HF model may be executed asynchronously
    // don't need to set component parallel mode since only queues the job
    if (sameModelInstance)
      hf_model.solution_level_index(truth_level_index());
    hf_model.evaluate_nowait(hi_fi_set);
    // store map from HF eval id to HierarchSurrModel id
    truthIdMap[hf_model.evaluation_id()] = surrModelEvalCntr;
  }
  if (lo_fi_eval && asynch_lo_fi) { // LF model may be executed asynchronously
    // don't need to set component parallel mode since only queues the job
    if (sameModelInstance)
      lf_model.solution_level_index(surrogate_level_index());
    lf_model.evaluate_nowait(lo_fi_set);
    // store map from LF eval id to HierarchSurrModel id
    surrIdMap[lf_model.evaluation_id()] = surrModelEvalCntr;
    // store variables set needed for correction
    if (responseMode == AUTO_CORRECTED_SURROGATE)
      rawVarsMap[surrModelEvalCntr] = currentVariables.copy();
  }

  // now launch any blocking evals
  if (hi_fi_eval && !asynch_hi_fi) { // execute HF synchronously & cache resp
    component_parallel_mode(TRUTH_MODEL);
    if (sameModelInstance)
      hf_model.solution_level_index(truth_level_index());
    hf_model.evaluate(hi_fi_set);
    // not part of rekey_synch(); can rekey to surrModelEvalCntr immediately
    cachedTruthRespMap[surrModelEvalCntr] = hf_model.current_response().copy();
  }
  if (lo_fi_eval && !asynch_lo_fi) { // execute LF synchronously & cache resp
    component_parallel_mode(SURROGATE_MODEL);
    if (sameModelInstance)
      lf_model.solution_level_index(surrogate_level_index());
    lf_model.evaluate(lo_fi_set);
    Response lo_fi_response(lf_model.current_response().copy());
    // correct LF response prior to caching
    if (responseMode == AUTO_CORRECTED_SURROGATE)
      // correct synch cases now (asynch cases get corrected in
      // derived_synchronize_aggregate*)
      recursive_apply(currentVariables, lo_fi_response);
    // cache corrected LF response for retrieval during synchronization.
    // not part of rekey_synch(); can rekey to surrModelEvalCntr immediately.
    cachedApproxRespMap[surrModelEvalCntr] = lo_fi_response;// deep copied above
  }
}


/** Blocking retrieval of asynchronous evaluations from LF model, HF
    model, or both (mixed case).  For the LF model portion, apply
    correction (if active) to each response in the array.
    derived_synchronize() is designed for the general case where
    derived_evaluate_nowait() may be inconsistent in its use of low
    fidelity evaluations, high fidelity evaluations, or both. */
const IntResponseMap& HierarchSurrModel::derived_synchronize()
{
  surrResponseMap.clear();

  if (sameModelInstance  || sameInterfaceInstance ||
      truthIdMap.empty() || surrIdMap.empty()) { // 1 queue: blocking synch
    IntResponseMap hf_resp_map_rekey, lf_resp_map_rekey;
    derived_synchronize_sequential(hf_resp_map_rekey, lf_resp_map_rekey, true);
    derived_synchronize_combine(hf_resp_map_rekey, lf_resp_map_rekey,
                                surrResponseMap);
  }
  else                               // competing queues: nonblocking synch
    derived_synchronize_competing();

  return surrResponseMap;
}


/** Nonblocking retrieval of asynchronous evaluations from LF model,
    HF model, or both (mixed case).  For the LF model portion, apply
    correction (if active) to each response in the map.
    derived_synchronize_nowait() is designed for the general case
    where derived_evaluate_nowait() may be inconsistent in its use of
    actual evals, approx evals, or both. */
const IntResponseMap& HierarchSurrModel::derived_synchronize_nowait()
{
  surrResponseMap.clear();

  IntResponseMap hf_resp_map_rekey, lf_resp_map_rekey;
  derived_synchronize_sequential(hf_resp_map_rekey, lf_resp_map_rekey, false);
  derived_synchronize_combine_nowait(hf_resp_map_rekey, lf_resp_map_rekey,
                                     surrResponseMap);

  return surrResponseMap;
}


void HierarchSurrModel::
derived_synchronize_sequential(IntResponseMap& hf_resp_map_rekey,
                               IntResponseMap& lf_resp_map_rekey, bool block)
{
  // --------------------------
  // synchronize HF model evals
  // --------------------------
  IntRespMCIter r_cit;
  if (!truthIdMap.empty()) { // synchronize HF evals
    component_parallel_mode(TRUTH_MODEL);
    rekey_synch(truth_model(), block, truthIdMap, hf_resp_map_rekey);
  }
  // add cached truth evals from:
  // (a) recovered HF asynch evals that could not be returned since LF
  //     eval portions were still pending, or
  // (b) synchronous HF evals performed within evaluate_nowait()
  hf_resp_map_rekey.insert(cachedTruthRespMap.begin(),
                           cachedTruthRespMap.end());
  cachedTruthRespMap.clear();

  // --------------------------
  // synchronize LF model evals
  // --------------------------
  if (!surrIdMap.empty()) { // synchronize LF evals
    component_parallel_mode(SURROGATE_MODEL);
    // Interface::rawResponseMap should _not_ be corrected directly since
    // rawResponseMap, beforeSynchCorePRPQueue, and data_pairs all share a
    // responseRep -> modifying rawResponseMap affects data_pairs.
    bool deep_copy = (responseMode == AUTO_CORRECTED_SURROGATE);
    rekey_synch(surrogate_model(), block, surrIdMap, lf_resp_map_rekey,
		deep_copy);
  }
  // add cached approx evals from:
  // (a) recovered LF asynch evals that could not be returned since HF
  //     eval portions were still pending, or
  // (b) synchronous LF evals performed within evaluate_nowait()
  lf_resp_map_rekey.insert(cachedApproxRespMap.begin(),
                           cachedApproxRespMap.end());
  cachedApproxRespMap.clear();
}


void HierarchSurrModel::derived_synchronize_competing()
{
  // in this case, we don't want to starve either LF or HF scheduling by
  // blocking on one or the other --> leverage derived_synchronize_nowait()
  IntResponseMap aggregated_map; // accumulate surrResponseMap returns
  while (!truthIdMap.empty() || !surrIdMap.empty()) {
    // partial_map is a reference to surrResponseMap, returned by _nowait()
    const IntResponseMap& partial_map = derived_synchronize_nowait();
    if (!partial_map.empty())
      aggregated_map.insert(partial_map.begin(), partial_map.end());
  }

  // Note: cached response maps and any LF/HF aggregations are managed
  // within derived_synchronize_nowait()

  std::swap(surrResponseMap, aggregated_map);
}


void HierarchSurrModel::
derived_synchronize_combine(const IntResponseMap& hf_resp_map,
                            IntResponseMap& lf_resp_map,
                            IntResponseMap& combined_resp_map)
{
  // ------------------------------
  // perform any LF/HF aggregations
  // ------------------------------
  // {hf,lf}_resp_map may be partial sets (partial surrogateFnIndices
  // in {UN,AUTO_}CORRECTED_SURROGATE) or full sets (MODEL_DISCREPANCY,
  // AGGREGATED_MODELS).

  IntRespMCIter hf_cit = hf_resp_map.begin(), lf_cit = lf_resp_map.begin();
  bool quiet_flag = (outputLevel < NORMAL_OUTPUT);
  switch (responseMode) {
  case MODEL_DISCREPANCY: {
    DiscrepancyCorrection& delta_corr = deltaCorr[fidelity_keys()];
    for (; hf_cit != hf_resp_map.end() && lf_cit != lf_resp_map.end();
	 ++hf_cit, ++lf_cit) {
      check_key(hf_cit->first, lf_cit->first);
      delta_corr.compute(hf_cit->second, lf_cit->second,
	combined_resp_map[hf_cit->first], quiet_flag);
    }
    break;
  }
  case AGGREGATED_MODELS:
    for (; hf_cit != hf_resp_map.end() && lf_cit != lf_resp_map.end();
	 ++hf_cit, ++lf_cit) {
      check_key(hf_cit->first, lf_cit->first);
      aggregate_response(hf_cit->second, lf_cit->second,
                         combined_resp_map[hf_cit->first]);
    }
    break;
  default: { // {UNCORRECTED,AUTO_CORRECTED,BYPASS}_SURROGATE modes
    if (lf_resp_map.empty()) {
      combined_resp_map = hf_resp_map;  // can't swap w/ const
      return;
    }
    if (responseMode == AUTO_CORRECTED_SURROGATE)
      compute_apply_delta(lf_resp_map);
    if (hf_resp_map.empty()) {
      std::swap(combined_resp_map, lf_resp_map);
      return;
    }
    // process combinations of HF and LF completions
    Response empty_resp;
    while (hf_cit != hf_resp_map.end() || lf_cit != lf_resp_map.end()) {
      // these have been rekeyed already to top-level surrModelEvalCntr:
      int hf_eval_id = (hf_cit == hf_resp_map.end()) ?
                       INT_MAX : hf_cit->first;
      int lf_eval_id = (lf_cit == lf_resp_map.end()) ?
                       INT_MAX : lf_cit->first;

      if (hf_eval_id < lf_eval_id) { // only HF available
        response_combine(hf_cit->second, empty_resp,
                         combined_resp_map[hf_eval_id]);
        ++hf_cit;
      }
      else if (lf_eval_id < hf_eval_id) { // only LF available
        response_combine(empty_resp, lf_cit->second,
                         combined_resp_map[lf_eval_id]);
        ++lf_cit;
      }
      else { // both LF and HF available
        response_combine(hf_cit->second, lf_cit->second,
                         combined_resp_map[hf_eval_id]);
        ++hf_cit;
        ++lf_cit;
      }
    }
    break;
  }
  }
}


void HierarchSurrModel::
derived_synchronize_combine_nowait(const IntResponseMap& hf_resp_map,
                                   IntResponseMap& lf_resp_map,
                                   IntResponseMap& combined_resp_map)
{
  // ------------------------------
  // perform any LF/HF aggregations
  // ------------------------------
  // {hf,lf}_resp_map may be partial sets (partial surrogateFnIndices
  // in {UN,AUTO_}CORRECTED_SURROGATE) or full sets (MODEL_DISCREPANCY).

  // Early return options avoid some overhead:
  if (lf_resp_map.empty() && surrIdMap.empty()) {// none completed, none pending
    combined_resp_map = hf_resp_map;  // can't swap w/ const
    return;
  }
  if (responseMode == AUTO_CORRECTED_SURROGATE)
    compute_apply_delta(lf_resp_map);
  if (hf_resp_map.empty() && truthIdMap.empty()) {//none completed, none pending
    std::swap(combined_resp_map, lf_resp_map);
    return;
  }

  // invert remaining entries (pending jobs) in truthIdMap and surrIdMap
  IntIntMap remain_truth_ids, remain_surr_ids;
  IntIntMCIter id_it;
  for (id_it=truthIdMap.begin(); id_it!=truthIdMap.end(); ++id_it)
    remain_truth_ids[id_it->second] = id_it->first;
  for (id_it=surrIdMap.begin();  id_it!=surrIdMap.end();  ++id_it)
    remain_surr_ids[id_it->second]  = id_it->first;

  // process any combination of HF and LF completions
  IntRespMCIter hf_cit = hf_resp_map.begin();
  IntRespMIter  lf_it  = lf_resp_map.begin();
  Response empty_resp;
  bool quiet_flag = (outputLevel < NORMAL_OUTPUT);
  DiscrepCorrMap::iterator dc_it;
  if (responseMode == MODEL_DISCREPANCY)
    dc_it = deltaCorr.find(fidelity_keys());
  while (hf_cit != hf_resp_map.end() || lf_it != lf_resp_map.end()) {
    // these have been rekeyed already to top-level surrModelEvalCntr:
    int hf_eval_id = (hf_cit == hf_resp_map.end()) ? INT_MAX : hf_cit->first;
    int lf_eval_id = (lf_it  == lf_resp_map.end()) ? INT_MAX : lf_it->first;
    // process LF/HF results or cache them for next pass
    if (hf_eval_id < lf_eval_id) { // only HF available
      switch (responseMode) {
      case MODEL_DISCREPANCY: case AGGREGATED_MODELS:
        // LF contribution is pending -> cache HF response
        cachedTruthRespMap[hf_eval_id] = hf_cit->second;
        break;
      default: // {UNCORRECTED,AUTO_CORRECTED,BYPASS}_SURROGATE modes
        if (remain_surr_ids.find(hf_eval_id) != remain_surr_ids.end())
          // LF contribution is pending -> cache HF response
          cachedTruthRespMap[hf_eval_id] = hf_cit->second;
        else // no LF component is pending -> HF contribution is sufficient
          response_combine(hf_cit->second, empty_resp,
                           surrResponseMap[hf_eval_id]);
        break;
      }
      ++hf_cit;
    }
    else if (lf_eval_id < hf_eval_id) { // only LF available
      switch (responseMode) {
      case MODEL_DISCREPANCY: case AGGREGATED_MODELS:
        // HF contribution is pending -> cache LF response
        cachedApproxRespMap[lf_eval_id] = lf_it->second;
        break;
      default: // {UNCORRECTED,AUTO_CORRECTED,BYPASS}_SURROGATE modes
        if (remain_truth_ids.find(lf_eval_id) != remain_truth_ids.end())
          // HF contribution is pending -> cache LF response
          cachedApproxRespMap[lf_eval_id] = lf_it->second;
        else // no HF component is pending -> LF contribution is sufficient
          response_combine(empty_resp, lf_it->second,
                           surrResponseMap[lf_eval_id]);
        break;
      }
      ++lf_it;
    }
    else { // both LF and HF available
      bool cache_for_pending_corr = false;
      switch (responseMode) {
      case MODEL_DISCREPANCY: {
        dc_it->second.compute(hf_cit->second, lf_it->second,
			      surrResponseMap[hf_eval_id], quiet_flag);
        break;
      }
      case AGGREGATED_MODELS:
        aggregate_response(hf_cit->second, lf_it->second,
                           surrResponseMap[hf_eval_id]);
        break;
      default: // {UNCORRECTED,AUTO_CORRECTED,BYPASS}_SURROGATE modes
        response_combine(hf_cit->second, lf_it->second,
                         surrResponseMap[hf_eval_id]);
        break;
      }
      ++hf_cit;
      ++lf_it;
    }
  }
}


void HierarchSurrModel::compute_apply_delta(IntResponseMap& lf_resp_map)
{
  // Incoming we have a completed LF evaluation that may be used to compute a
  // correction and may be the target of application of a correction.

  // First, test if a correction is previously available or can now be computed
  DiscrepancyCorrection& delta_corr = deltaCorr[fidelity_keys()];
  bool corr_comp = delta_corr.computed(), cache_for_pending_corr = false,
      quiet_flag = (outputLevel < NORMAL_OUTPUT);
  if (!corr_comp) {
    // compute a correction corresponding to the first entry in rawVarsMap
    IntVarsMCIter v_corr_cit = rawVarsMap.begin();
    if (v_corr_cit != rawVarsMap.end()) {
      // if corresponding LF response is complete, compute the delta
      IntRespMCIter lf_corr_cit = lf_resp_map.find(v_corr_cit->first);
      if (lf_corr_cit != lf_resp_map.end()) {
        delta_corr.compute(v_corr_cit->second,truthResponseRef[truthModelKey],
			   lf_corr_cit->second, quiet_flag);
        corr_comp = true;
      }
    }
  }

  // Next, apply the correction.  We cache an uncorrected eval when the
  // components necessary for correction are still pending (returning
  // corrected evals with the first available LF response would lead to
  // nondeterministic results).
  IntVarsMIter v_it; IntRespMIter lf_it; int lf_eval_id;
  for (lf_it=lf_resp_map.begin(); lf_it!=lf_resp_map.end(); ++lf_it) {
    lf_eval_id = lf_it->first;
    v_it = rawVarsMap.find(lf_eval_id);
    if (v_it != rawVarsMap.end()) {
      if (corr_comp) { // apply the correction to the LF response
	recursive_apply(v_it->second, lf_it->second);
        rawVarsMap.erase(v_it);
      }
      else // no new corrections can be applied -> cache uncorrected
        cachedApproxRespMap.insert(*lf_it);
    }
    // else correction already applied
  }
  // remove cached responses from lf_resp_map
  if (!corr_comp)
    for (lf_it =cachedApproxRespMap.begin();
         lf_it!=cachedApproxRespMap.end(); ++lf_it)
      lf_resp_map.erase(lf_it->first);
}


void HierarchSurrModel::
single_apply(const Variables& vars, Response& resp, const UShortArrayPair& keys)
{
  bool quiet_flag = (outputLevel < NORMAL_OUTPUT);
  bool apply_corr = true;
  DiscrepancyCorrection& delta_corr = deltaCorr[keys];
  if (!delta_corr.computed()) {
    std::map<UShortArray, Response>::iterator it
      = truthResponseRef.find(keys.second);
    if (it == truthResponseRef.end()) apply_corr = false; // not found
    else
      delta_corr.compute(vars, truthResponseRef[keys.second], resp, quiet_flag);
  }
  if (apply_corr)
    delta_corr.apply(vars, resp, quiet_flag);
}


void HierarchSurrModel::recursive_apply(const Variables& vars, Response& resp)
{
  switch (correctionMode) {
  case SINGLE_CORRECTION: case DEFAULT_CORRECTION: {
    UShortArrayPair corr_index(surrModelKey, truthModelKey);
    single_apply(vars, resp, corr_index);
    break;
  }
  case FULL_MODEL_FORM_CORRECTION: {
    // assume a consistent level index from surrModelKey
    size_t i, num_models = orderedModels.size();
    UShortArrayPair corr_index(surrModelKey, surrModelKey);
    for (i = surrModelKey.front(); i < num_models - 1; ++i) {
      corr_index.first[0] = i; corr_index.second[0] = i+1;
      single_apply(vars, resp, corr_index);
    }
    break;
  }
  case FULL_SOLUTION_LEVEL_CORRECTION: {
    // assume a consistent model index from surrModelKey[0]
    if (surrModelKey.size() < 2) {
      Cerr << "Error: FULL_SOLUTION_LEVEL_CORRECTION requires solution level "
	   << "within model key." << std::endl;
      abort_handler(MODEL_ERROR);
    }
    size_t i, num_levels = surrogate_model().solution_levels();
    UShortArrayPair corr_index(surrModelKey, surrModelKey);
    for (i = surrModelKey[1]; i < num_levels - 1; ++i) {
      corr_index.first[1] = i; corr_index.second[1] = i+1;
      single_apply(vars, resp, corr_index);
    }
    break;
  }
  case SEQUENCE_CORRECTION: // Apply sequence of discrepancy corrections
    for (size_t i = 0; i < corrSequence.size(); ++i)
      single_apply(vars, resp, corrSequence[i]);
    break;
  }
}


void HierarchSurrModel::resize_response(bool use_virtual_counts)
{
  size_t num_surr, num_truth;
  if (use_virtual_counts) { // allow models to consume lower-level aggregations
    num_surr  = surrogate_model().qoi();
    num_truth =     truth_model().qoi();
  }
  else { // raw counts align with currentResponse raw count
    num_surr  = surrogate_model().response_size();
    num_truth =     truth_model().response_size();
  }

  switch (responseMode) {
  case AGGREGATED_MODELS:
    numFns = num_surr + num_truth;  break;
  case MODEL_DISCREPANCY:
    if (num_surr != num_truth) {
      Cerr << "Error: mismatch in response sizes for MODEL_DISCREPANCY mode "
	   << "in HierarchSurrModel::resize_response()." << std::endl;
      abort_handler(MODEL_ERROR);
    }
    numFns = num_truth;  break;
  case BYPASS_SURROGATE:       case NO_SURROGATE:
    numFns = num_truth;  break;
  case UNCORRECTED_SURROGATE:  case AUTO_CORRECTED_SURROGATE:  default:
    numFns = num_surr;   break;
  }

  // gradient and Hessian settings are based on independent spec (not LF, HF)
  // --> preserve previous settings
  if (currentResponse.num_functions() != numFns) {
    currentResponse.reshape(numFns, currentVariables.cv(),
                            !currentResponse.function_gradients().empty(),
                            !currentResponse.function_hessians().empty());

    // update message lengths for send/receive of parallel jobs (normally
    // performed once in Model::init_communicators() just after construct time)
    //estimate_message_lengths();
    //
    // NOT NECESSARY: Model::synchronize() and Model::serve_run() delegate to
    // HierarchSurrModel::{derived_synchronize,serve_run}() which delegate to
    // synchronize() and serve_run() by the LF or HF model.
    // --> Jobs are never returned using messages containing the expanded
    //     Response object.  Expansion by combination only happens on
    //     iteratorCommRank 0 within derived_synchronize_combine{,_nowait}().
  }
}


void HierarchSurrModel::component_parallel_mode(short mode)
{
  // mode may be correct, but can't guarantee active parallel config is in sync
  //if (componentParallelMode == mode)
  //  return; // already in correct parallel mode

  // terminate previous serve mode (if active)
  UShortArray new_key;
  switch (mode) {
  case SURROGATE_MODEL: new_key =  surrModelKey; break;
  case TRUTH_MODEL:     new_key = truthModelKey; break;
  //default: // leave new_key empty
  }
  // TO DO: restarting servers for a change in soln control index w/o change
  // in model may be overkill (send of state vars in vars buffer sufficient?)
  bool restart = false;
  if (componentParallelMode != mode || componentParallelKey != new_key) {
    if (componentParallelMode) stop_model(componentParallelKey.front());
    restart = true;
  }

  // set ParallelConfiguration for new mode and retrieve new data
  if (mode == TRUTH_MODEL) { // new mode
    // activation delegated to HF model
  }
  else if (mode == SURROGATE_MODEL) { // new mode
    // activation delegated to LF model
  }

  // activate new serve mode (matches HierarchSurrModel::serve_run(pl_iter)).
  // These bcasts match the outer parallel context (pl_iter).
  if (restart && modelPCIter->mi_parallel_level_defined(miPLIndex)) {
    const ParallelLevel& mi_pl = modelPCIter->mi_parallel_level(miPLIndex);
    if (mi_pl.server_communicator_size() > 1) {
      parallelLib.bcast(mode, mi_pl);
      if (mode) { // send model index state corresponding to active mode
	MPIPackBuffer send_buff;
	if      (mode == SURROGATE_MODEL) send_buff << surrModelKey;
	else if (mode == TRUTH_MODEL)     send_buff << truthModelKey;
	send_buff << responseMode;
 	parallelLib.bcast(send_buff, mi_pl);
      }
    }
  }

  componentParallelMode = mode;  componentParallelKey = new_key;
}


void HierarchSurrModel::serve_run(ParLevLIter pl_iter, int max_eval_concurrency)
{
  set_communicators(pl_iter, max_eval_concurrency, false); // don't recurse

  // manage LF model and HF model servers, matching communication from
  // HierarchSurrModel::component_parallel_mode()
  // Note: could consolidate logic by bcasting componentParallelKey,
  //       except for special handling of responseMode for TRUTH_MODEL.
  componentParallelMode = 1;
  while (componentParallelMode) {
    parallelLib.bcast(componentParallelMode, *pl_iter); // outer context
    if (componentParallelMode) {
      // use a quick size estimation for recv buffer i/o size bcast
      MPIPackBuffer send_buff;
       // Note: LF approxDataKeys get appended but not model keys
      UShortArray model_key(2,0);
      send_buff << model_key << responseMode;
      int buffer_len = send_buff.size();
      // receive model state from HierarchSurrModel::component_parallel_mode()
      MPIUnpackBuffer recv_buffer(buffer_len);
      parallelLib.bcast(recv_buffer, *pl_iter);
      recv_buffer >> model_key >> responseMode;

      if (componentParallelMode == SURROGATE_MODEL) {
	// update model key
	surrogate_model_key(model_key); // set LF model + soln index
	// serve active LF model
	surrogate_model().serve_run(pl_iter, max_eval_concurrency);
	// Note: ignores erroneous BYPASS_SURROGATE
      }
      else if (componentParallelMode == TRUTH_MODEL) {
	// update model key
	truth_model_key(model_key); // set HF model + soln index
	// serve active HF model, employing correct iterator concurrency
	Model& hf_model = truth_model();
	switch (responseMode) {
	case UNCORRECTED_SURROGATE:
	  Cerr << "Error: setting parallel mode to TRUTH_MODEL is erroneous "
	       << "for a response mode of UNCORRECTED_SURROGATE." << std::endl;
	  abort_handler(-1);                                              break;
	case AUTO_CORRECTED_SURROGATE:
	  hf_model.serve_run(pl_iter, hf_model.derivative_concurrency()); break;
	case BYPASS_SURROGATE: case MODEL_DISCREPANCY: case AGGREGATED_MODELS:
	  hf_model.serve_run(pl_iter, max_eval_concurrency);              break;
	}
      }
    }
  }
}


void HierarchSurrModel::init_model(Model& model)
{
  // Set the low/high fidelity model variable descriptors with the variable
  // descriptors from currentVariables (eliminates the need to replicate
  // variable descriptors in the input file).  This only needs to be performed
  // once (as opposed to the other updates above).  However, performing this set
  // in the constructor does not propagate properly for multiple surrogates/
  // nestings since the sub-model construction (and therefore any sub-sub-model
  // constructions) must finish before calling any set functions on it.  That
  // is, after-the-fact updating in constructors only propagates one level,
  // whereas before-the-fact updating in compute/build functions propagates
  // multiple levels.
  if (!approxBuilds) {
    size_t num_cv  = currentVariables.cv(),  num_div = currentVariables.div(),
           num_drv = currentVariables.drv(), num_dsv = currentVariables.dsv();
    if (num_cv && num_cv == model.cv())
      model.continuous_variable_labels(
	currentVariables.continuous_variable_labels());
    if (num_div && num_div == model.div())
      model.discrete_int_variable_labels(
        currentVariables.discrete_int_variable_labels());
    if (num_drv && num_drv == model.drv())
      model.discrete_real_variable_labels(
        currentVariables.discrete_real_variable_labels());
    if (num_dsv && num_dsv == model.dsv())
      model.discrete_string_variable_labels(
        currentVariables.discrete_string_variable_labels());
  }

  // linear constraints
  if ( ( userDefinedConstraints.num_linear_ineq_constraints() || 
	 userDefinedConstraints.num_linear_eq_constraints() ) &&
       currentVariables.cv()  == model.cv()  &&
       currentVariables.div() == model.div() &&
       currentVariables.drv() == model.drv() ) {
    model.linear_ineq_constraint_coeffs(
      userDefinedConstraints.linear_ineq_constraint_coeffs());
    model.linear_ineq_constraint_lower_bounds(
      userDefinedConstraints.linear_ineq_constraint_lower_bounds());
    model.linear_ineq_constraint_upper_bounds(
      userDefinedConstraints.linear_ineq_constraint_upper_bounds());

    model.linear_eq_constraint_coeffs(
      userDefinedConstraints.linear_eq_constraint_coeffs());
    model.linear_eq_constraint_targets(
      userDefinedConstraints.linear_eq_constraint_targets());
  }

  // nonlinear constraints
  if (userDefinedConstraints.num_nonlinear_ineq_constraints()) {
    model.nonlinear_ineq_constraint_lower_bounds(
      userDefinedConstraints.nonlinear_ineq_constraint_lower_bounds());
    model.nonlinear_ineq_constraint_upper_bounds(
      userDefinedConstraints.nonlinear_ineq_constraint_upper_bounds());
  }
  if (userDefinedConstraints.num_nonlinear_eq_constraints())
    model.nonlinear_eq_constraint_targets(
      userDefinedConstraints.nonlinear_eq_constraint_targets());

  short active_view = currentVariables.view().first;
  if (active_view == RELAXED_ALL || active_view == MIXED_ALL)
    return;

  // update model with inactive currentVariables/userDefinedConstraints data.
  // For efficiency, we avoid doing this on every evaluation, instead calling
  // it from a pre-execution initialization context (initialize_mapping()).
  size_t num_icv  = currentVariables.icv(),  num_idiv = currentVariables.idiv(),
         num_idrv = currentVariables.idrv(), num_idsv = currentVariables.idsv();
  if (num_icv && num_icv == model.icv()) {
    model.inactive_continuous_variables(
      currentVariables.inactive_continuous_variables());
    model.inactive_continuous_lower_bounds(
      userDefinedConstraints.inactive_continuous_lower_bounds());
    model.inactive_continuous_upper_bounds(
      userDefinedConstraints.inactive_continuous_upper_bounds());
    if (!approxBuilds)
      model.inactive_continuous_variable_labels(
        currentVariables.inactive_continuous_variable_labels());
  }
  if (num_idiv && num_idiv == model.idiv()) {
    model.inactive_discrete_int_variables(
      currentVariables.inactive_discrete_int_variables());
    model.inactive_discrete_int_lower_bounds(
      userDefinedConstraints.inactive_discrete_int_lower_bounds());
    model.inactive_discrete_int_upper_bounds(
      userDefinedConstraints.inactive_discrete_int_upper_bounds());
    if (!approxBuilds)
      model.inactive_discrete_int_variable_labels(
        currentVariables.inactive_discrete_int_variable_labels());
  }
  if (num_idrv && num_idrv == model.idrv()) {
    model.inactive_discrete_real_variables(
      currentVariables.inactive_discrete_real_variables());
    model.inactive_discrete_real_lower_bounds(
      userDefinedConstraints.inactive_discrete_real_lower_bounds());
    model.inactive_discrete_real_upper_bounds(
      userDefinedConstraints.inactive_discrete_real_upper_bounds());
    if (!approxBuilds)
      model.inactive_discrete_real_variable_labels(
        currentVariables.inactive_discrete_real_variable_labels());
  }
  if (num_idsv && num_idsv == model.idsv()) {
    model.inactive_discrete_string_variables(
      currentVariables.inactive_discrete_string_variables());
    if (!approxBuilds)
      model.inactive_discrete_string_variable_labels(
        currentVariables.inactive_discrete_string_variable_labels());
  }
}


void HierarchSurrModel::update_model(Model& model)
{
  // update model with currentVariables/userDefinedConstraints data.  In the
  // hierarchical case, the variables view in LF/HF models correspond to the
  // currentVariables view.  Note: updating the bounds is not strictly necessary
  // in common usage for the HF model (when a single model evaluated only at the
  // TR center), but is needed for the LF model and could be relevant in cases
  // where the HF model involves additional surrogates/nestings.

  // active variable vals/bnds (active labels, inactive vals/bnds/labels, and
  // linear/nonlinear constraint coeffs/bnds updated in init_model())
  if (currentVariables.cv()) {
    model.continuous_variables(currentVariables.continuous_variables());
    model.continuous_lower_bounds(
      userDefinedConstraints.continuous_lower_bounds());
    model.continuous_upper_bounds(
      userDefinedConstraints.continuous_upper_bounds());
  }
  if (currentVariables.div()) {
    model.discrete_int_variables(currentVariables.discrete_int_variables());
    model.discrete_int_lower_bounds(
      userDefinedConstraints.discrete_int_lower_bounds());
    model.discrete_int_upper_bounds(
      userDefinedConstraints.discrete_int_upper_bounds());
  }
  if (currentVariables.drv()) {
    model.discrete_real_variables(currentVariables.discrete_real_variables());
    model.discrete_real_lower_bounds(
      userDefinedConstraints.discrete_real_lower_bounds());
    model.discrete_real_upper_bounds(
      userDefinedConstraints.discrete_real_upper_bounds());
  }
  if (currentVariables.dsv())
    model.discrete_string_variables(
      currentVariables.discrete_string_variables());
}


void HierarchSurrModel::update_from_model(Model& model)
{
  // update complement of active currentVariables using model data.

  // Note: this approach makes a strong assumption about non-active variable
  // consistency, which is limiting.  Better to perform an individual variable
  // mapping (e.g., solution control) when needed and allow for a different
  // ADV position.

  // active variable vals/bnds (active labels, inactive vals/bnds/labels, and
  // linear/nonlinear constraint coeffs/bnds updated in init_model())

  // *** TO DO: make this robust to differing inactive parameterizations using 
  // tag lookups.  Omit mappings for failed lookups.

  const Variables&   vars = model.current_variables();
  const Constraints& cons = model.user_defined_constraints();

  const RealVector& acv = vars.all_continuous_variables();
  StringMultiArrayConstView acv_labels = vars.all_continuous_variable_labels();
  const RealVector& acv_l_bnds = cons.all_continuous_lower_bounds();
  const RealVector& acv_u_bnds = cons.all_continuous_upper_bounds();
  StringMultiArrayConstView cv_acv_labels
    = currentVariables.all_continuous_variable_labels();
  size_t i, index, cv_begin = vars.cv_start(), num_cv = vars.cv(),
    cv_end = cv_begin + num_cv, num_acv = vars.acv();
  for (i=0; i<cv_begin; ++i) {
    index = find_index(cv_acv_labels, acv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_continuous_variable(acv[i], index);
      userDefinedConstraints.all_continuous_lower_bound(acv_l_bnds[i], index);
      userDefinedConstraints.all_continuous_upper_bound(acv_u_bnds[i], index);
    }
  }
  for (i=cv_end; i<num_acv; ++i) {
    index = find_index(cv_acv_labels, acv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_continuous_variable(acv[i], index);
      userDefinedConstraints.all_continuous_lower_bound(acv_l_bnds[i], index);
      userDefinedConstraints.all_continuous_upper_bound(acv_u_bnds[i], index);
    }
  }

  const IntVector& adiv = vars.all_discrete_int_variables();
  StringMultiArrayConstView adiv_labels
    = vars.all_discrete_int_variable_labels();
  const IntVector& adiv_l_bnds = cons.all_discrete_int_lower_bounds();
  const IntVector& adiv_u_bnds = cons.all_discrete_int_upper_bounds();
  StringMultiArrayConstView cv_adiv_labels
    = currentVariables.all_discrete_int_variable_labels();
  size_t div_begin = vars.div_start(), num_div = vars.div(),
    div_end = div_begin + num_div, num_adiv = vars.adiv();
  for (i=0; i<div_begin; ++i) {
    index = find_index(cv_adiv_labels, adiv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_discrete_int_variable(adiv[i], index);
      userDefinedConstraints.all_discrete_int_lower_bound(adiv_l_bnds[i],index);
      userDefinedConstraints.all_discrete_int_upper_bound(adiv_u_bnds[i],index);
    }
  }
  for (i=div_end; i<num_adiv; ++i) {
    index = find_index(cv_adiv_labels, adiv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_discrete_int_variable(adiv[i], index);
      userDefinedConstraints.all_discrete_int_lower_bound(adiv_l_bnds[i],index);
      userDefinedConstraints.all_discrete_int_upper_bound(adiv_u_bnds[i],index);
    }
  }

  size_t dsv_begin = vars.dsv_start(), num_dsv = vars.dsv(),
    dsv_end = dsv_begin + num_dsv, num_adsv = vars.adsv();
  StringMultiArrayConstView adsv = vars.all_discrete_string_variables();
  StringMultiArrayConstView adsv_labels
    = vars.all_discrete_string_variable_labels();
  StringMultiArrayConstView cv_adsv_labels
    = currentVariables.all_discrete_string_variable_labels();
  for (i=0; i<dsv_begin; ++i) {
    index = find_index(cv_adsv_labels, adsv_labels[i]);
    if (index != _NPOS)
      currentVariables.all_discrete_string_variable(adsv[i], index);
  }
  for (i=dsv_end; i<num_adsv; ++i) {
    index = find_index(cv_adsv_labels, adsv_labels[i]);
    if (index != _NPOS)
      currentVariables.all_discrete_string_variable(adsv[i], index);
  }

  const RealVector& adrv = vars.all_discrete_real_variables();
  StringMultiArrayConstView adrv_labels
    = vars.all_discrete_real_variable_labels();
  const RealVector& adrv_l_bnds = cons.all_discrete_real_lower_bounds();
  const RealVector& adrv_u_bnds = cons.all_discrete_real_upper_bounds();
  StringMultiArrayConstView cv_adrv_labels
    = currentVariables.all_discrete_real_variable_labels();
  size_t drv_begin = vars.drv_start(), num_drv = vars.drv(),
    drv_end = drv_begin + num_drv, num_adrv = vars.adrv();
  for (i=0; i<drv_begin; ++i) {
    index = find_index(cv_adrv_labels, adrv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_discrete_real_variable(adrv[i], index);
      userDefinedConstraints.all_discrete_real_lower_bound(adrv_l_bnds[i],
							   index);
      userDefinedConstraints.all_discrete_real_upper_bound(adrv_u_bnds[i],
							   index);
    }
  }
  for (i=drv_end; i<num_adrv; ++i) {
    index = find_index(cv_adrv_labels, adrv_labels[i]);
    if (index != _NPOS) {
      currentVariables.all_discrete_real_variable(adrv[i], index);
      userDefinedConstraints.all_discrete_real_lower_bound(adrv_l_bnds[i],
							   index);
      userDefinedConstraints.all_discrete_real_upper_bound(adrv_u_bnds[i],
							   index);
    }
  }
}

} // namespace Dakota
