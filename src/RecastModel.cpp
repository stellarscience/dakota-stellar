/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2020
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       RecastModel
//- Description: Implementation code for the RecastModel class
//- Owner:       Mike Eldred
//- Checked by:

#include "dakota_system_defs.hpp"
#include "RecastModel.hpp"
#include "EvaluationStore.hpp"

static const char rcsId[]="@(#) $Id: RecastModel.cpp 7029 2010-10-22 00:17:02Z mseldre $";


using namespace std;

namespace Dakota {

//#define DEBUG

// init static var
StringStringPairIntMap RecastModel::recastModelIdCounters;


/** Default recast model constructor.  Requires full definition of the
    transformation; if any mappings are NULL, they are assumed to
    remain so in later initialization or updates.  Parameter
    vars_comps_totals indicates the number of each type of variable {4
    types} x {3 domains} in the recast variable space.  Note:
    recast_secondary_offset is the start index for equality
    constraints, typically num nonlinear ineq constraints. */
RecastModel::
RecastModel(const Model& sub_model, const Sizet2DArray& vars_map_indices,
	    const SizetArray& vars_comps_totals, const BitArray& all_relax_di,
	    const BitArray& all_relax_dr, bool nonlinear_vars_mapping,
	    void (*variables_map)      (const Variables& recast_vars,
					Variables& sub_model_vars),
	    void (*set_map)            (const Variables& recast_vars,
					const ActiveSet& recast_set,
					ActiveSet& sub_model_set),
	    const Sizet2DArray& primary_resp_map_indices,
	    const Sizet2DArray& secondary_resp_map_indices,
	    size_t recast_secondary_offset, short recast_resp_order,
	    const BoolDequeArray& nonlinear_resp_mapping,
	    void (*primary_resp_map)   (const Variables& sub_model_vars,
					const Variables& recast_vars,
					const Response& sub_model_response,
					Response& recast_response),
	    void (*secondary_resp_map) (const Variables& sub_model_vars,
					const Variables& recast_vars,
					const Response& sub_model_response,
					Response& recast_response)):
  Model(LightWtBaseConstructor(), sub_model.problem_description_db(),
	sub_model.parallel_library()),
  subModel(sub_model), varsMapIndices(vars_map_indices),
  nonlinearVarsMapping(nonlinear_vars_mapping),
  primaryRespMapIndices(primary_resp_map_indices),
  secondaryRespMapIndices(secondary_resp_map_indices),
  nonlinearRespMapping(nonlinear_resp_mapping), recastModelEvalCntr(0),
  variablesMapping(variables_map), setMapping(set_map),
  primaryRespMapping(primary_resp_map),
  secondaryRespMapping(secondary_resp_map), invVarsMapping(NULL),
  invSetMapping(NULL), invPriRespMapping(NULL), invSecRespMapping(NULL)
{
  modelType = "recast"; 
  supportsEstimDerivs = false; // subModel estimates derivatives by default

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();

  // recasting of variables; only reshape if change in variable type counts
  bool reshape_vars = false;
  if (variablesMapping) {
    // reshape as dictated by variable type changes
    reshape_vars = 
      init_variables(vars_comps_totals, all_relax_di, all_relax_dr); 
  }
  else {
    // variables are not mapped: deep copy of vars to allow independence, but 
    // shallow copy of svd since types/labels/ids can be kept consistent
    currentVariables = subModel.current_variables().copy(); // shared svd
    numDerivVars = currentVariables.cv();
  }

  if (nonlinearRespMapping.size() != 
      primaryRespMapIndices.size() + secondaryRespMapIndices.size()) {
    Cerr << "Error: size mismatch in response mapping configuration." << endl;
    abort_handler(-1);
  }

  if (primaryRespMapping || secondaryRespMapping) {
    init_response(primaryRespMapIndices.size(), secondaryRespMapIndices.size(), 
		  recast_resp_order, reshape_vars);
  }
  else {
    currentResponse = subModel.current_response().copy();
    numFns = currentResponse.num_functions();
  }

  init_constraints(secondaryRespMapIndices.size(), 
		   recast_secondary_offset, reshape_vars);

  modelId = RecastModel::recast_model_id(root_model_id(), "RECAST");
}


/** This alternate constructor defers initialization of the function
    pointers until a separate call to initialize(), and accepts the
    minimum information needed to construct currentVariables,
    currentResponse, and userDefinedConstraints.  The resulting model
    is sufficiently complete for passing to an Iterator.  Parameter
    vars_comps_totals indicates the number of each type of variable {4
    types} x {3 domains} in the recast variable space. Note:
    recast_secondary_offset is the start index for equality
    constraints, typically num nonlinear ineq constraints. */
RecastModel::
RecastModel(const Model& sub_model, //size_t num_deriv_vars,
	    const SizetArray& vars_comps_totals, const BitArray& all_relax_di,
	    const BitArray& all_relax_dr,    size_t num_recast_primary_fns,
	    size_t num_recast_secondary_fns, size_t recast_secondary_offset,
	    short recast_resp_order):
  Model(LightWtBaseConstructor(), sub_model.problem_description_db(),
	sub_model.parallel_library()),
  subModel(sub_model), nonlinearVarsMapping(false), recastModelEvalCntr(0),
  variablesMapping(NULL), setMapping(NULL), primaryRespMapping(NULL),
  secondaryRespMapping(NULL), invVarsMapping(NULL), invSetMapping(NULL),
  invPriRespMapping(NULL), invSecRespMapping(NULL)
{
  modelType = "recast";
  supportsEstimDerivs = false; // subModel estimates derivatives by default

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();

  // initialize Variables, Response, and Constraints based on sizes
  init_sizes(vars_comps_totals, all_relax_di, all_relax_dr,
	     num_recast_primary_fns, num_recast_secondary_fns,
	     recast_secondary_offset, recast_resp_order);
  modelId = RecastModel::recast_model_id(root_model_id(), "RECAST");
}


RecastModel::RecastModel(ProblemDescDB& problem_db, const Model& sub_model):
  Model(BaseConstructor(), problem_db), subModel(sub_model),
  recastModelEvalCntr(0), variablesMapping(NULL), setMapping(NULL),
  primaryRespMapping(NULL), secondaryRespMapping(NULL), invVarsMapping(NULL),
  invSetMapping(NULL), invPriRespMapping(NULL), invSecRespMapping(NULL)
{
  modelType = "recast";
  supportsEstimDerivs = false; // subModel estimates derivatives by default

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();
  modelId = RecastModel::recast_model_id(root_model_id(), "RECAST");
}


RecastModel::RecastModel(const Model& sub_model):
  Model(LightWtBaseConstructor(), sub_model.problem_description_db(),
   	sub_model.parallel_library()),
  subModel(sub_model), recastModelEvalCntr(0), variablesMapping(NULL),
  setMapping(NULL), primaryRespMapping(NULL), secondaryRespMapping(NULL),
  invVarsMapping(NULL), invSetMapping(NULL), invPriRespMapping(NULL),
  invSecRespMapping(NULL)
{ 
  modelType = "recast";
  supportsEstimDerivs = false; // subModel estimates derivatives by default

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();
  numFns = sub_model.response_size();
  modelId = RecastModel::recast_model_id(root_model_id(), "RECAST");
}


void RecastModel::
init_sizes(const SizetArray& vars_comps_totals, const BitArray& all_relax_di,
	   const BitArray& all_relax_dr,    size_t num_recast_primary_fns,
	   size_t num_recast_secondary_fns, size_t recast_secondary_offset,
	   short recast_resp_order)
{
  // recasting of variables; only reshape if change in variable type counts
  bool reshape_vars =
    init_variables(vars_comps_totals, all_relax_di, all_relax_dr); 

  // Currently this reshape is subspace-specific and handled in derived classes
  //if (reshape_vars) // else full-space initialization is sufficient
  //  init_distribution(rv_types, vars_comps_totals); // reshape mvDist

  // recasting of response and constraints
  init_response(num_recast_primary_fns, num_recast_secondary_fns, 
		recast_resp_order, reshape_vars);

  init_constraints(num_recast_secondary_fns,
		   recast_secondary_offset, reshape_vars);
}


/** This function is used for late initialization of the recasting
    functions.  It is used in concert with the alternate constructor. */
void RecastModel::
init_maps(const Sizet2DArray& vars_map_indices,
	  bool nonlinear_vars_mapping,
	  void (*variables_map)      (const Variables& recast_vars,
				      Variables& sub_model_vars),
	  void (*set_map)            (const Variables& recast_vars,
				      const ActiveSet& recast_set,
				      ActiveSet& sub_model_set),
	  const Sizet2DArray& primary_resp_map_indices,
	  const Sizet2DArray& secondary_resp_map_indices,
	  const BoolDequeArray& nonlinear_resp_mapping,
	  void (*primary_resp_map)   (const Variables& sub_model_vars,
				      const Variables& recast_vars,
				      const Response& sub_model_response,
				      Response& recast_response),
	  void (*secondary_resp_map) (const Variables& sub_model_vars,
				      const Variables& recast_vars,
				      const Response& sub_model_response,
				      Response& recast_response))
{
  varsMapIndices          = vars_map_indices;
  nonlinearVarsMapping    = nonlinear_vars_mapping;
  variablesMapping        = variables_map;
  setMapping              = set_map;
  primaryRespMapIndices   = primary_resp_map_indices;
  secondaryRespMapIndices = secondary_resp_map_indices;
  nonlinearRespMapping    = nonlinear_resp_mapping;
  primaryRespMapping      = primary_resp_map;
  secondaryRespMapping    = secondary_resp_map;

  if (nonlinearRespMapping.size() != primaryRespMapIndices.size() +
      secondaryRespMapIndices.size()) {
    Cerr << "Error: size mismatch in response mapping configuration." << endl;
    abort_handler(-1);
  }
}


short RecastModel::response_order(const Model& sub_model)
{
  const Response& curr_resp = sub_model.current_response();

  short recast_resp_order = 1; // recast resp order to be same as original resp
  if (!curr_resp.function_gradients().empty()) recast_resp_order |= 2;
  if (!curr_resp.function_hessians().empty())  recast_resp_order |= 4;

  return recast_resp_order;
}


String RecastModel::recast_model_id(const String &root_id, const String &type) {
  auto key = std::make_pair(root_id, type);
  int id;
  if(recastModelIdCounters.find(key) == recastModelIdCounters.end())
    recastModelIdCounters[key] = id = 1;
  else
    id = ++recastModelIdCounters[key];
  return String("RECAST_") + root_id + "_" + type + "_" + std::to_string(id);
}


bool RecastModel::
init_variables(const SizetArray& vars_comps_totals,
	       const BitArray& all_relax_di, const BitArray& all_relax_dr)
{
  const Variables& sub_model_vars = subModel.current_variables();
  const SharedVariablesData&  svd = sub_model_vars.shared_data();

  // BMA: We actually don't allow the case of a change in
  // vars_comp_totals, but no mapping, but have to allow it here in
  // case mapping not yet provided.

  // if any change in variable types, will need a new SharedVariablesData
  bool vars_char_same = 
    ( vars_comps_totals.empty() || 
      svd.components_totals()         == vars_comps_totals ) &&
    ( all_relax_di.empty() ||
      svd.all_relaxed_discrete_int()  == all_relax_di )      &&
    ( all_relax_dr.empty() || 
      svd.all_relaxed_discrete_real() == all_relax_dr );

  // check change in character first as mapping may not yet be present...
  if (vars_char_same) {
    // variables are mapped but not resized: deep copy of vars and
    // same svd, since types may change in transformed space
    currentVariables = sub_model_vars.copy(true); // independent svd
  }
  else {
    // variables are resized; need new SVD regardless
    SharedVariablesData recast_svd(sub_model_vars.view(), vars_comps_totals,
				   all_relax_di, all_relax_dr);
    currentVariables = Variables(recast_svd);
  }

  // propagate number of active continuous vars to derivative vars
  numDerivVars = currentVariables.cv();

  return !vars_char_same; // return reshape_vars
}


void RecastModel::
init_response(size_t num_recast_primary_fns, size_t num_recast_secondary_fns,
	      short recast_resp_order, bool reshape_vars)
{
  numFns = num_recast_primary_fns + num_recast_secondary_fns;

  // recasting of response
  const Response& sub_model_resp = subModel.current_response();
  currentResponse = sub_model_resp.copy();

  bool grad_flag = (recast_resp_order & 2),
    hess_flag = (recast_resp_order & 4),
    sm_grad_flag = !sub_model_resp.function_gradients().empty(),
    sm_hess_flag = !sub_model_resp.function_hessians().empty();
  const Variables& sub_model_vars = subModel.current_variables();
  if ( sub_model_vars.cv()            != numDerivVars ||
       sub_model_resp.num_functions() != numFns       ||
       grad_flag != sm_grad_flag || hess_flag != sm_hess_flag )
    currentResponse.reshape(numFns, numDerivVars, grad_flag, hess_flag);
}


void RecastModel::
reshape_response(size_t num_recast_primary_fns, size_t num_recast_secondary_fns)
{
  numFns = num_recast_primary_fns + num_recast_secondary_fns;
  bool grad_flag = !currentResponse.function_gradients().empty();
  bool hess_flag = !currentResponse.function_hessians().empty();

  currentResponse.reshape(numFns, numDerivVars, grad_flag, hess_flag);
}


void RecastModel::
init_constraints(size_t num_recast_secondary_fns,
		 size_t recast_secondary_offset, bool reshape_vars)
{
  // recasting of constraints
  SharedVariablesData recast_svd = currentVariables.shared_data();
  const Constraints& sub_model_cons = subModel.user_defined_constraints();
  userDefinedConstraints = (reshape_vars) ?
    Constraints(recast_svd) : sub_model_cons.copy();

  // the recast_secondary_offset cannot in general be inferred from the
  // contributing fns in secondaryRespMapIndices (recast constraints may be
  // defined, e.g., with no contributing fns), and must therefore be passed.
  size_t num_recast_nln_ineq = recast_secondary_offset,
    num_recast_nln_eq = num_recast_secondary_fns - num_recast_nln_ineq;
  if ( num_recast_nln_ineq != sub_model_cons.num_nonlinear_ineq_constraints()
       || num_recast_nln_eq   != sub_model_cons.num_nonlinear_eq_constraints() )
    userDefinedConstraints.reshape(num_recast_nln_ineq, num_recast_nln_eq,
      sub_model_cons.num_linear_ineq_constraints(),
      sub_model_cons.num_linear_eq_constraints());
}


void RecastModel::inverse_mappings(
    void (*inv_vars_map)     (const Variables& recast_vars,
			      Variables& sub_model_vars),
    void (*inv_set_map)      (const Variables& recast_vars,
			      const ActiveSet& recast_set,
			      ActiveSet& sub_model_set),
    void (*inv_pri_resp_map) (const Variables& sub_model_vars,
			      const Variables& recast_vars,
			      const Response& sub_model_resp,
			      Response& recast_resp),
    void (*inv_sec_resp_map) (const Variables& sub_model_vars,
			      const Variables& recast_vars,
			      const Response& sub_model_resp,
			      Response& recast_resp))
{
  invVarsMapping    = inv_vars_map;     invSetMapping     = inv_set_map;
  invPriRespMapping = inv_pri_resp_map; invSecRespMapping = inv_sec_resp_map;
}


/** The RecastModel is evaluated by an Iterator for a recast problem
    formulation.  Therefore, the currentVariables, incoming active set,
    and output currentResponse all correspond to the recast inputs/outputs. */
void RecastModel::derived_evaluate(const ActiveSet& set)
{
  ++recastModelEvalCntr;

  // transform from recast (Iterator) to sub-model (user) variables
  transform_variables(currentVariables, subModel.current_variables());

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(currentVariables, set, sub_model_set);

  // evaluate the subModel in the original fn set definition.  Doing this here 
  // eliminates the need for eval tracking logic within the separate eval fns.
  subModel.evaluate(sub_model_set);

  // recast the subModel response ("user space") into the currentResponse
  // ("iterator space")
  currentResponse.active_set(set);
  if (primaryRespMapping || secondaryRespMapping)
    transform_response(currentVariables, subModel.current_variables(),
		       subModel.current_response(), currentResponse);
  else
    currentResponse.update(subModel.current_response());

#ifdef DEBUG
  Cout << "Recast variables:\n"   << currentVariables
       << "subModel variables:\n" << subModel.current_variables()
       << "subModel response:\n"  << subModel.current_response()
       << "Recast response:\n"    << currentResponse;
#endif
}


void RecastModel::derived_evaluate_nowait(const ActiveSet& set)
{
  ++recastModelEvalCntr;

  // transform from recast (Iterator) to sub-model (user) variables
  transform_variables(currentVariables, subModel.current_variables());

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(currentVariables, set, sub_model_set);

  // evaluate the subModel in the original fn set definition.  Doing this here 
  // eliminates the need for eval tracking logic within the separate eval fns.
  subModel.evaluate_nowait(sub_model_set);
  // in almost all cases, use of the subModel eval ids is sufficient, but
  // protect against the rare case where not all subModel evaluations being
  // scheduled were spawned from the RecastModel (e.g., a HierarchicalModel
  // that uses an ActiveSubspaceModel as LF and the original model as HF).
  recastIdMap[subModel.evaluation_id()] = recastModelEvalCntr;

  // bookkeep variables for use in primaryRespMapping/secondaryRespMapping
  if (primaryRespMapping || secondaryRespMapping) {
    recastSetMap[recastModelEvalCntr]  = set;
    recastVarsMap[recastModelEvalCntr] = currentVariables.copy();
    if (variablesMapping)
      subModelVarsMap[recastModelEvalCntr]
	= subModel.current_variables().copy();
  }
}


const IntResponseMap& RecastModel::derived_synchronize()
{
  recastResponseMap.clear();

  if (primaryRespMapping || secondaryRespMapping) {
    IntResponseMap resp_map_rekey;
    rekey_synch(subModel, true, recastIdMap, resp_map_rekey);
    transform_response_map(resp_map_rekey, recastResponseMap);
  }
  else
    rekey_synch(subModel, true, recastIdMap, recastResponseMap);

  return recastResponseMap;
}


const IntResponseMap& RecastModel::derived_synchronize_nowait()
{
  recastResponseMap.clear();

  if (primaryRespMapping || secondaryRespMapping) {
    IntResponseMap resp_map_rekey;
    rekey_synch(subModel, false, recastIdMap, resp_map_rekey);
    transform_response_map(resp_map_rekey, recastResponseMap);
  }
  else
    rekey_synch(subModel, false, recastIdMap, recastResponseMap);

  return recastResponseMap;
}


void RecastModel::
transform_response_map(const IntResponseMap& old_resp_map,
		       IntResponseMap& new_resp_map)
{
  IntRespMCIter r_cit; IntASMIter s_it; IntVarsMIter v_it, sm_v_it;
  for (r_cit=old_resp_map.begin(); r_cit!=old_resp_map.end(); ++r_cit) {
    int native_id = r_cit->first;
    s_it =  recastSetMap.find(native_id);
    v_it = recastVarsMap.find(native_id);
    if (variablesMapping) sm_v_it = subModelVarsMap.find(native_id);
    else                  sm_v_it = v_it;

    Response new_resp(currentResponse.copy()); // correct size, labels, etc.
    new_resp.active_set(s_it->second);
    transform_response(v_it->second, sm_v_it->second, r_cit->second, new_resp);
    new_resp_map[native_id] = new_resp;

    // cleanup
    recastSetMap.erase(s_it);  recastVarsMap.erase(v_it);
    if (variablesMapping) subModelVarsMap.erase(sm_v_it);
  }
}


void RecastModel::
transform_variables(const Variables& recast_vars, Variables& sub_model_vars)
{
  // typical flow: mapping from recast variables ("iterator space")
  // into the sub-model variables ("user space")
  if (variablesMapping) {
    assign_instance();
    variablesMapping(recast_vars, sub_model_vars);
  }
  else
    sub_model_vars.active_variables(recast_vars);
}


void RecastModel::
inverse_transform_variables(const Variables& sub_model_vars,
			    Variables& recast_vars)
{
  // atypical flow: mapping from sub-model variables ("user space")
  // into the recast variables ("iterator space")
  if (invVarsMapping) {
    assign_instance();
    invVarsMapping(sub_model_vars, recast_vars);
  }
  else
    recast_vars.active_variables(sub_model_vars);
}


void RecastModel::
transform_set(const Variables& recast_vars, const ActiveSet& recast_set,
	      ActiveSet& sub_model_set)
{
  // typical flow: mapping from recast set ("iterator space") into the
  // sub-model set ("user space")

  size_t i, j, num_recast_primary_fns = primaryRespMapIndices.size(),
    num_recast_secondary_fns = secondaryRespMapIndices.size(),
    num_recast_fns = num_recast_primary_fns + num_recast_secondary_fns;
  const ShortArray& recast_asv = recast_set.request_vector();
  if (recast_asv.size() != num_recast_fns) {
    Cerr << "Error: inconsistent asv sizing in RecastModel::transform_set().\n"
	 << "       recast asv size = " << recast_asv.size() << '\n'
	 << "       recast functions = " << num_recast_fns << endl;
    abort_handler(-1);
  }

  // Define default request vector and derivative vector mappings:
  // For the ASV, project each recast_asv request onto the contributing
  // set of functions within the sub_model_asv.  In the case of nonlinear
  // input/output mappings, the recast_asv request is augmented with
  // additional data requirements derived from chain rule differentiation.
  // The default sub-model DVV is just a copy of the recast DVV.
  ShortArray sub_model_asv(subModel.response_size(), 0);
  for (i=0; i<num_recast_fns; i++) {
    short asv_val = recast_asv[i];
    // For nonlinear variable mappings, gradient required to transform Hessian.
    // A single nonlinear variable mapping affects all function derivatives.
    if (nonlinearVarsMapping && (asv_val & 4))
      asv_val |= 2;
    // assign the asv_val to each contributing sub-model function
    const SizetArray& recast_fn_contributors = (i<num_recast_primary_fns) ?
      primaryRespMapIndices[i] :
      secondaryRespMapIndices[i-num_recast_primary_fns];
    size_t num_contributors = recast_fn_contributors.size();
    for (j=0; j<num_contributors; j++) {
      short sub_model_asv_val = asv_val;
      // Bit deletions: for NLS recasting for full Newton without LeastSq term
      // Hessians, could remove 4 bit based on {gradient,hessian}Type, but this
      // is better accomplished from an Iterator's configuration using the
      // setMapping plug-in below (e.g., see Optimizer::gnewton_set_recast()). 

      // Bit additions: for nonlinear resp mappings, derivatives require all
      // lower order data. The nonlinearity of each fn contribution is employed.
      if (nonlinearRespMapping[i][j]) {
	if (asv_val & 4)
	  sub_model_asv_val |= 3;
	else if (asv_val & 2)
	  sub_model_asv_val |= 1;
      }
      sub_model_asv[recast_fn_contributors[j]] |= sub_model_asv_val;
    }
  }
  sub_model_set.request_vector(sub_model_asv);
  sub_model_set.derivative_vector(recast_set.derivative_vector()); // copy

  // a setMapping (provided in the RecastModel ctor or initialize()) augments
  // the standard mappings.  Current examples include NonD::set_u_to_x_mapping,
  // NonDReliability::PMA2_set_mapping, and Optimizer::gauss_newton_set_recast.
  // This follows the standard mappings so that provided mappings don't get
  // overwritten by the standard logic.  However, this means that any provided
  // additions will not be automatically augmented by nonlinear mapping logic
  // above.  This should not be a significant problem, since the provided
  // additions have case-specific context whereas the logic above is generic.
  // It would be preferable if provided mappings focused on updating the
  // sub_model_set rather than generating it from recast_set.
  if (setMapping) {
    assign_instance();
    setMapping(recast_vars, recast_set, sub_model_set);
  }
}


void RecastModel::
inverse_transform_set(const Variables& sub_model_vars,
		      const ActiveSet& sub_model_set, ActiveSet& recast_set)
{
  // atypical flow: mapping from sub-model set ("user space") into the
  // recast set ("iterator space")

  /* TO DO: modify mapping below from forward to inverse

  size_t i, j, num_recast_primary_fns = primaryRespMapIndices.size(),
    num_recast_secondary_fns = secondaryRespMapIndices.size(),
    num_recast_fns = num_recast_primary_fns + num_recast_secondary_fns;
  const ShortArray& recast_asv = recast_set.request_vector();
  if (recast_asv.size() != num_recast_fns) {
    Cerr << "Error: inconsistent asv sizing in RecastModel::"
         << "inverse_transform_set()." << std::endl;
    abort_handler(-1);
  }

  // Define default request vector and derivative vector mappings:
  // For the ASV, project each recast_asv request onto the contributing
  // set of functions within the sub_model_asv.  In the case of nonlinear
  // input/output mappings, the recast_asv request is augmented with
  // additional data requirements derived from chain rule differentiation.
  // The default sub-model DVV is just a copy of the recast DVV.
  ShortArray sub_model_asv(subModel.response_size(), 0);
  for (i=0; i<num_recast_fns; i++) {
    short asv_val = recast_asv[i];
    // For nonlinear variable mappings, gradient required to transform Hessian.
    // A single nonlinear variable mapping affects all function derivatives.
    if (nonlinearVarsMapping && (asv_val & 4))
      asv_val |= 2;
    // assign the asv_val to each contributing sub-model function
    const SizetArray& recast_fn_contributors = (i<num_recast_primary_fns) ?
      primaryRespMapIndices[i] :
      secondaryRespMapIndices[i-num_recast_primary_fns];
    size_t num_contributors = recast_fn_contributors.size();
    for (j=0; j<num_contributors; j++) {
      short sub_model_asv_val = asv_val;
      // Bit deletions: for NLS recasting for full Newton without LeastSq term
      // Hessians, could remove 4 bit based on {gradient,hessian}Type, but this
      // is better accomplished from an Iterator's configuration using the
      // setMapping plug-in below (e.g., see Optimizer::gnewton_set_recast()). 

      // Bit additions: for nonlinear resp mappings, derivatives require all
      // lower order data. The nonlinearity of each fn contribution is employed.
      if (nonlinearRespMapping[i][j]) {
	if (asv_val & 4)
	  sub_model_asv_val |= 3;
	else if (asv_val & 2)
	  sub_model_asv_val |= 1;
      }
      sub_model_asv[recast_fn_contributors[j]] |= sub_model_asv_val;
    }
  }
  sub_model_set.request_vector(sub_model_asv);
  sub_model_set.derivative_vector(recast_set.derivative_vector()); // copy
  */

  // an invSetMapping (provided in inverse_mappings()) augments the standard
  // mappings above, such that the provided mappings don't get overwritten by
  // the standard logic.
  if (invSetMapping) {
    assign_instance();
    invSetMapping(sub_model_vars, sub_model_set, recast_set);
  }
}


void RecastModel::
transform_response(const Variables& recast_vars,
		   const Variables& sub_model_vars,
		   const Response& sub_model_resp, Response& recast_resp)
{
  // typical flow: mapping from sub-model response ("user space") into
  // the recast response ("iterator space")

  size_t num_recast_1_fns = primaryRespMapIndices.size();

  if (primaryRespMapping || secondaryRespMapping)
    assign_instance();

  if (primaryRespMapping)
    primaryRespMapping(sub_model_vars, recast_vars,
		       sub_model_resp, recast_resp);
  else // number of recast primary = number of sub-model primary
    recast_resp.update_partial(0, num_recast_1_fns, sub_model_resp, 0);

  if (secondaryRespMapping)
    secondaryRespMapping(sub_model_vars, recast_vars,
			 sub_model_resp, recast_resp);
  else {
    // number of recast secondary = number of sub-model secondary,
    // but primary offsets may differ
    size_t num_recast_2_fns = secondaryRespMapIndices.size(),
           num_sm_1_fns     = sub_model_resp.num_functions() - num_recast_2_fns;
    recast_resp.update_partial(num_recast_1_fns, num_recast_2_fns,
			       sub_model_resp, num_sm_1_fns);
  }
}


void RecastModel::
inverse_transform_response(const Variables& sub_model_vars,
			   const Variables& recast_vars,
			   const Response& recast_resp,
			   Response& sub_model_resp)
{
  // atypical flow: mapping from the recast response ("iterator space")
  // into the sub-model response ("user space")

  size_t num_recast_1_fns = primaryRespMapIndices.size();

  if (invPriRespMapping || invSecRespMapping)
    assign_instance();

  if (invPriRespMapping)
    invPriRespMapping(recast_vars, sub_model_vars, recast_resp, sub_model_resp);
  else // number of recast primary = number of sub-model primary
    sub_model_resp.update_partial(0, num_recast_1_fns, recast_resp, 0);

  if (invSecRespMapping)
    invSecRespMapping(recast_vars, sub_model_vars, recast_resp, sub_model_resp);
  else {
    // number of recast secondary = number of sub-model secondary,
    // but primary offsets may differ
    size_t num_recast_2_fns = secondaryRespMapIndices.size(),
           num_sm_1_fns     = sub_model_resp.num_functions() - num_recast_2_fns;
    sub_model_resp.update_partial(num_sm_1_fns, num_recast_2_fns,
				  recast_resp, num_recast_1_fns);
  }
}


void RecastModel::initialize_data_from_submodel()
{
  componentParallelMode = SUB_MODEL_MODE;
  outputLevel           = subModel.output_level();

  gradientType          = subModel.gradient_type();
  methodSource          = subModel.method_source();
  ignoreBounds          = subModel.ignore_bounds();
  centralHess	        = subModel.central_hess();
  intervalType          = subModel.interval_type();
  fdGradStepSize        = subModel.fd_gradient_step_size();
  fdGradStepType        = subModel.fd_gradient_step_type();
  gradIdAnalytic        = subModel.gradient_id_analytic();
  gradIdNumerical       = subModel.gradient_id_numerical();

  hessianType           = subModel.hessian_type();
  quasiHessType         = subModel.quasi_hessian_type();
  fdHessByFnStepSize    = subModel.fd_hessian_by_fn_step_size();
  fdHessByGradStepSize  = subModel.fd_hessian_by_grad_step_size();
  fdHessStepType        = subModel.fd_hessian_step_type();
  hessIdAnalytic        = subModel.hessian_id_analytic();
  hessIdNumerical       = subModel.hessian_id_numerical();
  hessIdQuasi           = subModel.hessian_id_quasi();

  scalingOpts           = subModel.scaling_options();
}


/** Update inactive values and labels in currentVariables and inactive
    bound constraints in userDefinedConstraints from variables and
    constraints data within subModel. */
void RecastModel::update_from_model(Model& model)
{
  // break up into pieces so that derived Recasts can override with subsets

  bool update_active_complement = update_variables_from_model(model);

  if (update_active_complement)
    update_variables_active_complement_from_model(model);

  update_response_from_model(model);
}


bool RecastModel::update_variables_from_model(Model& model)
{
  bool update_active_complement = true;
  if (invVarsMapping) { // inv mapping provided: sub-model -> recast
    assign_instance();

    // generally restricted to active variables
    invVarsMapping(model.current_variables(), currentVariables);

    // BMA TODO: there may be cases where we also want to update the
    // constraints and values, but there's currently no mechanism to
    // do so.  The client of a RecastModel must manage this.
  }
  else if (variablesMapping) { // only fwd mapping for recast -> sub-model
    // no reasonable default for active vars

    // can't just apply variables mapping to values/bounds, since need inverse
    // of variablesMapping to go from model vars to currentVariables

    // any label, uncertain variable distributions, and linear
    // constraint mappings must be performed explicitly

    // for partial mapping of variables that are unmodified by a variable
    // transformation, see NonDExpansion::initialize_expansion()

    // for mapping of distribution parameters that are part of a variable
    // transformation, see ProbabilityTransformModel::update_transformation()
  }
  else {
    update_active_complement = false; // can use all view updates below

    update_variable_values(model);
    update_variable_bounds(model);
    update_variable_labels(model);

    // uncertain variable distribution data
    // > deep copies were used previously for Pecos::DistributionParams
    //mvDist.update(model.multivariate_distribution());
    // Current approach: rep is shared
    // > tramples an mvDist construction from Model(BaseConstructor) ....
    // > populates mvDist for Model(LightWtBaseConstructor)
    // > reassignments protected by smart ptr management
    // Note: becomes less important w/ broader use of ProbabilityTransformModel
    mvDist = subModel.multivariate_distribution(); // shared rep

    update_linear_constraints(model);
  }

  return update_active_complement;
}


void RecastModel::update_variable_values(const Model& model)
{
  currentVariables.all_continuous_variables
    (model.all_continuous_variables());
  update_discrete_variable_values(model);
}


void RecastModel::update_discrete_variable_values(const Model& model)
{
  currentVariables.all_discrete_int_variables
    (model.all_discrete_int_variables());
  currentVariables.all_discrete_string_variables
    (model.all_discrete_string_variables());
  currentVariables.all_discrete_real_variables
    (model.all_discrete_real_variables());
}


void RecastModel::update_variable_bounds(const Model& model)
{
  userDefinedConstraints.all_continuous_lower_bounds
    (model.all_continuous_lower_bounds());
  userDefinedConstraints.all_continuous_upper_bounds
    (model.all_continuous_upper_bounds());
  update_discrete_variable_bounds(model);
}


void RecastModel::update_discrete_variable_bounds(const Model& model)
{
  userDefinedConstraints.all_discrete_int_lower_bounds
    (model.all_discrete_int_lower_bounds());
  userDefinedConstraints.all_discrete_int_upper_bounds
    (model.all_discrete_int_upper_bounds());
  userDefinedConstraints.all_discrete_real_lower_bounds
    (model.all_discrete_real_lower_bounds());
  userDefinedConstraints.all_discrete_real_upper_bounds
    (model.all_discrete_real_upper_bounds());
}


void RecastModel::update_variable_labels(const Model& model)
{
  currentVariables.all_continuous_variable_labels
    ( model.all_continuous_variable_labels());
  update_discrete_variable_labels(model);
}


void RecastModel::update_discrete_variable_labels(const Model& model)
{
  currentVariables.all_discrete_int_variable_labels
    (model.all_discrete_int_variable_labels());
  currentVariables.all_discrete_string_variable_labels
    (model.all_discrete_string_variable_labels());
  currentVariables.all_discrete_real_variable_labels
    (model.all_discrete_real_variable_labels());
}


void RecastModel::update_linear_constraints(const Model& model)
{
  if (model.num_linear_ineq_constraints()) {
    userDefinedConstraints.linear_ineq_constraint_coeffs
      (model.linear_ineq_constraint_coeffs());
    userDefinedConstraints.linear_ineq_constraint_lower_bounds
      (model.linear_ineq_constraint_lower_bounds());
    userDefinedConstraints.linear_ineq_constraint_upper_bounds
      (model.linear_ineq_constraint_upper_bounds());
  }
  if (model.num_linear_eq_constraints()) {
    userDefinedConstraints.linear_eq_constraint_coeffs
      (model.linear_eq_constraint_coeffs());
    userDefinedConstraints.linear_eq_constraint_targets
      (model.linear_eq_constraint_targets());
  }
}


void RecastModel::update_variables_active_complement_from_model(Model& model)
{
  size_t i, cv_begin = currentVariables.cv_start(),
    num_cv  = currentVariables.cv(), cv_end = cv_begin + num_cv,
    num_acv = currentVariables.acv();
  const RealVector& acv = model.all_continuous_variables();
  const RealVector& acv_l_bnds = model.all_continuous_lower_bounds();
  const RealVector& acv_u_bnds = model.all_continuous_upper_bounds();
  StringMultiArrayConstView acv_labels
    = model.all_continuous_variable_labels();
  for (i=0; i<cv_begin; ++i) {
    currentVariables.all_continuous_variable(acv[i], i);
    userDefinedConstraints.all_continuous_lower_bound(acv_l_bnds[i], i);
    userDefinedConstraints.all_continuous_upper_bound(acv_u_bnds[i], i);
    currentVariables.all_continuous_variable_label(acv_labels[i], i);
  }
  for (i=cv_end; i<num_acv; ++i) {
    currentVariables.all_continuous_variable(acv[i], i);
    userDefinedConstraints.all_continuous_lower_bound(acv_l_bnds[i], i);
    userDefinedConstraints.all_continuous_upper_bound(acv_u_bnds[i], i);
    currentVariables.all_continuous_variable_label(acv_labels[i], i);
  }

  size_t div_begin = currentVariables.div_start(),
    num_div  = currentVariables.div(), div_end = div_begin + num_div,
    num_adiv = currentVariables.adiv();
  const IntVector& adiv = model.all_discrete_int_variables();
  const IntVector& adiv_l_bnds = model.all_discrete_int_lower_bounds();
  const IntVector& adiv_u_bnds = model.all_discrete_int_upper_bounds();
  StringMultiArrayConstView adiv_labels
    = model.all_discrete_int_variable_labels();
  for (i=0; i<div_begin; ++i) {
    currentVariables.all_discrete_int_variable(adiv[i], i);
    userDefinedConstraints.all_discrete_int_lower_bound(adiv_l_bnds[i], i);
    userDefinedConstraints.all_discrete_int_upper_bound(adiv_u_bnds[i], i);
    currentVariables.all_discrete_int_variable_label(adiv_labels[i], i);
  }
  for (i=div_end; i<num_adiv; ++i) {
    currentVariables.all_discrete_int_variable(adiv[i], i);
    userDefinedConstraints.all_discrete_int_lower_bound(adiv_l_bnds[i], i);
    userDefinedConstraints.all_discrete_int_upper_bound(adiv_u_bnds[i], i);
    currentVariables.all_discrete_int_variable_label(adiv_labels[i], i);
  }

  size_t dsv_begin = currentVariables.dsv_start(),
    num_dsv  = currentVariables.dsv(), dsv_end = dsv_begin + num_dsv,
    num_adsv = currentVariables.adsv();
  StringMultiArrayConstView adsv = model.all_discrete_string_variables();
  StringMultiArrayConstView adsv_labels
    = model.all_discrete_string_variable_labels();
  for (i=0; i<dsv_begin; ++i) {
    currentVariables.all_discrete_string_variable(adsv[i], i);
    currentVariables.all_discrete_string_variable_label(adsv_labels[i], i);
  }
  for (i=dsv_end; i<num_adsv; ++i) {
    currentVariables.all_discrete_string_variable(adsv[i], i);
    currentVariables.all_discrete_string_variable_label(adsv_labels[i], i);
  }

  size_t drv_begin = currentVariables.drv_start(),
    num_drv  = currentVariables.drv(), drv_end = drv_begin + num_drv,
    num_adrv = currentVariables.adrv();
  const RealVector& adrv = model.all_discrete_real_variables();
  const RealVector& adrv_l_bnds = model.all_discrete_real_lower_bounds();
  const RealVector& adrv_u_bnds = model.all_discrete_real_upper_bounds();
  StringMultiArrayConstView adrv_labels
    = model.all_discrete_real_variable_labels();
  for (i=0; i<drv_begin; ++i) {
    currentVariables.all_discrete_real_variable(adrv[i], i);
    userDefinedConstraints.all_discrete_real_lower_bound(adrv_l_bnds[i], i);
    userDefinedConstraints.all_discrete_real_upper_bound(adrv_u_bnds[i], i);
    currentVariables.all_discrete_real_variable_label(adrv_labels[i], i);
  }
  for (i=drv_end; i<num_adrv; ++i) {
    currentVariables.all_discrete_real_variable(adrv[i], i);
    userDefinedConstraints.all_discrete_real_lower_bound(adrv_l_bnds[i], i);
    userDefinedConstraints.all_discrete_real_upper_bound(adrv_u_bnds[i], i);
    currentVariables.all_discrete_real_variable_label(adrv_labels[i], i);
  }
}


void RecastModel::update_response_from_model(Model& model)
{
  if (primaryRespMapping) {
    // response mappings are in opposite direction from variables
    // mappings, so primaryRespMapping could potentially be used to
    // update currentResponse from model primary fns
  }
  else {
    // primary response function weights
    primaryRespFnWts = model.primary_response_fn_weights();
    // primary response function sense (min or max)
    primaryRespFnSense = model.primary_response_fn_sense();

    // primary response function labels
    const StringArray& sm_resp_labels = model.response_labels();
    size_t i, num_primary = numFns 
      - userDefinedConstraints.num_nonlinear_eq_constraints()
      - userDefinedConstraints.num_nonlinear_ineq_constraints();
    for (i=0; i<num_primary; i++)
      currentResponse.shared_data().function_label(sm_resp_labels[i], i);
  }

  if (secondaryRespMapping) {
    // response mappings are in opposite direction from variables
    // mappings, so secondaryRespMapping could potentially be used to
    // update currentResponse from model secondary fns
  }
  else {
    update_secondary_response(model);
  }
}


void RecastModel::update_secondary_response(const Model& model)
{
  // secondary response function labels
  const StringArray& sm_resp_labels = model.response_labels();
  size_t i,
    num_nln_con = userDefinedConstraints.num_nonlinear_eq_constraints() +
    userDefinedConstraints.num_nonlinear_ineq_constraints(),
    num_primary    = numFns - num_nln_con,
    num_sm_primary = model.response_size() - num_nln_con;
  for (i=0; i<num_nln_con; i++)
    currentResponse.shared_data().function_label
      (sm_resp_labels[num_sm_primary+i], num_primary+i);

  // nonlinear constraint bounds/targets
  if (model.num_nonlinear_ineq_constraints()) {
    userDefinedConstraints.nonlinear_ineq_constraint_lower_bounds
      (model.nonlinear_ineq_constraint_lower_bounds());
    userDefinedConstraints.nonlinear_ineq_constraint_upper_bounds
      (model.nonlinear_ineq_constraint_upper_bounds());
  }
  if (model.num_nonlinear_eq_constraints())
    userDefinedConstraints.nonlinear_eq_constraint_targets
      (model.nonlinear_eq_constraint_targets());
}


const RealVector& RecastModel::error_estimates()
{
  // linear mappings are fine (see NestedModel), but general nonlinear mappings
  // are problematic.
  if (primaryRespMapping || secondaryRespMapping) {

    // preclude nonlinear mappings and multi-component mappings for now.
    // Note: a linear multi-component mapping can be supported by NestedModel::
    //       iterator_error_estimation(), because the mapping coeffs are known.

    size_t i, num_recast_fns = nonlinearRespMapping.size();
    for (i=0; i<num_recast_fns; ++i) {
      const BoolDeque& nln_resp_map_i = nonlinearRespMapping[i];
      if (nln_resp_map_i.size() > 1 ||
	  std::find(nln_resp_map_i.begin(), nln_resp_map_i.end(), true) !=
	            nln_resp_map_i.end()) {
	Cerr << "Error: error estimation not currently supported for Recast"
	     << "Model with nonlinear or multi-component response mapping."
	     << std::endl;
	abort_handler(MODEL_ERROR);
      }
    }

    // push errors through linear single-component mapping (individual scaling)

    // make dummy responses for use with transform_response()
    const Response& sm_resp = subModel.current_response();
    ActiveSet sm_set = sm_resp.active_set(),
          recast_set = currentResponse.active_set();
    sm_set.request_values(1); recast_set.request_values(1);
    Response sm_error_est(sm_resp.shared_data(), sm_set),
         recast_error_est(currentResponse.shared_data(), recast_set);
    // transform the error estimates as Response::functionValues
    sm_error_est.function_values(subModel.error_estimates());
    if (outputLevel >= DEBUG_OUTPUT) // distinguish ScalingModel debug blocks
      Cout << "Transforming Error Estimates:\n";
    transform_response(currentVariables, subModel.current_variables(),
		       sm_error_est, recast_error_est);
    mappedErrorEstimates = recast_error_est.function_values();
    return mappedErrorEstimates;
  }
  else
    return subModel.error_estimates();
}


void RecastModel::resize_response_mapping()
{
  // reshaping primaryRespMapIndices, secondaryRespMapIndices
  size_t num_curr_fns       = response_size(),
         num_curr_secondary = num_secondary_fns(),
         num_curr_primary   = num_curr_fns - num_curr_secondary,
         num_sm_fns         = subModel.response_size(),
         num_sm_secondary   = subModel.num_secondary_fns(),
         num_sm_primary     = num_sm_fns - num_sm_secondary,
         num_replicates, num_ind, offset, i, j, k;

  primaryRespMapIndices.resize(num_sm_primary);
  secondaryRespMapIndices.resize(num_sm_secondary);
  nonlinearRespMapping.resize(num_sm_fns);

  // the number of mappings to replicate corresponds to the number
  // of data sets being aggregated in, e.g., AGGREGATED_MODELS mode

  if (num_sm_primary > num_curr_primary) { // inflate existing mappings
    num_replicates = num_sm_primary / num_curr_primary;
    if (num_sm_primary % num_curr_primary) {
      Cerr << "Error: non-integer multiplier for response aggregation in "
	   << "RecastModel::resize_response_mapping()" << std::endl;
      abort_handler(MODEL_ERROR);
    }
    for (i=0; i<num_curr_primary; ++i) {
      // mirror pattern within new response indices, but offset index mapping
      SizetArray& primary_ind_i = primaryRespMapIndices[i];
      num_ind = primary_ind_i.size();
      for (j=1; j<num_replicates; ++j) {
	offset = j * num_curr_primary;
	SizetArray& primary_ind_j = primaryRespMapIndices[i + offset];
	primary_ind_j.resize(num_ind);
	for (k=0; k<num_ind; ++k)
	  primary_ind_j[k] = primary_ind_i[k] + offset;
      }
    }
  }
  if (num_sm_secondary > num_curr_secondary) { // inflate existing mappings
    num_replicates = num_sm_secondary / num_curr_secondary;
    if (num_sm_secondary % num_curr_secondary) {
      Cerr << "Error: non-integer multiplier for response aggregation in "
	   << "RecastModel::resize_response_mapping()" << std::endl;
      abort_handler(MODEL_ERROR);
    }
    for (i=0; i<num_curr_secondary; ++i) {
      // mirror pattern within new response indices, but offset index mapping
      SizetArray& secondary_ind_i = secondaryRespMapIndices[i];
      num_ind = secondary_ind_i.size();
      for (j=1; j<num_replicates; ++j) {
	offset = j * num_curr_secondary;
	SizetArray& secondary_ind_j = secondaryRespMapIndices[i + offset];
	secondary_ind_j.resize(num_ind);
	for (k=0; k<num_ind; ++k)
	  secondary_ind_j[k]   = secondary_ind_i[k] + offset;
      }
    }
  }
  if (num_sm_fns > num_curr_fns) { // inflate existing mappings
    num_replicates = num_sm_fns / num_curr_fns;
    if (num_sm_fns % num_curr_fns) {
      Cerr << "Error: non-integer multiplier for response aggregation in "
	   << "RecastModel::resize_response_mapping()" << std::endl;
      abort_handler(MODEL_ERROR);
    }
    for (i=0; i<num_curr_fns; ++i) {
      // mirror pattern within new response indices, but offset index mapping
      BoolDeque& nonlin_resp_map_i = nonlinearRespMapping[i];
      for (j=1; j<num_replicates; ++j)
	nonlinearRespMapping[i + j * num_curr_fns] = nonlin_resp_map_i;
    }
  }
}


bool RecastModel::
db_lookup(const Variables& search_vars, const ActiveSet& search_set,
	  Response& found_resp)
{
  // transform from recast (Iterator) to sub-model (user) variables;
  // making copy to avoid modifying submodel state during the lookup
  Variables sub_model_vars(subModel.current_variables().copy());
  transform_variables(search_vars, sub_model_vars);

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(search_vars, search_set, sub_model_set);

  // invoke submodel lookup; making copy to avoid modifying submodel state 
  // during the lookup
  Response sub_model_resp(subModel.current_response().copy());
  // sub_model_resp must have right ASV so lookup's update will pull right data
  sub_model_resp.active_set(sub_model_set);
  bool eval_found
    = subModel.db_lookup(sub_model_vars, sub_model_set, sub_model_resp);
  if (!eval_found)
    return false;

  // recast the subModel response ("user space") into the "iterator space"
  found_resp.active_set(search_set);
  if (primaryRespMapping || secondaryRespMapping)
    transform_response(search_vars, sub_model_vars, sub_model_resp, found_resp);
  else
    found_resp.update(sub_model_resp);

  return eval_found;
}


void RecastModel::assign_instance()
{ } // no static instance pointer to assign at base (default is no-op)


String RecastModel::root_model_id() {
  return subModel.root_model_id();
}


ActiveSet RecastModel::default_active_set() {
  // The "base class" implementation assumes that supportsEstimDerivs is false
  // and that gradients/hessians, if available, are computed by a submodel and
  // hence can be provided by this model.
  ActiveSet set;
  set.derivative_vector(currentVariables.all_continuous_variable_ids());
  bool has_deriv_vars = set.derivative_vector().size() != 0;
  ShortArray asv(numFns, 1);
  if(has_deriv_vars) {
    if(gradientType != "none")// && (gradientType == "analytic" || supportsEstimDerivs))
        for(auto &a : asv)
          a |=  2;

    if(hessianType != "none") // && (hessianType == "analytic" || supportsEstimDerivs))
        for(auto &a : asv)
          a |=  4;
  }
  set.request_vector(asv);
  return set;
}


void RecastModel::declare_sources() {
  evaluationsDB.declare_source(modelId, modelType, subModel.model_id(), subModel.model_type());
}

} // namespace Dakota
