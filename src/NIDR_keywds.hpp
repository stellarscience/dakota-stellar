
namespace Dakota {

static KeyWord
	kw_1[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_stm(augment_utype,postRunInputFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_stm(augment_utype,postRunInputFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_stm(augment_utype,postRunInputFormat_TABULAR_IFACE_ID)}
		},
	kw_2[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_stm(utype,postRunInputFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_1,0.,0.,0,N_stm(utype,postRunInputFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_stm(utype,postRunInputFormat_TABULAR_NONE)}
		},
	kw_3[2] = {
		{"input",11,3,1,0,kw_2,0.,0.,0,N_stm(str,postRunInput)},
		{"output",11,0,2,0,0,0.,0.,0,N_stm(str,postRunOutput)}
		},
	kw_4[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_stm(augment_utype,preRunOutputFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_stm(augment_utype,preRunOutputFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_stm(augment_utype,preRunOutputFormat_TABULAR_IFACE_ID)}
		},
	kw_5[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_stm(utype,preRunOutputFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_4,0.,0.,0,N_stm(utype,preRunOutputFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_stm(utype,preRunOutputFormat_TABULAR_NONE)}
		},
	kw_6[2] = {
		{"input",11,0,1,0,0,0.,0.,0,N_stm(str,preRunInput)},
		{"output",11,3,2,0,kw_5,0.,0.,0,N_stm(str,preRunOutput)}
		},
	kw_7[1] = {
		{"stop_restart",0x29,0,1,0,0,0.,0.,0,N_stm(int,stopRestart)}
		},
	kw_8[3] = {
		{"all",8,0,1,1,0,0.,0.,0,N_stm(utype,interfEvalsSelection_INTERF_EVAL_STORE_ALL)},
		{"none",8,0,1,1,0,0.,0.,0,N_stm(utype,interfEvalsSelection_INTERF_EVAL_STORE_NONE)},
		{"simulation",8,0,1,1,0,0.,0.,0,N_stm(utype,interfEvalsSelection_INTERF_EVAL_STORE_SIMULATION)}
		},
	kw_9[4] = {
		{"all",8,0,1,1,0,0.,0.,0,N_stm(utype,modelEvalsSelection_MODEL_EVAL_STORE_ALL)},
		{"all_methods",8,0,1,1,0,0.,0.,0,N_stm(utype,modelEvalsSelection_MODEL_EVAL_STORE_ALL_METHODS)},
		{"none",8,0,1,1,0,0.,0.,0,N_stm(utype,modelEvalsSelection_MODEL_EVAL_STORE_NONE)},
		{"top_method",8,0,1,1,0,0.,0.,0,N_stm(utype,modelEvalsSelection_MODEL_EVAL_STORE_TOP_METHOD)}
		},
	kw_10[2] = {
		{"interface_selection",8,3,2,0,kw_8},
		{"model_selection",8,4,1,0,kw_9}
		},
	kw_11[3] = {
		{"hdf5",8,2,3,0,kw_10,0.,0.,0,N_stm(augment_utype,resultsOutputFormat_RESULTS_OUTPUT_HDF5)},
		{"results_output_file",11,0,1,0,0,0.,0.,0,N_stm(str,resultsOutputFile)},
		{"text",8,0,2,0,0,0.,0.,0,N_stm(augment_utype,resultsOutputFormat_RESULTS_OUTPUT_TEXT)}
		},
	kw_12[2] = {
		{"input",11,0,1,0,0,0.,0.,0,N_stm(str,runInput)},
		{"output",11,0,2,0,0,0.,0.,0,N_stm(str,runOutput)}
		},
	kw_13[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_stm(augment_utype,tabularFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_stm(augment_utype,tabularFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_stm(augment_utype,tabularFormat_TABULAR_IFACE_ID)}
		},
	kw_14[5] = {
		{"annotated",8,0,2,0,0,0.,0.,0,N_stm(utype,tabularFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,2,0,kw_13,0.,0.,0,N_stm(utype,tabularFormat_TABULAR_NONE)},
		{"freeform",8,0,2,0,0,0.,0.,0,N_stm(utype,tabularFormat_TABULAR_NONE)},
		{"tabular_data_file",11,0,1,0,0,0.,0.,0,N_stm(str,tabularDataFile)},
		{"tabular_graphics_file",3,0,1,0,0,0.,0.,-1,N_stm(str,tabularDataFile)}
		},
	kw_15[15] = {
		{"check",8,0,9,0,0,0.,0.,0,N_stm(true,checkFlag)},
		{"error_file",11,0,3,0,0,0.,0.,0,N_stm(str,errorFile)},
		{"graphics",8,0,8,0,0,0.,0.,0,N_stm(true,graphicsFlag)},
		{"method_pointer",3,0,13,0,0,0.,0.,10,N_stm(str,topMethodPointer)},
		{"output_file",11,0,2,0,0,0.,0.,0,N_stm(str,outputFile)},
		{"output_precision",0x29,0,6,0,0,0.,0.,0,N_stm(int,outputPrecision)},
		{"post_run",8,2,12,0,kw_3,0.,0.,0,N_stm(true,postRunFlag)},
		{"pre_run",8,2,10,0,kw_6,0.,0.,0,N_stm(true,preRunFlag)},
		{"read_restart",11,1,4,0,kw_7,0.,0.,0,N_stm(str,readRestart)},
		{"results_output",8,3,7,0,kw_11,0.,0.,0,N_stm(true,resultsOutputFlag)},
		{"run",8,2,11,0,kw_12,0.,0.,0,N_stm(true,runFlag)},
		{"tabular_data",8,5,1,0,kw_14,0.,0.,0,N_stm(true,tabularDataFlag)},
		{"tabular_graphics_data",0,5,1,0,kw_14,0.,0.,-1,N_stm(true,tabularDataFlag)},
		{"top_method_pointer",11,0,13,0,0,0.,0.,0,N_stm(str,topMethodPointer)},
		{"write_restart",11,0,5,0,0,0.,0.,0,N_stm(str,writeRestart)}
		},
	kw_16[1] = {
		{"processors_per_analysis",0x19,0,1,0,0,0.,0.,0,N_ifm(int,procsPerAnalysis)}
		},
	kw_17[8] = {
		{"copy_files",15,0,5,0,0,0.,0.,0,N_ifm(strL,copyFiles)},
		{"dir_save",0,0,3,0,0,0.,0.,2,N_ifm(true,dirSave)},
		{"dir_tag",0,0,2,0,0,0.,0.,2,N_ifm(true,dirTag)},
		{"directory_save",8,0,3,0,0,0.,0.,0,N_ifm(true,dirSave)},
		{"directory_tag",8,0,2,0,0,0.,0.,0,N_ifm(true,dirTag)},
		{"link_files",15,0,4,0,0,0.,0.,0,N_ifm(strL,linkFiles)},
		{"named",11,0,1,0,0,0.,0.,0,N_ifm(str,workDir)},
		{"replace",8,0,6,0,0,0.,0.,0,N_ifm(true,templateReplace)}
		},
	kw_18[10] = {
		{"allow_existing_results",8,0,8,0,0,0.,0.,0,N_ifm(true,allowExistingResultsFlag)},
		{"aprepro",8,0,6,0,0,0.,0.,0,N_ifm(true,apreproFlag)},
		{"dprepro",0,0,6,0,0,0.,0.,-1,N_ifm(true,apreproFlag)},
		{"file_save",8,0,4,0,0,0.,0.,0,N_ifm(true,fileSaveFlag)},
		{"file_tag",8,0,3,0,0,0.,0.,0,N_ifm(true,fileTagFlag)},
		{"labeled",8,0,5,0,0,0.,0.,0,N_ifm(type,resultsFileFormat_LABELED_RESULTS)},
		{"parameters_file",11,0,1,0,0,0.,0.,0,N_ifm(str,parametersFile)},
		{"results_file",11,0,2,0,0,0.,0.,0,N_ifm(str,resultsFile)},
		{"verbatim",8,0,9,0,0,0.,0.,0,N_ifm(true,verbatimFlag)},
		{"work_directory",8,8,7,0,kw_17,0.,0.,0,N_ifm(true,useWorkdir)}
		},
	kw_19[1] = {
		{"numpy",8,0,1,0,0,0.,0.,0,N_ifm(true,numpyFlag)}
		},
	kw_20[8] = {
		{"copy_files",15,0,5,0,0,0.,0.,0,N_ifm(strL,copyFiles)},
		{"dir_save",0,0,3,0,0,0.,0.,2,N_ifm(true,dirSave)},
		{"dir_tag",0,0,2,0,0,0.,0.,2,N_ifm(true,dirTag)},
		{"directory_save",8,0,3,0,0,0.,0.,0,N_ifm(true,dirSave)},
		{"directory_tag",8,0,2,0,0,0.,0.,0,N_ifm(true,dirTag)},
		{"link_files",15,0,4,0,0,0.,0.,0,N_ifm(strL,linkFiles)},
		{"named",11,0,1,0,0,0.,0.,0,N_ifm(str,workDir)},
		{"replace",8,0,6,0,0,0.,0.,0,N_ifm(true,templateReplace)}
		},
	kw_21[10] = {
		{"allow_existing_results",8,0,8,0,0,0.,0.,0,N_ifm(true,allowExistingResultsFlag)},
		{"aprepro",8,0,6,0,0,0.,0.,0,N_ifm(true,apreproFlag)},
		{"dprepro",0,0,6,0,0,0.,0.,-1,N_ifm(true,apreproFlag)},
		{"file_save",8,0,4,0,0,0.,0.,0,N_ifm(true,fileSaveFlag)},
		{"file_tag",8,0,3,0,0,0.,0.,0,N_ifm(true,fileTagFlag)},
		{"labeled",8,0,5,0,0,0.,0.,0,N_ifm(type,resultsFileFormat_LABELED_RESULTS)},
		{"parameters_file",11,0,1,0,0,0.,0.,0,N_ifm(str,parametersFile)},
		{"results_file",11,0,2,0,0,0.,0.,0,N_ifm(str,resultsFile)},
		{"verbatim",8,0,9,0,0,0.,0.,0,N_ifm(true,verbatimFlag)},
		{"work_directory",8,8,7,0,kw_20,0.,0.,0,N_ifm(true,useWorkdir)}
		},
	kw_22[11] = {
		{"analysis_components",15,0,4,0,0,0.,0.,0,N_ifm(str2D,analysisComponents)},
		{"direct",8,1,3,1,kw_16,0.,0.,0,N_ifm(type,interfaceType_TEST_INTERFACE)},
		{"fork",8,10,3,1,kw_18,0.,0.,0,N_ifm(type,interfaceType_FORK_INTERFACE)},
		{"grid",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_GRID_INTERFACE)},
		{"input_filter",11,0,1,0,0,0.,0.,0,N_ifm(str,inputFilter)},
		{"matlab",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_MATLAB_INTERFACE)},
		{"output_filter",11,0,2,0,0,0.,0.,0,N_ifm(str,outputFilter)},
		{"pybind11",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_PYBIND11_INTERFACE)},
		{"python",8,1,3,1,kw_19,0.,0.,0,N_ifm(type,interfaceType_PYTHON_INTERFACE)},
		{"scilab",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_SCILAB_INTERFACE)},
		{"system",8,10,3,1,kw_21,0.,0.,0,N_ifm(type,interfaceType_SYSTEM_INTERFACE)}
		},
	kw_23[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_ifm(type,analysisScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_ifm(type,analysisScheduling_PEER_SCHEDULING)}
		},
	kw_24[2] = {
		{"dynamic",8,0,1,1,0,0.,0.,0,N_ifm(type,asynchLocalEvalScheduling_DYNAMIC_SCHEDULING)},
		{"static",8,0,1,1,0,0.,0.,0,N_ifm(type,asynchLocalEvalScheduling_STATIC_SCHEDULING)}
		},
	kw_25[3] = {
		{"analysis_concurrency",0x19,0,3,0,0,0.,0.,0,N_ifm(int,asynchLocalAnalysisConcurrency)},
		{"evaluation_concurrency",0x19,0,1,0,0,0.,0.,0,N_ifm(int,asynchLocalEvalConcurrency)},
		{"local_evaluation_scheduling",8,2,2,0,kw_24}
		},
	kw_26[1] = {
		{"size",0x19,0,1,0,0,0.,0.,0,N_ifm(int,asynchLocalEvalConcurrency)}
		},
	kw_27[1] = {
		{"cache_tolerance",10,0,1,0,0,0.,0.,0,N_ifm(Real,nearbyEvalCacheTol)}
		},
	kw_28[4] = {
		{"active_set_vector",8,0,1,0,0,0.,0.,0,N_ifm(false,activeSetVectorFlag)},
		{"evaluation_cache",8,0,2,0,0,0.,0.,0,N_ifm(false,evalCacheFlag)},
		{"restart_file",8,0,4,0,0,0.,0.,0,N_ifm(false,restartFileFlag)},
		{"strict_cache_equality",8,1,3,0,kw_27,0.,0.,0,N_ifm(true,nearbyEvalCacheFlag)}
		},
	kw_29[2] = {
		{"dynamic",8,0,1,1,0,0.,0.,0,N_ifm(type,evalScheduling_PEER_DYNAMIC_SCHEDULING)},
		{"static",8,0,1,1,0,0.,0.,0,N_ifm(type,evalScheduling_PEER_STATIC_SCHEDULING)}
		},
	kw_30[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_ifm(type,evalScheduling_MASTER_SCHEDULING)},
		{"peer",8,2,1,1,kw_29}
		},
	kw_31[4] = {
		{"abort",8,0,1,1,0,0.,0.,0,N_ifm(lit,failAction_abort)},
		{"continuation",8,0,1,1,0,0.,0.,0,N_ifm(lit,failAction_continuation)},
		{"recover",14,0,1,1,0,0.,0.,0,N_ifm(Rlit,TYPE_DATA_failAction_recover)},
		{"retry",9,0,1,1,0,0.,0.,0,N_ifm(ilit,TYPE_DATA_failAction_retry)}
		},
	kw_32[12] = {
		{"algebraic_mappings",11,0,3,0,0,0.,0.,0,N_ifm(str,algebraicMappings)},
		{"analysis_drivers",15,11,2,0,kw_22,0.,0.,0,N_ifm(strL,analysisDrivers)},
		{"analysis_scheduling",8,2,11,0,kw_23},
		{"analysis_servers",0x19,0,10,0,0,0.,0.,0,N_ifm(int,analysisServers)},
		{"asynchronous",8,3,6,0,kw_25,0.,0.,0,N_ifm(true, asynchFlag)},
		{"batch",8,1,6,0,kw_26,0.,0.,0,N_ifm(true,batchEvalFlag)},
		{"deactivate",8,4,5,0,kw_28},
		{"evaluation_scheduling",8,2,8,0,kw_30},
		{"evaluation_servers",0x19,0,7,0,0,0.,0.,0,N_ifm(int,evalServers)},
		{"failure_capture",8,4,4,0,kw_31},
		{"id_interface",11,0,1,0,0,0.,0.,0,N_ifm(str,idInterface)},
		{"processors_per_evaluation",0x19,0,9,0,0,0.,0.,0,N_ifm(int,procsPerEval)}
		},
	kw_33[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_34[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_35[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_34,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_36[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_37[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_38[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_39[3] = {
		{"offline_pilot",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_OFFLINE_PILOT)},
		{"online_pilot",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_ONLINE_PILOT)},
		{"pilot_projection",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_PILOT_PROJECTION)}
		},
	kw_40[23] = {
		{"acv_adaptive",0,0,1,1,0,0.,0.,3,N_mdm(utype,subMethod_SUBMETHOD_ACV_KL)},
		{"acv_independent_sampling",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_ACV_IS)},
		{"acv_is",0,0,1,1,0,0.,0.,-1,N_mdm(utype,subMethod_SUBMETHOD_ACV_IS)},
		{"acv_kl",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_ACV_KL)},
		{"acv_mf",0,0,1,1,0,0.,0.,1,N_mdm(utype,subMethod_SUBMETHOD_ACV_MF)},
		{"acv_multifidelity",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_ACV_MF)},
		{"convergence_tolerance",10,0,10,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,14,0,kw_33},
		{"export_sample_sequence",8,3,9,0,kw_35,0.,0.,0,N_mdm(true,exportSampleSeqFlag)},
		{"final_moments",8,3,13,0,kw_36},
		{"fixed_seed",8,0,7,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"initial_samples",5,0,2,0,0,0.,0.,5,N_mdm(szarray,pilotSamples)},
		{"max_function_evaluations",0x29,0,12,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,16,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"nip",8,0,5,0,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_NIP)},
		{"pilot_samples",13,0,2,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"rng",8,2,15,0,kw_37},
		{"sample_type",8,2,8,0,kw_38},
		{"seed_sequence",13,0,6,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"solution_mode",8,3,3,0,kw_39},
		{"sqp",8,0,5,0,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_SQP)},
		{"truth_fixed_by_pilot",8,0,4,0,0,0.,0.,0,N_mdm(true,truthPilotConstraint)}
		},
	kw_41[4] = {
		{"constant_liar",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_cl)},
		{"distance_penalty",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_distance)},
		{"naive",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_naive)},
		{"topology",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_topology)}
		},
	kw_42[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_43[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_44[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_43,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_45[3] = {
		{"distance",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_distance)},
		{"gradient",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_gradient)},
		{"predicted_variance",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_predicted_variance)}
		},
	kw_46[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_47[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_48[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_47,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_49[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_50[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_51[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_50}
		},
	kw_52[2] = {
		{"compute",8,3,2,0,kw_51},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_53[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_54[19] = {
		{"batch_selection",8,4,5,0,kw_41},
		{"distribution",8,2,14,0,kw_42},
		{"export_approx_points_file",11,3,8,0,kw_44,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,8,0,kw_44,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"fitness_metric",8,3,4,0,kw_45},
		{"gen_reliability_levels",14,1,13,0,kw_46,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,7,0,kw_48,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,7,0,kw_48,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,9,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,16,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,12,0,kw_49,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,6,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"response_levels",14,2,11,0,kw_52,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,15,0,kw_53},
		{"samples",1,0,1,0,0,0.,0.,-8,N_mdm(int,numSamples)},
		{"samples_on_emulator",9,0,3,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_55[7] = {
		{"merit1",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit1)},
		{"merit1_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit1_smooth)},
		{"merit2",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2)},
		{"merit2_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2_smooth)},
		{"merit2_squared",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2_squared)},
		{"merit_max",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit_max)},
		{"merit_max_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit_max_smooth)}
		},
	kw_56[2] = {
		{"blocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_BLOCKING_SYNCHRONIZATION)},
		{"nonblocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_NONBLOCKING_SYNCHRONIZATION)}
		},
	kw_57[13] = {
		{"constraint_penalty",10,0,7,0,0,0.,0.,0,N_mdm(Real,constrPenalty)},
		{"constraint_tolerance",10,0,9,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"contraction_factor",10,0,2,0,0,0.,0.,0,N_mdm(Real,contractStepLength)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initStepLength)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"merit_function",8,7,6,0,kw_55},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"smoothing_factor",10,0,8,0,0,0.,0.,0,N_mdm(Real,smoothFactor)},
		{"solution_accuracy",2,0,4,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,4,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"synchronization",8,2,5,0,kw_56},
		{"variable_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_58[1] = {
		{"hyperprior_betas",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,hyperPriorBetas)}
		},
	kw_59[5] = {
		{"both",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_BOTH)},
		{"hyperprior_alphas",14,1,2,0,kw_58,0.,0.,0,N_mdm(RealDL,hyperPriorAlphas)},
		{"one",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_ONE)},
		{"per_experiment",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_PER_EXPER)},
		{"per_response",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_PER_RESP)}
		},
	kw_60[1] = {
		{"confidence_intervals",8,0,1,0,0,0.,0.,0,N_mdm(true,chainDiagnosticsCI)}
		},
	kw_61[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_62[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_61,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_63[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_62,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_62,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_64[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_65[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_66[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_67[3] = {
		{"adapted",8,2,1,1,kw_66,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_68[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_69[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_70[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_71[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_72[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_73[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_74[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_68,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_68,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_69,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_70,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_71,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_71,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_70,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_72,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_73,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_73,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_75[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_76[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_77[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_76,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_78[6] = {
		{"basis_type",8,3,2,0,kw_67},
		{"collocation_ratio",10,21,3,1,kw_74,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_75,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_77,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_77,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_79[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_80[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_79,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_81[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_80,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_80,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_82[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_83[2] = {
		{"dimension_adaptive",8,3,1,1,kw_82},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_84[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_85[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_86[16] = {
		{"allocation_control",8,1,3,0,kw_64},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_65},
		{"expansion_order_sequence",13,6,5,1,kw_78,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_81,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_81,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_81,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_83,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_84,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_85,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_87[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_88[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_89[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_90[3] = {
		{"dimension_adaptive",8,2,1,1,kw_89},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_91[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_92[2] = {
		{"dimension_adaptive",8,2,1,1,kw_91},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_93[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_94[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_95[13] = {
		{"allocation_control",8,1,3,0,kw_87},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_88},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_90,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_92,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_93,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_94,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_96[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_97[2] = {
		{"estimator_variance",8,1,1,1,kw_96,0.,0.,0,N_mdm(type,multilevAllocControl_ESTIMATOR_VARIANCE)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_RIP_SAMPLING)}
		},
	kw_98[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_99[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_100[3] = {
		{"adapted",8,2,1,1,kw_99,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_101[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_102[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_103[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_104[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_105[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_106[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_107[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_101,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_101,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_102,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_103,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_104,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_104,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_103,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_105,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_106,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_106,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_108[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_109[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_110[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_109,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_111[6] = {
		{"basis_type",8,3,2,0,kw_100},
		{"collocation_ratio",10,21,3,1,kw_107,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_108,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_110,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_110,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_112[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_113[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_112,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_114[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_113,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_113,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_115[12] = {
		{"allocation_control",8,2,1,0,kw_97},
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,2,0,kw_98},
		{"expansion_order_sequence",13,6,3,1,kw_111,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,3,1,kw_114,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,3,1,kw_114,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,3,1,kw_114,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_116[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_117[3] = {
		{"adapted",8,2,1,1,kw_116,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_118[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_119[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_120[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_121[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_122[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_123[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_124[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_118,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_118,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_119,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_120,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_121,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_121,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_120,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_122,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_123,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_123,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_125[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_126[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_127[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_128[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_129[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_130[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_131[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_125,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_125,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_126,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_127,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_128,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_128,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_127,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_129,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_130,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_130,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_132[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_133[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_134[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_133,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_135[8] = {
		{"basis_type",8,3,2,0,kw_117},
		{"collocation_points",9,19,3,1,kw_124,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,19,3,1,kw_131,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,2,3,1,kw_132,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_134,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_134,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_136[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_137[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_136,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_138[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_137,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_137,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_139[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_140[2] = {
		{"dimension_adaptive",8,3,1,1,kw_139},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_141[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_142[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_143[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_135,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_138,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_138,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_138,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_140,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_141,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_142,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_144[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_145[3] = {
		{"dimension_adaptive",8,2,1,1,kw_144},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_146[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_147[2] = {
		{"dimension_adaptive",8,2,1,1,kw_146},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_148[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_149[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_150[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_145,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_147,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_148,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_149,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_151[7] = {
		{"gaussian_process",8,6,1,1,kw_63},
		{"kriging",0,6,1,1,kw_63,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_86,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,13,1,1,kw_95,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,12,1,1,kw_115,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_143,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_150,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
		},
	kw_152[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_153[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_152,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_154[11] = {
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"chains",0x29,0,3,0,0,3.,0.,0,N_mdm(int,numChains)},
		{"crossover_chain_pairs",0x29,0,5,0,0,0.,0.,0,N_mdm(int,crossoverChainPairs)},
		{"emulator",8,7,8,0,kw_151},
		{"export_chain_points_file",11,3,10,0,kw_153,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"gr_threshold",0x1a,0,6,0,0,0.,0.,0,N_mdm(Real,grThreshold)},
		{"jump_step",0x29,0,7,0,0,0.,0.,0,N_mdm(int,jumpStep)},
		{"num_cr",0x29,0,4,0,0,1.,0.,0,N_mdm(int,numCR)},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"standardized_space",8,0,9,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_155[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_IFACE_ID)}
		},
	kw_156[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_155,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_NONE)}
		},
	kw_157[7] = {
		{"batch_size",0x29,0,4,0,0,1.,0.,0,N_mdm(int,batchSize)},
		{"import_candidate_points_file",11,3,5,0,kw_156,0.,0.,0,N_mdm(str,importCandPtsFile)},
		{"initial_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,numSamples)},
		{"ksg2",8,0,6,0,0,0.,0.,0,N_mdm(true,mutualInfoKSG2)},
		{"max_hifi_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(int,maxHifiEvals)},
		{"num_candidates",0x19,0,2,2,0,0.,0.,0,N_mdm(sizet,numCandidates)},
		{"samples",1,0,1,1,0,0.,0.,-4,N_mdm(int,numSamples)}
		},
	kw_158[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_159[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_158,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_160[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_161[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_160,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_162[1] = {
		{"update_period",9,0,1,0,0,0.,0.,0,N_mdm(int,proposalCovUpdatePeriod)}
		},
	kw_163[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_164[1] = {
		{"multiplier",0x1a,0,1,0,0,0.,0.,0,N_mdm(Real,priorPropCovMult)}
		},
	kw_165[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_166[4] = {
		{"derivatives",8,1,1,1,kw_162,0.,0.,0,N_mdm(lit,proposalCovType_derivatives)},
		{"filename",11,2,1,1,kw_163,0.,0.,0,N_mdm(str,proposalCovFile)},
		{"prior",8,1,1,1,kw_164,0.,0.,0,N_mdm(lit,proposalCovType_prior)},
		{"values",14,2,1,1,kw_165,0.,0.,0,N_mdm(RealDL,proposalCovData)}
		},
	kw_167[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_168[17] = {
		{"adaptive_metropolis",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_adaptive_metropolis)},
		{"build_samples",9,0,4,2,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"delayed_rejection",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_delayed_rejection)},
		{"dram",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_dram)},
		{"export_chain_points_file",11,3,9,0,kw_159,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"gpmsa_normalize",8,0,8,0,0,0.,0.,0,N_mdm(true,gpmsaNormalize)},
		{"import_build_points_file",11,3,5,0,kw_161,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,3,5,0,kw_161,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"logit_transform",8,0,7,0,0,0.,0.,0,N_mdm(true,logitTransform)},
		{"metropolis_hastings",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_metropolis_hastings)},
		{"options_file",11,0,12,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"proposal_covariance",8,4,11,0,kw_166,0.,0.,0,N_mdm(lit,proposalCovType_user)},
		{"rng",8,2,3,0,kw_167},
		{"samples",1,0,1,1,0,0.,0.,-12,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"standardized_space",8,0,6,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_169[1] = {
		{"trend_order",0x29,0,1,0,0,0.,0.,0,N_mdm(shint,polynomialOrder)}
		},
	kw_170[1] = {
		{"basis_order",0x29,0,1,0,0,0.,0.,0,N_mdm(shint,polynomialOrder)}
		},
	kw_171[3] = {
		{"gaussian_process",8,1,1,1,kw_169,0.,0.,0,N_mdm(lit,modelDiscrepancyType_global_kriging)},
		{"kriging",0,1,1,1,kw_169,0.,0.,-1,N_mdm(lit,modelDiscrepancyType_global_kriging)},
		{"polynomial",8,1,1,1,kw_170,0.,0.,0,N_mdm(lit,modelDiscrepancyType_global_polynomial)}
		},
	kw_172[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_IFACE_ID)}
		},
	kw_173[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_172,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_NONE)}
		},
	kw_174[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_IFACE_ID)}
		},
	kw_175[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_174,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_NONE)}
		},
	kw_176[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_IFACE_ID)}
		},
	kw_177[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_176,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_NONE)}
		},
	kw_178[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_IFACE_ID)}
		},
	kw_179[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_178,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_NONE)}
		},
	kw_180[7] = {
		{"discrepancy_type",8,3,1,0,kw_171},
		{"export_corrected_model_file",11,3,6,0,kw_173,0.,0.,0,N_mdm(str,exportCorrModelFile)},
		{"export_corrected_variance_file",11,3,7,0,kw_175,0.,0.,0,N_mdm(str,exportCorrVarFile)},
		{"export_discrepancy_file",11,3,5,0,kw_177,0.,0.,0,N_mdm(str,exportDiscrepFile)},
		{"import_prediction_configs",11,3,4,0,kw_179,0.,0.,0,N_mdm(str,importPredConfigs)},
		{"num_prediction_configs",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numPredConfigs)},
		{"prediction_configs",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,predictionConfigList)}
		},
	kw_181[3] = {
		{"evidence_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,evidenceSamples)},
		{"laplace_approx",8,0,3,0,0,0.,0.,0,N_mdm(true,modelEvidLaplace)},
		{"mc_approx",8,0,1,0,0,0.,0.,0,N_mdm(true,modelEvidMC)}
		},
	kw_182[1] = {
		{"update_period",9,0,1,0,0,0.,0.,0,N_mdm(int,proposalCovUpdatePeriod)}
		},
	kw_183[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_184[1] = {
		{"multiplier",0x1a,0,1,0,0,0.,0.,0,N_mdm(Real,priorPropCovMult)}
		},
	kw_185[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_186[4] = {
		{"derivatives",8,1,1,1,kw_182,0.,0.,0,N_mdm(lit,proposalCovType_derivatives)},
		{"filename",11,2,1,1,kw_183,0.,0.,0,N_mdm(str,proposalCovFile)},
		{"prior",8,1,1,1,kw_184,0.,0.,0,N_mdm(lit,proposalCovType_prior)},
		{"values",14,2,1,1,kw_185,0.,0.,0,N_mdm(RealDL,proposalCovData)}
		},
	kw_187[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_188[9] = {
		{"adaptive_metropolis",8,0,4,0,0,0.,0.,0,N_mdm(lit,mcmcType_adaptive_metropolis)},
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"delayed_rejection",8,0,4,0,0,0.,0.,0,N_mdm(lit,mcmcType_delayed_rejection)},
		{"dram",8,0,4,0,0,0.,0.,0,N_mdm(lit,mcmcType_dram)},
		{"metropolis_hastings",8,0,4,0,0,0.,0.,0,N_mdm(lit,mcmcType_metropolis_hastings)},
		{"proposal_covariance",8,4,5,0,kw_186,0.,0.,0,N_mdm(lit,proposalCovType_user)},
		{"rng",8,2,3,0,kw_187},
		{"samples",1,0,1,1,0,0.,0.,-6,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_189[1] = {
		{"ksg2",8,0,1,0,0,0.,0.,0,N_mdm(true,mutualInfoKSG2)}
		},
	kw_190[3] = {
		{"kde",8,0,3,0,0,0.,0.,0,N_mdm(true,posteriorStatsKDE)},
		{"kl_divergence",8,0,1,0,0,0.,0.,0,N_mdm(true,posteriorStatsKL)},
		{"mutual_info",8,1,2,0,kw_189,0.,0.,0,N_mdm(true,posteriorStatsMutual)}
		},
	kw_191[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_192[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_193[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_192,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_194[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_193,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_193,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_195[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_196[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_197[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_198[3] = {
		{"adapted",8,2,1,1,kw_197,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_199[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_200[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_201[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_202[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_203[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_204[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_205[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_199,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_199,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_200,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_201,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_202,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_202,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_201,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_203,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_204,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_204,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_206[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_207[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_208[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_207,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_209[6] = {
		{"basis_type",8,3,2,0,kw_198},
		{"collocation_ratio",10,21,3,1,kw_205,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_206,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_208,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_208,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_210[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_211[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_210,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_212[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_211,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_211,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_213[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_214[2] = {
		{"dimension_adaptive",8,3,1,1,kw_213},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_215[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_216[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_217[16] = {
		{"allocation_control",8,1,3,0,kw_195},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_196},
		{"expansion_order_sequence",13,6,5,1,kw_209,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_212,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_212,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_212,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_214,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_215,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_216,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_218[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_219[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_220[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_221[3] = {
		{"dimension_adaptive",8,2,1,1,kw_220},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_222[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_223[2] = {
		{"dimension_adaptive",8,2,1,1,kw_222},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_224[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_225[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_226[13] = {
		{"allocation_control",8,1,3,0,kw_218},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_219},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_221,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_223,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_224,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_225,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_227[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_228[2] = {
		{"estimator_variance",8,1,1,1,kw_227,0.,0.,0,N_mdm(type,multilevAllocControl_ESTIMATOR_VARIANCE)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_RIP_SAMPLING)}
		},
	kw_229[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_230[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_231[3] = {
		{"adapted",8,2,1,1,kw_230,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_232[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_233[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_234[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_235[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_236[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_237[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_238[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_232,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_232,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_233,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_234,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_235,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_235,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_234,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_236,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_237,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_237,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_239[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_240[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_241[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_240,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_242[6] = {
		{"basis_type",8,3,2,0,kw_231},
		{"collocation_ratio",10,21,3,1,kw_238,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_239,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_241,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_241,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_243[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_244[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_243,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_245[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_244,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_244,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_246[12] = {
		{"allocation_control",8,2,1,0,kw_228},
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,2,0,kw_229},
		{"expansion_order_sequence",13,6,3,1,kw_242,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,3,1,kw_245,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,3,1,kw_245,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,3,1,kw_245,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_247[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_248[3] = {
		{"adapted",8,2,1,1,kw_247,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_249[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_250[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_251[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_252[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_253[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_254[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_255[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_249,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_249,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_250,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_251,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_252,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_252,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_251,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_253,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_254,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_254,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_256[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_257[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_258[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_259[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_260[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_261[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_262[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_256,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_256,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_257,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_258,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_259,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_259,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_258,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_260,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_261,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_261,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_263[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_264[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_265[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_264,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_266[8] = {
		{"basis_type",8,3,2,0,kw_248},
		{"collocation_points",9,19,3,1,kw_255,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,19,3,1,kw_262,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,2,3,1,kw_263,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_265,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_265,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_267[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_268[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_267,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_269[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_268,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_268,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_270[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_271[2] = {
		{"dimension_adaptive",8,3,1,1,kw_270},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_272[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_273[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_274[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_266,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_269,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_269,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_269,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_271,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_272,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_273,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_275[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_276[3] = {
		{"dimension_adaptive",8,2,1,1,kw_275},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_277[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_278[2] = {
		{"dimension_adaptive",8,2,1,1,kw_277},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_279[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_280[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_281[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_276,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_278,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_279,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_280,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_282[7] = {
		{"gaussian_process",8,6,1,1,kw_194},
		{"kriging",0,6,1,1,kw_194,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_217,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,13,1,1,kw_226,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,12,1,1,kw_246,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_274,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_281,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
		},
	kw_283[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_284[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_283,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_285[3] = {
		{"nip",8,0,1,1,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_NIP)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_NONE)},
		{"sqp",8,0,1,1,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_SQP)}
		},
	kw_286[1] = {
		{"update_period",9,0,1,0,0,0.,0.,0,N_mdm(int,proposalCovUpdatePeriod)}
		},
	kw_287[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_288[1] = {
		{"multiplier",0x1a,0,1,0,0,0.,0.,0,N_mdm(Real,priorPropCovMult)}
		},
	kw_289[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_290[4] = {
		{"derivatives",8,1,1,1,kw_286,0.,0.,0,N_mdm(lit,proposalCovType_derivatives)},
		{"filename",11,2,1,1,kw_287,0.,0.,0,N_mdm(str,proposalCovFile)},
		{"prior",8,1,1,1,kw_288,0.,0.,0,N_mdm(lit,proposalCovType_prior)},
		{"values",14,2,1,1,kw_289,0.,0.,0,N_mdm(RealDL,proposalCovData)}
		},
	kw_291[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_292[16] = {
		{"adaptive_metropolis",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_adaptive_metropolis)},
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"delayed_rejection",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_delayed_rejection)},
		{"dram",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_dram)},
		{"emulator",8,7,4,0,kw_282},
		{"export_chain_points_file",11,3,7,0,kw_284,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"logit_transform",8,0,6,0,0,0.,0.,0,N_mdm(true,logitTransform)},
		{"metropolis_hastings",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_metropolis_hastings)},
		{"multilevel",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_multilevel)},
		{"options_file",11,0,11,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"pre_solve",8,3,9,0,kw_285},
		{"proposal_covariance",8,4,10,0,kw_290,0.,0.,0,N_mdm(lit,proposalCovType_user)},
		{"rng",8,2,3,0,kw_291},
		{"samples",1,0,1,1,0,0.,0.,-12,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"standardized_space",8,0,5,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_293[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,dataDistCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,dataDistCovInputType_matrix)}
		},
	kw_294[2] = {
		{"covariance",14,2,2,2,kw_293,0.,0.,0,N_mdm(RealDL,dataDistCovariance)},
		{"means",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,dataDistMeans)}
		},
	kw_295[2] = {
		{"gaussian",8,2,1,1,kw_294},
		{"obs_data_filename",11,0,1,1,0,0.,0.,0,N_mdm(str,dataDistFile)}
		},
	kw_296[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_297[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_296,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_298[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_297,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_297,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_299[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_300[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_301[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_302[3] = {
		{"adapted",8,2,1,1,kw_301,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_303[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_304[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_305[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_306[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_307[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_308[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_309[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_303,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_303,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_304,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_305,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_306,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_306,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_305,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_307,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_308,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_308,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_310[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_311[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_312[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_311,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_313[6] = {
		{"basis_type",8,3,2,0,kw_302},
		{"collocation_ratio",10,21,3,1,kw_309,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_310,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_312,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_312,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_314[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_315[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_314,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_316[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_315,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_315,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_317[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_318[2] = {
		{"dimension_adaptive",8,3,1,1,kw_317},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_319[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_320[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_321[16] = {
		{"allocation_control",8,1,3,0,kw_299},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_300},
		{"expansion_order_sequence",13,6,5,1,kw_313,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_316,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_316,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_316,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_318,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_319,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_320,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_322[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_323[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_324[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_325[3] = {
		{"dimension_adaptive",8,2,1,1,kw_324},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_326[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_327[2] = {
		{"dimension_adaptive",8,2,1,1,kw_326},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_328[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_329[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_330[13] = {
		{"allocation_control",8,1,3,0,kw_322},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,4,0,kw_323},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_325,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_327,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_328,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_329,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_331[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_332[2] = {
		{"estimator_variance",8,1,1,1,kw_331,0.,0.,0,N_mdm(type,multilevAllocControl_ESTIMATOR_VARIANCE)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_RIP_SAMPLING)}
		},
	kw_333[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_334[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_335[3] = {
		{"adapted",8,2,1,1,kw_334,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_336[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_337[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_338[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_339[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_340[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_341[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_342[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_336,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_336,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_337,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_338,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_339,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_339,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_338,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_340,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_341,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_341,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_343[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_344[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_345[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_344,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_346[6] = {
		{"basis_type",8,3,2,0,kw_335},
		{"collocation_ratio",10,21,3,1,kw_342,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_343,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_345,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_345,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_347[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_348[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_347,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_349[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_348,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_348,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_350[12] = {
		{"allocation_control",8,2,1,0,kw_332},
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,2,0,kw_333},
		{"expansion_order_sequence",13,6,3,1,kw_346,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,3,1,kw_349,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,3,1,kw_349,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,3,1,kw_349,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_351[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_352[3] = {
		{"adapted",8,2,1,1,kw_351,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_353[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_354[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_355[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_356[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_357[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_358[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_359[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_353,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_353,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_354,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_355,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_356,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_356,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_355,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_357,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_358,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_358,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_360[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_361[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_362[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_363[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_364[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_365[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_366[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_360,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_360,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_361,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_362,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_363,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_363,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_362,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_364,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_365,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_365,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_367[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_368[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_369[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_368,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_370[8] = {
		{"basis_type",8,3,2,0,kw_352},
		{"collocation_points",9,19,3,1,kw_359,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,19,3,1,kw_366,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,2,3,1,kw_367,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_369,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_369,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_371[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_372[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_371,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_373[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_372,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_372,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_374[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_375[2] = {
		{"dimension_adaptive",8,3,1,1,kw_374},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_376[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_377[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_378[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_370,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_373,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_373,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_373,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_375,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_376,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_377,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_379[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_380[3] = {
		{"dimension_adaptive",8,2,1,1,kw_379},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_381[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_382[2] = {
		{"dimension_adaptive",8,2,1,1,kw_381},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_383[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_384[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_385[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_380,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_382,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_383,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_384,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_386[7] = {
		{"gaussian_process",8,6,1,1,kw_298},
		{"kriging",0,6,1,1,kw_298,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_321,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,13,1,1,kw_330,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,12,1,1,kw_350,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_378,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_385,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
		},
	kw_387[1] = {
		{"posterior_density_export_filename",11,0,1,0,0,0.,0.,0,N_mdm(str,posteriorDensityExportFilename)}
		},
	kw_388[1] = {
		{"posterior_samples_export_filename",11,0,1,0,0,0.,0.,0,N_mdm(str,posteriorSamplesExportFilename)}
		},
	kw_389[8] = {
		{"data_distribution",8,2,5,2,kw_295},
		{"emulator",8,7,3,0,kw_386},
		{"evaluate_posterior_density",8,1,8,0,kw_387,0.,0.,0,N_mdm(true,evaluatePosteriorDensity)},
		{"generate_posterior_samples",8,1,7,0,kw_388,0.,0.,0,N_mdm(true,generatePosteriorSamples)},
		{"posterior_samples_import_filename",11,0,6,0,0,0.,0.,0,N_mdm(str,posteriorSamplesImportFilename)},
		{"pushforward_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,numPushforwardSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"standardized_space",8,0,4,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_390[18] = {
		{"burn_in_samples",9,0,4,0,0,0.,0.,0,N_mdm(int,burnInSamples)},
		{"calibrate_error_multipliers",8,5,3,0,kw_59},
		{"chain_diagnostics",8,1,6,0,kw_60,0.,0.,0,N_mdm(true,chainDiagnostics)},
		{"convergence_tolerance",10,0,11,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"dream",8,11,1,1,kw_154,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_DREAM)},
		{"experimental_design",8,7,2,0,kw_157,0.,0.,0,N_mdm(true,adaptExpDesign)},
		{"gpmsa",8,17,1,1,kw_168,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_GPMSA)},
		{"max_iterations",0x29,0,12,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_discrepancy",8,7,8,0,kw_180,0.,0.,0,N_mdm(true,calModelDiscrepancy)},
		{"model_evidence",8,3,7,0,kw_181,0.,0.,0,N_mdm(true,modelEvidence)},
		{"model_pointer",11,0,13,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"muq",8,9,1,1,kw_188,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_MUQ)},
		{"posterior_stats",8,3,5,0,kw_190},
		{"probability_levels",14,1,10,0,kw_191,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"queso",8,16,1,1,kw_292,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_QUESO)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"sub_sampling_period",9,0,9,0,0,0.,0.,0,N_mdm(int,subSamplingPeriod)},
		{"wasabi",8,8,1,1,kw_389,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_WASABI)}
		},
	kw_391[1] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,modelPointer)}
		},
	kw_392[3] = {
		{"method_name",11,1,1,1,kw_391,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"scaling",8,0,2,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_393[4] = {
		{"deltas_per_variable",5,0,2,2,0,0.,0.,3,N_mdm(ivec,stepsPerVariable)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"step_vector",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,stepVector)},
		{"steps_per_variable",13,0,2,2,0,0.,0.,0,N_mdm(ivec,stepsPerVariable)}
		},
	kw_394[11] = {
		{"beta_solver_name",11,0,1,1,0,0.,0.,0,N_mdm(str,betaSolverName)},
		{"convergence_tolerance",10,0,7,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,5,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,4,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,2,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,2,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_395[12] = {
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,6,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,10,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,4,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,5,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,3,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,3,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_396[2] = {
		{"all_dimensions",8,0,1,1,0,0.,0.,0,N_mdm(lit,boxDivision_all_dimensions)},
		{"major_dimension",8,0,1,1,0,0.,0.,0,N_mdm(lit,boxDivision_major_dimension)}
		},
	kw_397[16] = {
		{"constraint_penalty",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"convergence_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"division",8,2,1,0,kw_396},
		{"global_balance_parameter",10,0,2,0,0,0.,0.,0,N_mdm(Real,globalBalanceParam)},
		{"local_balance_parameter",10,0,3,0,0,0.,0.,0,N_mdm(Real,localBalanceParam)},
		{"max_boxsize_limit",10,0,4,0,0,0.,0.,0,N_mdm(Real,maxBoxSize)},
		{"max_function_evaluations",0x29,0,13,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"min_boxsize_limit",10,0,5,0,0,0.,0.,0,N_mdm(Real,minBoxSize)},
		{"misc_options",15,0,10,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,8,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,9,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,7,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,7,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_398[3] = {
		{"blend",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_blend)},
		{"two_point",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_two_point)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_uniform)}
		},
	kw_399[2] = {
		{"linear_rank",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_linear_rank)},
		{"merit_function",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_proportional)}
		},
	kw_400[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_401[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_402[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_403[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_404[5] = {
		{"non_adaptive",8,0,2,0,0,0.,0.,0,N_mdm(false,mutationAdaptive)},
		{"offset_cauchy",8,2,1,1,kw_401,0.,0.,0,N_mdm(lit,mutationType_offset_cauchy)},
		{"offset_normal",8,2,1,1,kw_402,0.,0.,0,N_mdm(lit,mutationType_offset_normal)},
		{"offset_uniform",8,2,1,1,kw_403,0.,0.,0,N_mdm(lit,mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_405[4] = {
		{"chc",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_chc)},
		{"elitist",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_elitist)},
		{"new_solutions_generated",9,0,2,0,0,0.,0.,0,N_mdm(int,newSolnsGenerated)},
		{"random",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_random)}
		},
	kw_406[19] = {
		{"constraint_penalty",10,0,9,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"convergence_tolerance",10,0,15,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"crossover_rate",10,0,5,0,0,0.,0.,0,N_mdm(Real,crossoverRate)},
		{"crossover_type",8,3,6,0,kw_398},
		{"fitness_type",8,2,3,0,kw_399},
		{"initialization_type",8,3,2,0,kw_400},
		{"max_function_evaluations",0x29,0,16,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,14,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,13,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,18,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_rate",10,0,7,0,0,0.,0.,0,N_mdm(Real,mutationRate)},
		{"mutation_type",8,5,8,0,kw_404},
		{"population_size",0x19,0,1,0,0,0.,0.,0,N_mdm(int,populationSize)},
		{"replacement_type",8,4,4,0,kw_405},
		{"scaling",8,0,17,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,11,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,12,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,10,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,10,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_407[3] = {
		{"adaptive_pattern",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_adaptive)},
		{"basic_pattern",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_simple)},
		{"multi_step",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_multi_step)}
		},
	kw_408[2] = {
		{"coordinate",8,0,1,1,0,0.,0.,0,N_mdm(lit,patternBasis_coordinate)},
		{"simplex",8,0,1,1,0,0.,0.,0,N_mdm(lit,patternBasis_simplex)}
		},
	kw_409[2] = {
		{"blocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_BLOCKING_SYNCHRONIZATION)},
		{"nonblocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_NONBLOCKING_SYNCHRONIZATION)}
		},
	kw_410[22] = {
		{"constant_penalty",8,0,1,0,0,0.,0.,0,N_mdm(true,constantPenalty)},
		{"constraint_penalty",10,0,10,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"contraction_factor",10,0,9,0,0,0.,0.,0,N_mdm(Real,contractFactor)},
		{"convergence_tolerance",10,0,18,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"expand_after_success",9,0,3,0,0,0.,0.,0,N_mdm(int,expandAfterSuccess)},
		{"exploratory_moves",8,3,7,0,kw_407},
		{"initial_delta",10,0,11,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,19,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,17,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,16,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,21,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"no_expansion",8,0,2,0,0,0.,0.,0,N_mdm(false,expansionFlag)},
		{"pattern_basis",8,2,4,0,kw_408},
		{"scaling",8,0,20,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,14,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,15,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,13,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,13,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"stochastic",8,0,5,0,0,0.,0.,0,N_mdm(true,randomizeOrderFlag)},
		{"synchronization",8,2,8,0,kw_409},
		{"total_pattern_size",9,0,6,0,0,0.,0.,0,N_mdm(int,totalPatternSize)},
		{"variable_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_411[18] = {
		{"constant_penalty",8,0,4,0,0,0.,0.,0,N_mdm(true,constantPenalty)},
		{"constraint_penalty",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"contract_after_failure",9,0,1,0,0,0.,0.,0,N_mdm(int,contractAfterFail)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,contractFactor)},
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"expand_after_success",9,0,3,0,0,0.,0.,0,N_mdm(int,expandAfterSuccess)},
		{"initial_delta",10,0,7,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,15,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,13,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,12,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,17,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"no_expansion",8,0,2,0,0,0.,0.,0,N_mdm(false,expansionFlag)},
		{"scaling",8,0,16,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,10,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"show_misc_options",8,0,11,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,9,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,9,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_412[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_413[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_414[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_415[15] = {
		{"box_behnken",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_BOX_BEHNKEN)},
		{"central_composite",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_CENTRAL_COMPOSITE)},
		{"fixed_seed",8,0,4,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"grid",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_GRID)},
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_LHS)},
		{"main_effects",8,0,5,0,0,0.,0.,0,N_mdm(true,mainEffectsFlag)},
		{"model_pointer",11,0,9,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"oa_lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_OA_LHS)},
		{"oas",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_OAS)},
		{"quality_metrics",8,0,6,0,0,0.,0.,0,N_mdm(true,volQualityFlag)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_RANDOM)},
		{"samples",9,0,2,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"symbols",9,0,8,0,0,0.,0.,0,N_mdm(int,numSymbols)},
		{"variance_based_decomp",8,1,7,0,kw_414,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_416[7] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"options_file",11,0,6,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"solution_accuracy",2,0,5,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,5,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_417[3] = {
		{"max_function_evaluations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,2,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_418[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_419[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_420[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_421[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_422[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_423[2] = {
		{"blocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_BLOCKING_SYNCHRONIZATION)},
		{"nonblocking",8,0,1,1,0,0.,0.,0,N_mdm(type,evalSynchronize_NONBLOCKING_SYNCHRONIZATION)}
		},
	kw_424[2] = {
		{"exploration",0x29,0,1,0,0,0.,0.,0,N_mdm(int,batchSizeExplore)},
		{"synchronization",8,2,2,0,kw_423}
		},
	kw_425[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_426[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_425,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_427[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_428[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_427}
		},
	kw_429[2] = {
		{"export_model",8,2,1,0,kw_428,0.,0.,0,N_mdm(true,exportSurrogate)},
		{"options_file",11,0,2,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)}
		},
	kw_430[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_431[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_430}
		},
	kw_432[1] = {
		{"export_model",8,2,1,0,kw_431,0.,0.,0,N_mdm(true,exportSurrogate)}
		},
	kw_433[3] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"experimental",8,2,1,1,kw_429,0.,0.,0,N_mdm(type,emulatorType_EXPGP_EMULATOR)},
		{"surfpack",8,1,1,1,kw_432,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_434[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_435[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_434,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_436[14] = {
		{"batch_size",0x29,2,3,0,kw_424,1.,0.,0,N_mdm(int,batchSize)},
		{"convergence_tolerance",10,0,5,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"export_approx_points_file",11,3,10,0,kw_426,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,10,0,kw_426,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,3,7,0,kw_433},
		{"import_build_points_file",11,4,9,0,kw_435,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,9,0,kw_435,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"kriging",0,3,7,0,kw_433,0.,0.,-4},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"use_derivatives",8,0,8,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"x_conv_tol",10,0,6,0,0,0.,0.,0,N_mdm(Real,xConvTol)}
		},
	kw_437[3] = {
		{"grid",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_grid)},
		{"halton",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_halton)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_random)}
		},
	kw_438[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_439[10] = {
		{"fixed_seed",8,0,3,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"latinize",8,0,4,0,0,0.,0.,0,N_mdm(true,latinizeFlag)},
		{"max_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"num_trials",9,0,8,0,0,0.,0.,0,N_mdm(int,numTrials)},
		{"quality_metrics",8,0,5,0,0,0.,0.,0,N_mdm(true,volQualityFlag)},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"trial_type",8,3,7,0,kw_437},
		{"variance_based_decomp",8,1,6,0,kw_438,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_440[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_441[12] = {
		{"fixed_sequence",8,0,6,0,0,0.,0.,0,N_mdm(true,fixedSequenceFlag)},
		{"halton",8,0,1,1,0,0.,0.,0,N_mdm(utype,methodName_FSU_HALTON)},
		{"hammersley",8,0,1,1,0,0.,0.,0,N_mdm(utype,methodName_FSU_HAMMERSLEY)},
		{"latinize",8,0,2,0,0,0.,0.,0,N_mdm(true,latinizeFlag)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"prime_base",13,0,9,0,0,0.,0.,0,N_mdm(ivec,primeBase)},
		{"quality_metrics",8,0,3,0,0,0.,0.,0,N_mdm(true,volQualityFlag)},
		{"samples",9,0,5,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sequence_leap",13,0,8,0,0,0.,0.,0,N_mdm(ivec,sequenceLeap)},
		{"sequence_start",13,0,7,0,0,0.,0.,0,N_mdm(ivec,sequenceStart)},
		{"variance_based_decomp",8,1,4,0,kw_440,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_442[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_443[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_444[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_443,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_445[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_446[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_447[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_448[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_447,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_449[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_450[1] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)}
		},
	kw_451[5] = {
		{"increment_max_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_ORDER_ADVANCEMENT)},
		{"increment_max_rank",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_RANK_ADVANCEMENT)},
		{"increment_max_rank_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_RANK_ORDER_ADVANCEMENT)},
		{"increment_start_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_START_ORDER_ADVANCEMENT)},
		{"increment_start_rank",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_START_RANK_ADVANCEMENT)}
		},
	kw_452[1] = {
		{"uniform",8,5,1,1,kw_451,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_453[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_454[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_455[1] = {
		{"l2_penalty",10,0,1,1,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)}
		},
	kw_456[2] = {
		{"ls",8,0,1,1,0,0.,0.,0,N_mdm(type,regressionType_FT_LS)},
		{"rls2",8,1,1,1,kw_455,0.,0.,0,N_mdm(type,regressionType_FT_RLS2)}
		},
	kw_457[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_458[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_459[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_458}
		},
	kw_460[2] = {
		{"compute",8,4,2,0,kw_459},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_461[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_462[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_463[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_464[47] = {
		{"adapt_order",8,0,15,0,0,0.,0.,0,N_mdm(true,adaptOrder)},
		{"adapt_rank",8,0,20,0,0,0.,0.,0,N_mdm(true,adaptRank)},
		{"arithmetic_tolerance",10,0,13,0,0,0.,0.,0,N_mdm(Real,statsRoundingTol)},
		{"collocation_points",9,0,11,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,0,11,1,0,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,35,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,33,0,kw_442},
		{"export_approx_points_file",11,3,37,0,kw_444,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,37,0,kw_444,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,28,0,kw_445},
		{"fixed_seed",8,0,39,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,35,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,32,0,kw_446,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,36,0,kw_448,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"kick_order",0x19,0,16,0,0,0.,0.,0,N_mdm(ushint,kickOrder)},
		{"kick_rank",0x19,0,21,0,0,0.,0.,0,N_mdm(sizet,kickRank)},
		{"max_cross_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(int,maxCrossIterations)},
		{"max_cv_order_candidates",0x29,0,18,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"max_cv_rank_candidates",0x29,0,23,0,0,0.,0.,0,N_mdm(sizet,maxCVRankCandidates)},
		{"max_order",0x29,0,17,0,0,0.,0.,0,N_mdm(ushint,maxOrder)},
		{"max_rank",0x29,0,22,0,0,0.,0.,0,N_mdm(sizet,maxRank)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"max_solver_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"metric_scale",8,2,4,0,kw_449},
		{"model_pointer",11,0,40,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"order",0x21,1,14,0,kw_450,0.,0.,17,N_mdm(ushint,startOrder)},
		{"p_refinement",8,1,1,0,kw_452,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"probability_levels",14,1,30,0,kw_453,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,27,0,kw_454},
		{"rank",0x21,0,19,0,0,0.,0.,14,N_mdm(sizet,startRank)},
		{"regression_type",8,2,5,0,kw_456},
		{"reliability_levels",14,1,31,0,kw_457,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,29,0,kw_460,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"response_scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"rng",8,2,26,0,kw_461},
		{"rounding_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,solverRoundingTol)},
		{"sample_refinement",0,4,27,0,kw_454,0.,0.,-8},
		{"sample_type",8,2,25,0,kw_462},
		{"samples",1,0,24,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,24,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,38,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"solver_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,solverTol)},
		{"start_order",0x29,1,14,0,kw_450,0.,0.,0,N_mdm(ushint,startOrder)},
		{"start_rank",0x29,0,19,0,0,0.,0.,0,N_mdm(sizet,startRank)},
		{"tensor_grid",8,0,10,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"variance_based_decomp",8,2,34,0,kw_463,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_465[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_466[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_467[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_466,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_468[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_469[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_470[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_469,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_471[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_472[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_473[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_472}
		},
	kw_474[2] = {
		{"compute",8,3,2,0,kw_473},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_475[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_476[15] = {
		{"build_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,10,0,kw_465},
		{"export_approx_points_file",11,3,5,0,kw_467,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,5,0,kw_467,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gen_reliability_levels",14,1,9,0,kw_468,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,4,0,kw_470,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_470,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,8,0,kw_471,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,7,0,kw_474,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_475},
		{"samples",1,0,1,0,0,0.,0.,-12,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,3,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_477[4] = {
		{"max_function_evaluations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,3,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_478[4] = {
		{"max_function_evaluations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,3,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_479[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_480[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_481[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_480,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_482[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_483[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_484[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_483,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_485[7] = {
		{"export_approx_points_file",11,3,4,0,kw_481,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_481,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,1,0,kw_482},
		{"import_build_points_file",11,4,3,0,kw_484,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_484,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,2,1,0,kw_482,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_486[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_487[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_488[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_489[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_488}
		},
	kw_490[2] = {
		{"compute",8,3,2,0,kw_489},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_491[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_492[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_493[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_492,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_494[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_495[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_496[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_495,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_497[7] = {
		{"export_approx_points_file",11,3,4,0,kw_493,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_493,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,1,0,kw_494},
		{"import_build_points_file",11,4,3,0,kw_496,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_496,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,2,1,0,kw_494,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_498[12] = {
		{"distribution",8,2,7,0,kw_479},
		{"ea",8,0,3,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EA)},
		{"ego",8,7,3,0,kw_485,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGO)},
		{"gen_reliability_levels",14,1,6,0,kw_486,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lhs",8,0,3,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_LHS)},
		{"model_pointer",11,0,9,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,5,0,kw_487,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,4,0,kw_490,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,8,0,kw_491},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sbo",8,7,3,0,kw_497,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SBO)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_499[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_500[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_499,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_501[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_502[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_501}
		},
	kw_503[2] = {
		{"export_model",8,2,1,0,kw_502,0.,0.,0,N_mdm(true,exportSurrogate)},
		{"options_file",11,0,2,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)}
		},
	kw_504[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_505[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_504}
		},
	kw_506[1] = {
		{"export_model",8,2,1,0,kw_505,0.,0.,0,N_mdm(true,exportSurrogate)}
		},
	kw_507[3] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"experimental",8,2,1,1,kw_503,0.,0.,0,N_mdm(type,emulatorType_EXPGP_EMULATOR)},
		{"surfpack",8,1,1,1,kw_506,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_508[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_509[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_508,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_510[7] = {
		{"export_approx_points_file",11,3,4,0,kw_500,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_500,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,3,1,0,kw_507},
		{"import_build_points_file",11,4,3,0,kw_509,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_509,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,3,1,0,kw_507,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_511[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_512[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_513[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_512,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_514[1] = {
		{"options_file",11,0,1,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)}
		},
	kw_515[3] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"experimental",8,1,1,1,kw_514,0.,0.,0,N_mdm(type,emulatorType_EXPGP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_516[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_517[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_516,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_518[7] = {
		{"export_approx_points_file",11,3,4,0,kw_513,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_513,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,3,1,0,kw_515},
		{"import_build_points_file",11,4,3,0,kw_517,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_517,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,3,1,0,kw_515,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_519[11] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"ea",8,0,6,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EA)},
		{"ego",8,7,6,0,kw_510,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGO)},
		{"lhs",8,0,6,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_LHS)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"rng",8,2,7,0,kw_511},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sbo",8,7,6,0,kw_518,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SBO)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_520[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_521[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_522[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_521}
		},
	kw_523[2] = {
		{"export_model",8,2,1,0,kw_522,0.,0.,0,N_mdm(true,exportSurrogate)},
		{"options_file",11,0,2,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)}
		},
	kw_524[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_525[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_524,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_526[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_527[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_528[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_527,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_529[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_530[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_531[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_530}
		},
	kw_532[2] = {
		{"compute",8,3,2,0,kw_531},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_533[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_534[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_535[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mdm(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_534}
		},
	kw_536[1] = {
		{"export_model",8,2,1,0,kw_535,0.,0.,0,N_mdm(true,exportSurrogate)}
		},
	kw_537[22] = {
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"dakota",8,0,3,0,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"distribution",8,2,12,0,kw_520},
		{"experimental",8,2,3,0,kw_523,0.,0.,0,N_mdm(type,emulatorType_EXPGP_EMULATOR)},
		{"export_approx_points_file",11,3,5,0,kw_525,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,5,0,kw_525,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gen_reliability_levels",14,1,11,0,kw_526,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,4,0,kw_528,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_528,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,13,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,10,0,kw_529,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,9,0,kw_532,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,8,0,kw_533},
		{"seed",0x19,0,7,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"surfpack",8,1,3,0,kw_536,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)},
		{"u_gaussian_process",8,0,2,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGRA_U)},
		{"u_kriging",0,0,2,1,0,0.,0.,-1,N_mdm(utype,subMethod_SUBMETHOD_EGRA_U)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"x_gaussian_process",8,0,2,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGRA_X)},
		{"x_kriging",0,0,2,1,0,0.,0.,-1,N_mdm(utype,subMethod_SUBMETHOD_EGRA_X)}
		},
	kw_538[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_539[1] = {
		{"model_pointer_list",15,0,1,0,0,0.,0.,0,N_mdm(strL,hybridModelPointers)}
		},
	kw_540[5] = {
		{"iterator_scheduling",8,2,3,0,kw_538},
		{"iterator_servers",0x19,0,2,0,0,0.,0.,0,N_mdm(int,iteratorServers)},
		{"method_name_list",15,1,1,1,kw_539,0.,0.,0,N_mdm(strL,hybridMethodNames)},
		{"method_pointer_list",15,0,1,1,0,0.,0.,0,N_mdm(strL,hybridMethodPointers)},
		{"processors_per_iterator",0x19,0,4,0,0,0.,0.,0,N_mdm(int,procsPerIterator)}
		},
	kw_541[1] = {
		{"global_model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,hybridGlobalModelPointer)}
		},
	kw_542[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_543[1] = {
		{"local_model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,hybridLocalModelPointer)}
		},
	kw_544[8] = {
		{"global_method_name",11,1,1,1,kw_541,0.,0.,0,N_mdm(str,hybridGlobalMethodName)},
		{"global_method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,hybridGlobalMethodPointer)},
		{"iterator_scheduling",8,2,5,0,kw_542},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(int,iteratorServers)},
		{"local_method_name",11,1,2,2,kw_543,0.,0.,0,N_mdm(str,hybridLocalMethodName)},
		{"local_method_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,hybridLocalMethodPointer)},
		{"local_search_probability",10,0,3,0,0,0.,0.,0,N_mdm(Real,hybridLSProb)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(int,procsPerIterator)}
		},
	kw_545[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_546[1] = {
		{"model_pointer_list",15,0,1,0,0,0.,0.,0,N_mdm(strL,hybridModelPointers)}
		},
	kw_547[5] = {
		{"iterator_scheduling",8,2,3,0,kw_545},
		{"iterator_servers",0x19,0,2,0,0,0.,0.,0,N_mdm(int,iteratorServers)},
		{"method_name_list",15,1,1,1,kw_546,0.,0.,0,N_mdm(strL,hybridMethodNames)},
		{"method_pointer_list",15,0,1,1,0,0.,0.,0,N_mdm(strL,hybridMethodPointers)},
		{"processors_per_iterator",0x19,0,4,0,0,0.,0.,0,N_mdm(int,procsPerIterator)}
		},
	kw_548[5] = {
		{"collaborative",8,5,1,1,kw_540,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_COLLABORATIVE)},
		{"coupled",0,8,1,1,kw_544,0.,0.,1,N_mdm(utype,subMethod_SUBMETHOD_EMBEDDED)},
		{"embedded",8,8,1,1,kw_544,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EMBEDDED)},
		{"sequential",8,5,1,1,kw_547,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SEQUENTIAL)},
		{"uncoupled",0,5,1,1,kw_547,0.,0.,-1,N_mdm(utype,subMethod_SUBMETHOD_SEQUENTIAL)}
		},
	kw_549[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_550[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_551[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_552[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_553[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_552}
		},
	kw_554[2] = {
		{"compute",8,3,2,0,kw_553},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_555[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_556[15] = {
		{"adapt_import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"convergence_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,10,0,kw_549},
		{"gen_reliability_levels",14,1,9,0,kw_550,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"initial_samples",1,0,1,0,0,0.,0.,8,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"mm_adapt_import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,8,0,kw_551,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,4,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"response_levels",14,2,7,0,kw_554,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_555},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_557[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_IFACE_ID)}
		},
	kw_558[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,pstudyFileActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_557,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_NONE)}
		},
	kw_559[3] = {
		{"import_points_file",11,4,1,1,kw_558,0.,0.,0,N_mdm(str,pstudyFilename)},
		{"list_of_points",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,listOfPoints)},
		{"model_pointer",11,0,2,0,0,0.,0.,0,N_mdm(str,modelPointer)}
		},
	kw_560[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_561[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_562[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_563[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_564[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_563}
		},
	kw_565[2] = {
		{"compute",8,3,2,0,kw_564},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_566[7] = {
		{"distribution",8,2,5,0,kw_560},
		{"gen_reliability_levels",14,1,4,0,kw_561,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"nip",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NIP)},
		{"probability_levels",14,1,3,0,kw_562,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,2,0,kw_565,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"sqp",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SQP)}
		},
	kw_567[4] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"nip",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NIP)},
		{"sqp",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SQP)}
		},
	kw_568[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_569[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_570[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_571[5] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_572[4] = {
		{"first_order",8,0,1,1,0,0.,0.,0,N_mdm(lit,reliabilityIntegration_first_order)},
		{"probability_refinement",8,5,2,0,kw_571},
		{"sample_refinement",0,5,2,0,kw_571,0.,0.,-1},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mdm(lit,reliabilityIntegration_second_order)}
		},
	kw_573[12] = {
		{"integration",8,4,3,0,kw_572},
		{"nip",8,0,2,0,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_NIP)},
		{"no_approx",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NO_APPROX)},
		{"sqp",8,0,2,0,0,0.,0.,0,N_mdm(utype,optSubProbSolver_SUBMETHOD_SQP)},
		{"u_multi_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_QMEA_U)},
		{"u_taylor_mean",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_AMV_U)},
		{"u_taylor_mpp",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_AMV_PLUS_U)},
		{"u_two_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_TANA_U)},
		{"x_multi_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_QMEA_X)},
		{"x_taylor_mean",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_AMV_X)},
		{"x_taylor_mpp",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_AMV_PLUS_X)},
		{"x_two_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_TANA_X)}
		},
	kw_574[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_575[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_576[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_577[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_576}
		},
	kw_578[2] = {
		{"compute",8,4,2,0,kw_577},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_579[10] = {
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,6,0,kw_568},
		{"final_moments",8,3,9,0,kw_569},
		{"gen_reliability_levels",14,1,5,0,kw_570,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mpp_search",8,12,1,0,kw_573},
		{"probability_levels",14,1,3,0,kw_574,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"reliability_levels",14,1,4,0,kw_575,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,2,0,kw_578,0.,0.,0,N_mdm(resplevs,responseLevels)}
		},
	kw_580[2] = {
		{"inform_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,useSurrogate_inform_search)},
		{"optimize",8,0,1,1,0,0.,0.,0,N_mdm(lit,useSurrogate_optimize)}
		},
	kw_581[14] = {
		{"display_all_evaluations",8,0,9,0,0,0.,0.,0,N_mdm(true,showAllEval)},
		{"display_format",11,0,6,0,0,0.,0.,0,N_mdm(str,displayFormat)},
		{"function_precision",10,0,3,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"history_file",11,0,5,0,0,0.,0.,0,N_mdm(str,historyFile)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initMeshSize)},
		{"max_function_evaluations",0x29,0,12,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,14,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"neighbor_order",0x19,0,8,0,0,0.,0.,0,N_mdm(int,neighborOrder)},
		{"scaling",8,0,13,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,4,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"use_surrogate",8,2,10,0,kw_580},
		{"variable_neighborhood_search",10,0,7,0,0,0.,0.,0,N_mdm(Real,vns)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,minMeshSize)}
		},
	kw_582[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_583[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_584[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_583,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_585[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_586[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_587[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_588[3] = {
		{"offline_pilot",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_OFFLINE_PILOT)},
		{"online_pilot",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_ONLINE_PILOT)},
		{"pilot_projection",8,0,1,1,0,0.,0.,0,N_mdm(type,ensembleSampSolnMode_PILOT_PROJECTION)}
		},
	kw_589[14] = {
		{"convergence_tolerance",10,0,7,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,11,0,kw_582},
		{"export_sample_sequence",8,3,6,0,kw_584,0.,0.,0,N_mdm(true,exportSampleSeqFlag)},
		{"final_moments",8,3,10,0,kw_585},
		{"fixed_seed",8,0,2,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"initial_samples",5,0,3,0,0,0.,0.,4,N_mdm(szarray,pilotSamples)},
		{"max_function_evaluations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,13,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"pilot_samples",13,0,3,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"rng",8,2,12,0,kw_586},
		{"sample_type",8,2,5,0,kw_587},
		{"seed_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"solution_mode",8,3,4,0,kw_588}
		},
	kw_590[2] = {
		{"optimization",8,0,2,0,0,0.,0.,0,N_mdm(true,useTargetVarianceOptimizationFlag)},
		{"scalarization_response_mapping",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,scalarizationRespCoeffs)}
		},
	kw_591[1] = {
		{"optimization",8,0,1,0,0,0.,0.,0,N_mdm(true,useTargetVarianceOptimizationFlag)}
		},
	kw_592[1] = {
		{"optimization",8,0,1,0,0,0.,0.,0,N_mdm(true,useTargetVarianceOptimizationFlag)}
		},
	kw_593[4] = {
		{"mean",8,0,1,1,0,0.,0.,0,N_mdm(type,allocationTarget_TARGET_MEAN)},
		{"scalarization",8,2,1,1,kw_590,0.,0.,0,N_mdm(type,allocationTarget_TARGET_SCALARIZATION)},
		{"standard_deviation",8,1,1,1,kw_591,0.,0.,0,N_mdm(type,allocationTarget_TARGET_SIGMA)},
		{"variance",8,1,1,1,kw_592,0.,0.,0,N_mdm(type,allocationTarget_TARGET_VARIANCE)}
		},
	kw_594[2] = {
		{"cost_constraint",8,0,1,1,0,0.,0.,0,N_mdm(type,convergenceToleranceTarget_CONVERGENCE_TOLERANCE_TARGET_COST_CONSTRAINT)},
		{"variance_constraint",8,0,1,1,0,0.,0.,0,N_mdm(type,convergenceToleranceTarget_CONVERGENCE_TOLERANCE_TARGET_VARIANCE_CONSTRAINT)}
		},
	kw_595[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(type,convergenceToleranceType_CONVERGENCE_TOLERANCE_TYPE_ABSOLUTE)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(type,convergenceToleranceType_CONVERGENCE_TOLERANCE_TYPE_RELATIVE)}
		},
	kw_596[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_597[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_598[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_597,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_599[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_600[2] = {
		{"max",8,0,1,1,0,0.,0.,0,N_mdm(type,qoiAggregation_QOI_AGGREGATION_MAX)},
		{"sum",8,0,1,1,0,0.,0.,0,N_mdm(type,qoiAggregation_QOI_AGGREGATION_SUM)}
		},
	kw_601[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_602[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_603[17] = {
		{"allocation_target",8,4,6,0,kw_593},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"convergence_tolerance_target",8,2,10,0,kw_594},
		{"convergence_tolerance_type",8,2,9,0,kw_595},
		{"distribution",8,2,14,0,kw_596},
		{"export_sample_sequence",8,3,5,0,kw_598,0.,0.,0,N_mdm(true,exportSampleSeqFlag)},
		{"final_moments",8,3,13,0,kw_599},
		{"fixed_seed",8,0,2,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"initial_samples",5,0,3,0,0,0.,0.,4,N_mdm(szarray,pilotSamples)},
		{"max_function_evaluations",0x29,0,12,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,16,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"pilot_samples",13,0,3,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"qoi_aggregation",8,2,7,0,kw_600},
		{"rng",8,2,15,0,kw_601},
		{"sample_type",8,2,4,0,kw_602},
		{"seed_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)}
		},
	kw_604[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_605[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_606[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_605,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_607[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_608[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_609[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_610[13] = {
		{"convergence_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,10,0,kw_604},
		{"export_sample_sequence",8,3,5,0,kw_606,0.,0.,0,N_mdm(true,exportSampleSeqFlag)},
		{"final_moments",8,3,9,0,kw_607},
		{"fixed_seed",8,0,2,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"initial_samples",5,0,3,0,0,0.,0.,4,N_mdm(szarray,pilotSamples)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"pilot_samples",13,0,3,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"rng",8,2,11,0,kw_608},
		{"sample_type",8,2,4,0,kw_609},
		{"seed_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)}
		},
	kw_611[3] = {
		{"metric_tracker",8,0,1,1,0,0.,0.,0,N_mdm(lit,convergenceType_metric_tracker)},
		{"num_generations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,2,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_612[2] = {
		{"num_offspring",0x19,0,2,0,0,0.,0.,0,N_mdm(sizet,numOffspring)},
		{"num_parents",0x19,0,1,0,0,0.,0.,0,N_mdm(sizet,numParents)}
		},
	kw_613[5] = {
		{"crossover_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_crossoverType_null_crossover)},
		{"multi_point_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_binary)},
		{"multi_point_parameterized_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_parameterized_binary)},
		{"multi_point_real",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_real)},
		{"shuffle_random",8,2,1,1,kw_612,0.,0.,0,N_mdm(litc,TYPE_DATA_crossoverType_shuffle_random)}
		},
	kw_614[2] = {
		{"domination_count",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_domination_count)},
		{"layer_rank",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_layer_rank)}
		},
	kw_615[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_616[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_617[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_618[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_619[6] = {
		{"bit_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_bit_random)},
		{"mutation_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_mutationType_null_mutation)},
		{"offset_cauchy",8,1,1,1,kw_616,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_cauchy)},
		{"offset_normal",8,1,1,1,kw_617,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_normal)},
		{"offset_uniform",8,1,1,1,kw_618,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_620[1] = {
		{"num_designs",0x29,0,1,0,0,2.,0.,0,N_mdm(sizet,numDesigns)}
		},
	kw_621[3] = {
		{"distance",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_distance)},
		{"max_designs",14,1,1,1,kw_620,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_max_designs)},
		{"radial",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_radial)}
		},
	kw_622[1] = {
		{"orthogonal_distance",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_postProcessorType_distance_postprocessor)}
		},
	kw_623[2] = {
		{"shrinkage_fraction",10,0,1,0,0,0.,0.,0,N_mdm(Real01,shrinkagePercent)},
		{"shrinkage_percentage",2,0,1,0,0,0.,0.,-1,N_mdm(Real01,shrinkagePercent)}
		},
	kw_624[4] = {
		{"below_limit",10,2,1,1,kw_623,0.,0.,0,N_mdm(litp,TYPE_DATA_replacementType_below_limit)},
		{"elitist",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_elitist)},
		{"roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_roulette_wheel)},
		{"unique_roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_unique_roulette_wheel)}
		},
	kw_625[17] = {
		{"convergence_tolerance",10,0,16,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"convergence_type",8,3,4,0,kw_611},
		{"crossover_type",8,5,13,0,kw_613},
		{"fitness_type",8,2,1,0,kw_614},
		{"initialization_type",8,3,12,0,kw_615},
		{"log_file",11,0,10,0,0,0.,0.,0,N_mdm(str,logFile)},
		{"max_function_evaluations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,17,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_type",8,6,14,0,kw_619},
		{"niching_type",8,3,3,0,kw_621},
		{"population_size",0x29,0,9,0,0,0.,0.,0,N_mdm(int,populationSize)},
		{"postprocessor_type",8,1,5,0,kw_622},
		{"print_each_pop",8,0,11,0,0,0.,0.,0,N_mdm(true,printPopFlag)},
		{"replacement_type",8,4,2,0,kw_624},
		{"scaling",8,0,8,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,15,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_626[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_627[1] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,subModelPointer)}
		},
	kw_628[1] = {
		{"seed",9,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_629[7] = {
		{"iterator_scheduling",8,2,5,0,kw_626},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(int,iteratorServers)},
		{"method_name",11,1,1,1,kw_627,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(int,procsPerIterator)},
		{"random_starts",9,1,2,0,kw_628,0.,0.,0,N_mdm(int,concurrentRandomJobs)},
		{"starting_points",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,concurrentParameterSets)}
		},
	kw_630[2] = {
		{"model_pointer",11,0,2,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"partitions",13,0,1,1,0,0.,0.,0,N_mdm(usharray,varPartitions)}
		},
	kw_631[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_632[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_633[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_634[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_635[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_634,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_636[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_637[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_638[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_639[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_638,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_640[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_641[1] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)}
		},
	kw_642[5] = {
		{"increment_max_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_ORDER_ADVANCEMENT)},
		{"increment_max_rank",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_RANK_ADVANCEMENT)},
		{"increment_max_rank_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_MAX_RANK_ORDER_ADVANCEMENT)},
		{"increment_start_order",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_START_ORDER_ADVANCEMENT)},
		{"increment_start_rank",8,0,1,1,0,0.,0.,0,N_mdm(type,c3AdvanceType_START_RANK_ADVANCEMENT)}
		},
	kw_643[1] = {
		{"uniform",8,5,1,1,kw_642,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_644[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_645[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_646[1] = {
		{"l2_penalty",10,0,1,1,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)}
		},
	kw_647[2] = {
		{"ls",8,0,1,1,0,0.,0.,0,N_mdm(type,regressionType_FT_LS)},
		{"rls2",8,1,1,1,kw_646,0.,0.,0,N_mdm(type,regressionType_FT_RLS2)}
		},
	kw_648[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_649[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_650[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_649}
		},
	kw_651[2] = {
		{"compute",8,4,2,0,kw_650},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_652[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_653[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_654[2] = {
		{"active",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_ACTIVE_EXPANSION_STATS)},
		{"combined",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_COMBINED_EXPANSION_STATS)}
		},
	kw_655[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_656[51] = {
		{"adapt_order",8,0,19,0,0,0.,0.,0,N_mdm(true,adaptOrder)},
		{"adapt_rank",8,0,24,0,0,0.,0.,0,N_mdm(true,adaptRank)},
		{"allocation_control",8,1,6,0,kw_631},
		{"arithmetic_tolerance",10,0,9,0,0,0.,0.,0,N_mdm(Real,statsRoundingTol)},
		{"collocation_points_sequence",13,0,16,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,0,17,0,0,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,39,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,7,0,kw_632},
		{"distribution",8,2,37,0,kw_633},
		{"export_approx_points_file",11,3,41,0,kw_635,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,41,0,kw_635,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,32,0,kw_636},
		{"fixed_seed",8,0,43,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,39,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,36,0,kw_637,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,40,0,kw_639,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"kick_order",0x19,0,20,0,0,0.,0.,0,N_mdm(ushint,kickOrder)},
		{"kick_rank",0x19,0,25,0,0,0.,0.,0,N_mdm(sizet,kickRank)},
		{"max_cross_iterations",0x29,0,12,0,0,0.,0.,0,N_mdm(int,maxCrossIterations)},
		{"max_cv_order_candidates",0x29,0,22,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"max_cv_rank_candidates",0x29,0,27,0,0,0.,0.,0,N_mdm(sizet,maxCVRankCandidates)},
		{"max_order",0x29,0,21,0,0,0.,0.,0,N_mdm(ushint,maxOrder)},
		{"max_rank",0x29,0,26,0,0,0.,0.,0,N_mdm(sizet,maxRank)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"max_solver_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"metric_scale",8,2,4,0,kw_640},
		{"model_pointer",11,0,44,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"order_sequence",5,1,18,0,kw_641,0.,0.,18,N_mdm(usharray,startOrderSeq)},
		{"p_refinement",8,1,1,0,kw_643,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"pilot_samples",5,0,16,0,0,0.,0.,-26,N_mdm(szarray,collocationPointsSeq)},
		{"probability_levels",14,1,34,0,kw_644,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,31,0,kw_645},
		{"rank_sequence",5,0,23,0,0,0.,0.,14,N_mdm(szarray,startRankSeq)},
		{"regression_type",8,2,10,0,kw_647},
		{"reliability_levels",14,1,35,0,kw_648,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,33,0,kw_651,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"response_scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"rng",8,2,30,0,kw_652},
		{"rounding_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,solverRoundingTol)},
		{"sample_refinement",0,4,31,0,kw_645,0.,0.,-8},
		{"sample_type",8,2,29,0,kw_653},
		{"samples",1,0,28,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,28,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed_sequence",13,0,42,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"solver_tolerance",10,0,13,0,0,0.,0.,0,N_mdm(Real,solverTol)},
		{"start_order_sequence",13,1,18,0,kw_641,0.,0.,0,N_mdm(usharray,startOrderSeq)},
		{"start_rank_sequence",13,0,23,0,0,0.,0.,0,N_mdm(szarray,startRankSeq)},
		{"statistics_mode",8,2,5,0,kw_654},
		{"tensor_grid",8,0,15,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"variance_based_decomp",8,2,38,0,kw_655,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_657[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_658[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_659[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_660[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_661[3] = {
		{"adapted",8,2,1,1,kw_660,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_662[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_663[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_664[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_665[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_666[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_667[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_668[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_662,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_662,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_663,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_664,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_665,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_665,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_664,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_666,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_667,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_667,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_669[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_670[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_671[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_670,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_672[6] = {
		{"basis_type",8,3,2,0,kw_661},
		{"collocation_ratio",10,21,3,1,kw_668,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_669,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_671,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_671,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_673[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_674[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_673,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_675[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_676[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_677[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_678[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_677,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_679[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_680[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_679,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_681[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_680,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_680,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_682[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_683[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_684[2] = {
		{"dimension_adaptive",8,3,1,1,kw_683},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_685[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_686[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_687[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_688[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_689[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_690[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_689}
		},
	kw_691[2] = {
		{"compute",8,4,2,0,kw_690},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_692[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_693[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_694[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_695[2] = {
		{"active",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_ACTIVE_EXPANSION_STATS)},
		{"combined",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_COMBINED_EXPANSION_STATS)}
		},
	kw_696[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_697[38] = {
		{"allocation_control",8,1,6,0,kw_657},
		{"askey",8,0,9,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,23,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,7,0,kw_658},
		{"distribution",8,2,21,0,kw_659},
		{"expansion_order_sequence",13,6,8,1,kw_672,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_approx_points_file",11,3,25,0,kw_674,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,11,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,25,0,kw_674,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,16,0,kw_675},
		{"fixed_seed",8,0,27,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,23,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,20,0,kw_676,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,24,0,kw_678,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"least_interpolation",0,6,8,1,kw_681,0.,0.,6,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"metric_scale",8,2,4,0,kw_682},
		{"model_pointer",11,0,28,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,10,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,8,1,kw_681,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,8,1,kw_681,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_684,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"probability_levels",14,1,18,0,kw_685,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,15,0,kw_686},
		{"quadrature_order_sequence",13,3,8,1,kw_687,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"reliability_levels",14,1,19,0,kw_688,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,17,0,kw_691,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,14,0,kw_692},
		{"sample_refinement",0,4,15,0,kw_686,0.,0.,-5},
		{"sample_type",8,2,13,0,kw_693},
		{"samples",1,0,12,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,12,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed_sequence",13,0,26,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"sparse_grid_level_sequence",13,5,8,1,kw_694,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"statistics_mode",8,2,5,0,kw_695},
		{"variance_based_decomp",8,2,22,0,kw_696,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,9,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_698[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_GREEDY_REFINEMENT)}
		},
	kw_699[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_700[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_701[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_702[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_701,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_703[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_704[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_705[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_706[3] = {
		{"dimension_adaptive",8,2,1,1,kw_705},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_707[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_708[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_707,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_709[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_710[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_711[2] = {
		{"dimension_adaptive",8,2,1,1,kw_710},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_712[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_713[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_714[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_715[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_716[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_717[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_716}
		},
	kw_718[2] = {
		{"compute",8,4,2,0,kw_717},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_719[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_720[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_721[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_722[2] = {
		{"active",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_ACTIVE_EXPANSION_STATS)},
		{"combined",8,0,1,1,0,0.,0.,0,N_mdm(type,statsMetricMode_COMBINED_EXPANSION_STATS)}
		},
	kw_723[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_724[35] = {
		{"allocation_control",8,1,6,0,kw_698},
		{"askey",8,0,9,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,7,0,kw_699},
		{"distribution",8,2,20,0,kw_700},
		{"export_approx_points_file",11,3,24,0,kw_702,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,24,0,kw_702,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,15,0,kw_703},
		{"fixed_seed",8,0,26,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,19,0,kw_704,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"h_refinement",8,3,1,0,kw_706,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"import_approx_points_file",11,4,23,0,kw_708,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"metric_scale",8,2,4,0,kw_709},
		{"model_pointer",11,0,27,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"p_refinement",8,2,1,0,kw_711,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,9,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"probability_levels",14,1,17,0,kw_712,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,14,0,kw_713},
		{"quadrature_order_sequence",13,3,8,1,kw_714,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"reliability_levels",14,1,18,0,kw_715,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,16,0,kw_718,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,13,0,kw_719},
		{"sample_refinement",0,4,14,0,kw_713,0.,0.,-5},
		{"sample_type",8,2,12,0,kw_720},
		{"samples",1,0,11,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,11,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed_sequence",13,0,25,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"sparse_grid_level_sequence",13,7,8,1,kw_721,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"statistics_mode",8,2,5,0,kw_722},
		{"use_derivatives",8,0,10,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"variance_based_decomp",8,2,21,0,kw_723,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,9,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_725[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_726[2] = {
		{"estimator_variance",8,1,1,1,kw_725,0.,0.,0,N_mdm(type,multilevAllocControl_ESTIMATOR_VARIANCE)},
		{"rank_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_RANK_SAMPLING)}
		},
	kw_727[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_728[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_729[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_730[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_729,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_731[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_732[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_733[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_734[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_733,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_735[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_736[1] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)}
		},
	kw_737[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_738[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_739[1] = {
		{"l2_penalty",10,0,1,1,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)}
		},
	kw_740[2] = {
		{"ls",8,0,1,1,0,0.,0.,0,N_mdm(type,regressionType_FT_LS)},
		{"rls2",8,1,1,1,kw_739,0.,0.,0,N_mdm(type,regressionType_FT_RLS2)}
		},
	kw_741[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_742[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_743[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_742}
		},
	kw_744[2] = {
		{"compute",8,4,2,0,kw_743},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_745[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_746[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_747[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_748[49] = {
		{"adapt_order",8,0,17,0,0,0.,0.,0,N_mdm(true,adaptOrder)},
		{"adapt_rank",8,0,22,0,0,0.,0.,0,N_mdm(true,adaptRank)},
		{"allocation_control",8,2,2,0,kw_726},
		{"arithmetic_tolerance",10,0,7,0,0,0.,0.,0,N_mdm(Real,statsRoundingTol)},
		{"collocation_points_sequence",13,0,14,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,0,15,0,0,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,37,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,5,0,kw_727},
		{"distribution",8,2,35,0,kw_728},
		{"export_approx_points_file",11,3,39,0,kw_730,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,39,0,kw_730,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,30,0,kw_731},
		{"fixed_seed",8,0,41,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,37,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,34,0,kw_732,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,38,0,kw_734,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"kick_order",0x19,0,18,0,0,0.,0.,0,N_mdm(ushint,kickOrder)},
		{"kick_rank",0x19,0,23,0,0,0.,0.,0,N_mdm(sizet,kickRank)},
		{"max_cross_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(int,maxCrossIterations)},
		{"max_cv_order_candidates",0x29,0,20,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"max_cv_rank_candidates",0x29,0,25,0,0,0.,0.,0,N_mdm(sizet,maxCVRankCandidates)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_order",0x29,0,19,0,0,0.,0.,0,N_mdm(ushint,maxOrder)},
		{"max_rank",0x29,0,24,0,0,0.,0.,0,N_mdm(sizet,maxRank)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"metric_scale",8,2,4,0,kw_735},
		{"model_pointer",11,0,42,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"order_sequence",5,1,16,0,kw_736,0.,0.,17,N_mdm(usharray,startOrderSeq)},
		{"pilot_samples",5,0,14,0,0,0.,0.,-25,N_mdm(szarray,collocationPointsSeq)},
		{"probability_levels",14,1,32,0,kw_737,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,29,0,kw_738},
		{"rank_sequence",5,0,21,0,0,0.,0.,14,N_mdm(szarray,startRankSeq)},
		{"regression_type",8,2,8,0,kw_740},
		{"reliability_levels",14,1,33,0,kw_741,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,31,0,kw_744,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"response_scaling",8,0,12,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"rng",8,2,28,0,kw_745},
		{"rounding_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,solverRoundingTol)},
		{"sample_refinement",0,4,29,0,kw_738,0.,0.,-8},
		{"sample_type",8,2,27,0,kw_746},
		{"samples",1,0,26,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,26,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed_sequence",13,0,40,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"solver_tolerance",10,0,11,0,0,0.,0.,0,N_mdm(Real,solverTol)},
		{"start_order_sequence",13,1,16,0,kw_736,0.,0.,0,N_mdm(usharray,startOrderSeq)},
		{"start_rank_sequence",13,0,21,0,0,0.,0.,0,N_mdm(szarray,startRankSeq)},
		{"tensor_grid",8,0,13,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"variance_based_decomp",8,2,36,0,kw_747,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_749[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_750[2] = {
		{"estimator_variance",8,1,1,1,kw_749,0.,0.,0,N_mdm(type,multilevAllocControl_ESTIMATOR_VARIANCE)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevAllocControl_RIP_SAMPLING)}
		},
	kw_751[3] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"paired",0,0,1,1,0,0.,0.,-1,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_752[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_753[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_754[3] = {
		{"adapted",8,2,1,1,kw_753,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_755[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_756[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_757[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_758[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_759[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_760[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_761[21] = {
		{"basis_pursuit",8,0,2,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,2,0,kw_755,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,2,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,2,0,kw_755,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"collocation_points_sequence",13,0,1,0,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"cross_validation",8,2,3,0,kw_756,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,2,0,kw_757,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,2,0,kw_758,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,2,0,kw_758,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,2,0,kw_757,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,2,0,kw_759,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,2,0,kw_760,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,2,0,kw_760,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"pilot_samples",5,0,1,0,0,0.,0.,-10,N_mdm(szarray,collocationPointsSeq)},
		{"ratio_order",10,0,4,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,8,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,8,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,7,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_762[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_763[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_764[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_763,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_765[6] = {
		{"basis_type",8,3,2,0,kw_754},
		{"collocation_ratio",10,21,3,1,kw_761,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,2,3,1,kw_762,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_764,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_764,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_766[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_767[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_766,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_768[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_769[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_770[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_771[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_770,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_772[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_773[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_772,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_774[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_773,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_773,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_775[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_776[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_777[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_778[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_779[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_780[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_779}
		},
	kw_781[2] = {
		{"compute",8,4,2,0,kw_780},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_782[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_783[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_784[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_785[34] = {
		{"allocation_control",8,2,2,0,kw_750},
		{"askey",8,0,7,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,21,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,3,5,0,kw_751},
		{"distribution",8,2,19,0,kw_752},
		{"expansion_order_sequence",13,6,6,1,kw_765,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_approx_points_file",11,3,23,0,kw_767,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,9,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,23,0,kw_767,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,14,0,kw_768},
		{"fixed_seed",8,0,25,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,21,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,18,0,kw_769,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,22,0,kw_771,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"least_interpolation",0,6,6,1,kw_774,0.,0.,6,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"metric_scale",8,2,4,0,kw_775},
		{"model_pointer",11,0,26,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,8,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,6,1,kw_774,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,6,1,kw_774,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"probability_levels",14,1,16,0,kw_776,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,13,0,kw_777},
		{"reliability_levels",14,1,17,0,kw_778,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,15,0,kw_781,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,12,0,kw_782},
		{"sample_refinement",0,4,13,0,kw_777,0.,0.,-4},
		{"sample_type",8,2,11,0,kw_783},
		{"samples",1,0,10,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,10,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed_sequence",13,0,24,0,0,0.,0.,0,N_mdm(szarray,randomSeedSeq)},
		{"variance_based_decomp",8,2,20,0,kw_784,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,7,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_786[9] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"min_boxsize_limit",10,0,2,0,0,0.,0.,0,N_mdm(Real,minBoxSize)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,7,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"solution_accuracy",2,0,1,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,1,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"volume_boxsize_limit",10,0,3,0,0,0.,0.,0,N_mdm(Real,volBoxSize)}
		},
	kw_787[15] = {
		{"absolute_conv_tol",10,0,2,0,0,0.,0.,0,N_mdm(Real,absConvTol)},
		{"convergence_tolerance",10,0,10,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"covariance",9,0,8,0,0,0.,0.,0,N_mdm(int,covarianceType)},
		{"false_conv_tol",10,0,6,0,0,0.,0.,0,N_mdm(Real,falseConvTol)},
		{"function_precision",10,0,1,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"initial_trust_radius",10,0,7,0,0,0.,0.,0,N_mdm(Real,initTRRadius)},
		{"max_function_evaluations",0x29,0,13,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"regression_diagnostics",8,0,9,0,0,0.,0.,0,N_mdm(true,regressDiag)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"singular_conv_tol",10,0,4,0,0,0.,0.,0,N_mdm(Real,singConvTol)},
		{"singular_radius",10,0,5,0,0,0.,0.,0,N_mdm(Real,singRadius)},
		{"speculative",8,0,12,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"x_conv_tol",10,0,3,0,0,0.,0.,0,N_mdm(Real,xConvTol)}
		},
	kw_788[5] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_789[10] = {
		{"constraint_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"function_precision",10,0,2,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"linesearch_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,lineSearchTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,7,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"verify_level",9,0,1,0,0,0.,0.,0,N_mdm(int,verifyLevel)}
		},
	kw_790[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_791[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_792[2] = {
		{"global",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_global)},
		{"local",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_local)}
		},
	kw_793[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_794[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_795[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_794}
		},
	kw_796[2] = {
		{"compute",8,3,2,0,kw_795},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_797[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_798[11] = {
		{"build_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,8,0,kw_790},
		{"gen_reliability_levels",14,1,7,0,kw_791,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lipschitz",8,2,3,0,kw_792},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,6,0,kw_793,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,5,0,kw_796,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,9,0,kw_797},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,4,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_799[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_800[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_801[3] = {
		{"adapted",8,2,1,1,kw_800,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_802[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_803[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_804[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_805[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_806[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_807[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_808[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_802,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_802,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_803,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_804,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_805,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_805,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_804,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_806,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_807,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_807,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_809[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_810[2] = {
		{"max_cv_order_candidates",0x29,0,2,0,0,0.,0.,0,N_mdm(ushint,maxCVOrderCandidates)},
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_811[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_812[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_813[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_814[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_815[19] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_809,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_809,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,2,2,0,kw_810,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_811,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_812,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_812,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_811,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_813,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxSolverIterations)},
		{"omp",0,1,1,0,kw_814,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_814,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"response_scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,respScalingFlag)},
		{"reuse_points",8,0,7,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,7,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_816[2] = {
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_817[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_818[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_817,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_819[7] = {
		{"basis_type",8,3,2,0,kw_801},
		{"collocation_points",9,19,3,1,kw_808,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,19,3,1,kw_815,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,2,3,1,kw_816,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_818,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_818,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_820[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_821[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_820,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_822[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_823[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_824[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_825[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_824,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_826[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_827[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_826,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_828[6] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_827,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_827,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_829[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_830[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_831[2] = {
		{"dimension_adaptive",8,3,1,1,kw_830},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_832[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_833[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_834[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_835[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_836[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_837[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_836}
		},
	kw_838[2] = {
		{"compute",8,4,2,0,kw_837},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_839[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_840[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_841[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_842[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_843[37] = {
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"cubature_integrand",9,0,5,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,20,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,18,0,kw_799},
		{"expansion_order",9,7,5,1,kw_819,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_approx_points_file",11,3,22,0,kw_821,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,22,0,kw_821,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,13,0,kw_822},
		{"fixed_seed",8,0,24,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,20,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,17,0,kw_823,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,21,0,kw_825,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"import_expansion_file",11,0,5,1,0,0.,0.,0,N_mdm(str,importExpansionFile)},
		{"least_interpolation",0,6,5,1,kw_828,0.,0.,6,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"metric_scale",8,2,4,0,kw_829},
		{"model_pointer",11,0,25,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_828,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_828,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_831,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"probability_levels",14,1,15,0,kw_832,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,12,0,kw_833},
		{"quadrature_order",9,3,5,1,kw_834,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"reliability_levels",14,1,16,0,kw_835,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,14,0,kw_838,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_839},
		{"sample_refinement",0,4,12,0,kw_833,0.,0.,-5},
		{"sample_type",8,2,10,0,kw_840},
		{"samples",1,0,9,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,9,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,23,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"sparse_grid_level",9,5,5,1,kw_841,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"variance_based_decomp",8,2,19,0,kw_842,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_844[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_845[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_846[2] = {
		{"global",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_global)},
		{"local",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_local)}
		},
	kw_847[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_848[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_849[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_848}
		},
	kw_850[2] = {
		{"compute",8,3,2,0,kw_849},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_851[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_852[11] = {
		{"build_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,8,0,kw_844},
		{"gen_reliability_levels",14,1,7,0,kw_845,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lipschitz",8,2,3,0,kw_846},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,6,0,kw_847,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,5,0,kw_850,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,9,0,kw_851},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,4,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_853[2] = {
		{"candidate_designs",0x19,0,1,0,0,0.,0.,0,N_mdm(sizet,numCandidateDesigns)},
		{"leja_oversample_ratio",10,0,1,0,0,0.,0.,0,N_mdm(Real,collocationRatio)}
		},
	kw_854[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_855[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_856[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_857[1] = {
		{"percent_variance_explained",10,0,1,0,0,0.,0.,0,N_mdm(Real,percentVarianceExplained)}
		},
	kw_858[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_859[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_860[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_861[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_860}
		},
	kw_862[2] = {
		{"compute",8,4,2,0,kw_861},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_863[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_864[4] = {
		{"incremental_lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"incremental_random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)},
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_865[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_866[5] = {
		{"confidence_level",10,0,2,0,0,0.,0.,0,N_mdm(Real,wilksConfidenceLevel)},
		{"one_sided_lower",8,0,3,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_ONE_SIDED_LOWER)},
		{"one_sided_upper",8,0,4,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_ONE_SIDED_UPPER)},
		{"order",9,0,1,0,0,0.,0.,0,N_mdm(ushint,wilksOrder)},
		{"two_sided",8,0,5,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_TWO_SIDED)}
		},
	kw_867[19] = {
		{"backfill",8,0,8,0,0,0.,0.,0,N_mdm(true,backfillFlag)},
		{"d_optimal",8,2,6,0,kw_853,0.,0.,0,N_mdm(true,dOptimal)},
		{"distribution",8,2,16,0,kw_854},
		{"final_moments",8,3,11,0,kw_855},
		{"fixed_seed",8,0,3,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"gen_reliability_levels",14,1,15,0,kw_856,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"initial_samples",1,0,1,0,0,0.,0.,9,N_mdm(int,numSamples)},
		{"model_pointer",11,0,18,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"principal_components",8,1,9,0,kw_857,0.,0.,0,N_mdm(true,pcaFlag)},
		{"probability_levels",14,1,13,0,kw_858,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,5,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"reliability_levels",14,1,14,0,kw_859,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,12,0,kw_862,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,17,0,kw_863},
		{"sample_type",8,4,4,0,kw_864},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"variance_based_decomp",8,1,7,0,kw_865,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wilks",8,5,10,0,kw_866,0.,0.,0,N_mdm(true,wilksFlag)}
		},
	kw_868[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_869[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_870[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_869,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_871[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_872[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_873[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_874[3] = {
		{"dimension_adaptive",8,2,1,1,kw_873},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_875[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_876[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_875,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_877[2] = {
		{"absolute",8,0,1,1,0,0.,0.,0,N_mdm(false,relativeConvMetric)},
		{"relative",8,0,1,1,0,0.,0.,0,N_mdm(true,relativeConvMetric)}
		},
	kw_878[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_879[2] = {
		{"dimension_adaptive",8,2,1,1,kw_878},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_880[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_881[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_882[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_883[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_884[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_885[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_884}
		},
	kw_886[2] = {
		{"compute",8,4,2,0,kw_885},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_887[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_888[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_889[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_890[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_891[32] = {
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,19,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,17,0,kw_868},
		{"export_approx_points_file",11,3,21,0,kw_870,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,21,0,kw_870,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,12,0,kw_871},
		{"fixed_seed",8,0,23,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,19,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,16,0,kw_872,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"h_refinement",8,3,1,0,kw_874,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"import_approx_points_file",11,4,20,0,kw_876,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxRefineIterations)},
		{"metric_scale",8,2,4,0,kw_877},
		{"model_pointer",11,0,24,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"p_refinement",8,2,1,0,kw_879,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"probability_levels",14,1,14,0,kw_880,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,11,0,kw_881},
		{"quadrature_order",9,3,5,1,kw_882,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"reliability_levels",14,1,15,0,kw_883,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,13,0,kw_886,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,10,0,kw_887},
		{"sample_refinement",0,4,11,0,kw_881,0.,0.,-5},
		{"sample_type",8,2,9,0,kw_888},
		{"samples",1,0,8,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,8,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,22,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"sparse_grid_level",9,7,5,1,kw_889,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"variance_based_decomp",8,2,18,0,kw_890,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_892[5] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"misc_options",15,0,1,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_893[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_894[5] = {
		{"max_function_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"trust_region",8,6,1,0,kw_893,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final}
		},
	kw_895[10] = {
		{"constraint_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"function_precision",10,0,2,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"linesearch_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,lineSearchTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,7,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"verify_level",9,0,1,0,0,0.,0.,0,N_mdm(int,verifyLevel)}
		},
	kw_896[8] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,6,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_step",10,0,1,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,7,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,5,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_897[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_898[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_899[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_897},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_898},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_900[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_901[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_902[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_900},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_901},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_903[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_904[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_905[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_903},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_904},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_906[6] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_scheme_size",9,0,1,0,0,0.,0.,0,N_mdm(int,searchSchemeSize)}
		},
	kw_907[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_908[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_909[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_907},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_908},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_910[5] = {
		{"debug",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_DEBUG_OUTPUT)},
		{"normal",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_NORMAL_OUTPUT)},
		{"quiet",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_QUIET_OUTPUT)},
		{"silent",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_SILENT_OUTPUT)},
		{"verbose",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_VERBOSE_OUTPUT)}
		},
	kw_911[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_912[2] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,subModelPointer)},
		{"opt_model_pointer",3,0,1,0,0,0.,0.,-1,N_mdm(str,subModelPointer)}
		},
	kw_913[1] = {
		{"seed",9,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_914[10] = {
		{"iterator_scheduling",8,2,5,0,kw_911},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(int,iteratorServers)},
		{"method_name",11,2,1,1,kw_912,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"multi_objective_weight_sets",6,0,3,0,0,0.,0.,5,N_mdm(RealDL,concurrentParameterSets)},
		{"opt_method_name",3,2,1,1,kw_912,0.,0.,-3,N_mdm(str,subMethodName)},
		{"opt_method_pointer",3,0,1,1,0,0.,0.,-3,N_mdm(str,subMethodPointer)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(int,procsPerIterator)},
		{"random_weight_sets",9,1,2,0,kw_913,0.,0.,0,N_mdm(int,concurrentRandomJobs)},
		{"weight_sets",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,concurrentParameterSets)}
		},
	kw_915[4] = {
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"partitions",13,0,1,0,0,0.,0.,0,N_mdm(usharray,varPartitions)},
		{"samples",9,0,2,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_916[7] = {
		{"converge_order",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_CONVERGE_ORDER)},
		{"converge_qoi",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_CONVERGE_QOI)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"estimate_order",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_ESTIMATE_ORDER)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"refinement_rate",10,0,2,0,0,0.,0.,0,N_mdm(Real,refinementRate)}
		},
	kw_917[7] = {
		{"constraint_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"gradient_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"options_file",11,0,5,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_918[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_919[6] = {
		{"max_function_evaluations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"trust_region",8,6,2,0,kw_918,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final}
		},
	kw_920[2] = {
		{"num_generations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,1,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_921[2] = {
		{"num_generations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,1,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_922[2] = {
		{"average_fitness_tracker",8,2,1,1,kw_920,0.,0.,0,N_mdm(lit,convergenceType_average_fitness_tracker)},
		{"best_fitness_tracker",8,2,1,1,kw_921,0.,0.,0,N_mdm(lit,convergenceType_best_fitness_tracker)}
		},
	kw_923[2] = {
		{"num_offspring",0x19,0,2,0,0,0.,0.,0,N_mdm(sizet,numOffspring)},
		{"num_parents",0x19,0,1,0,0,0.,0.,0,N_mdm(sizet,numParents)}
		},
	kw_924[5] = {
		{"crossover_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_crossoverType_null_crossover)},
		{"multi_point_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_binary)},
		{"multi_point_parameterized_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_parameterized_binary)},
		{"multi_point_real",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_real)},
		{"shuffle_random",8,2,1,1,kw_923,0.,0.,0,N_mdm(litc,TYPE_DATA_crossoverType_shuffle_random)}
		},
	kw_925[2] = {
		{"constraint_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Realp,constraintTolerance)},
		{"merit_function",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_merit_function)}
		},
	kw_926[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_927[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_928[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_929[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_930[6] = {
		{"bit_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_bit_random)},
		{"mutation_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_mutationType_null_mutation)},
		{"offset_cauchy",8,1,1,1,kw_927,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_cauchy)},
		{"offset_normal",8,1,1,1,kw_928,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_normal)},
		{"offset_uniform",8,1,1,1,kw_929,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_931[4] = {
		{"elitist",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_elitist)},
		{"favor_feasible",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_favor_feasible)},
		{"roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_roulette_wheel)},
		{"unique_roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_unique_roulette_wheel)}
		},
	kw_932[15] = {
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"convergence_type",8,2,3,0,kw_922},
		{"crossover_type",8,5,11,0,kw_924},
		{"fitness_type",8,2,1,0,kw_925},
		{"initialization_type",8,3,10,0,kw_926},
		{"log_file",11,0,8,0,0,0.,0.,0,N_mdm(str,logFile)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(sizet,maxFunctionEvals)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_type",8,6,12,0,kw_930},
		{"population_size",0x29,0,7,0,0,0.,0.,0,N_mdm(int,populationSize)},
		{"print_each_pop",8,0,9,0,0,0.,0.,0,N_mdm(true,printPopFlag)},
		{"replacement_type",8,4,2,0,kw_931},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,13,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_933[8] = {
		{"approx_method_name",3,0,1,1,0,0.,0.,4,N_mdm(str,subMethodName)},
		{"approx_method_pointer",3,0,1,1,0,0.,0.,4,N_mdm(str,subMethodPointer)},
		{"approx_model_pointer",3,0,2,2,0,0.,0.,4,N_mdm(str,modelPointer)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"method_name",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"model_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"replace_points",8,0,3,0,0,0.,0.,0,N_mdm(true,surrBasedGlobalReplacePts)}
		},
	kw_934[2] = {
		{"filter",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalAcceptLogic_FILTER)},
		{"tr_ratio",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalAcceptLogic_TR_RATIO)}
		},
	kw_935[7] = {
		{"augmented_lagrangian_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_AUGMENTED_LAGRANGIAN_OBJECTIVE)},
		{"lagrangian_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_LAGRANGIAN_OBJECTIVE)},
		{"linearized_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_LINEARIZED_CONSTRAINTS)},
		{"no_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_NO_CONSTRAINTS)},
		{"original_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_ORIGINAL_CONSTRAINTS)},
		{"original_primary",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_ORIGINAL_PRIMARY)},
		{"single_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_SINGLE_OBJECTIVE)}
		},
	kw_936[1] = {
		{"homotopy",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalConstrRelax_HOMOTOPY)}
		},
	kw_937[4] = {
		{"adaptive_penalty_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_ADAPTIVE_PENALTY_MERIT)},
		{"augmented_lagrangian_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_AUGMENTED_LAGRANGIAN_MERIT)},
		{"lagrangian_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_LAGRANGIAN_MERIT)},
		{"penalty_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_PENALTY_MERIT)}
		},
	kw_938[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_939[16] = {
		{"acceptance_logic",8,2,7,0,kw_934},
		{"approx_method_name",3,0,1,1,0,0.,0.,9,N_mdm(str,subMethodName)},
		{"approx_method_pointer",3,0,1,1,0,0.,0.,9,N_mdm(str,subMethodPointer)},
		{"approx_model_pointer",3,0,2,2,0,0.,0.,9,N_mdm(str,modelPointer)},
		{"approx_subproblem",8,7,5,0,kw_935},
		{"constraint_relax",8,1,8,0,kw_936},
		{"constraint_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,11,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(sizet,maxIterations)},
		{"merit_function",8,4,6,0,kw_937},
		{"method_name",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"model_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"soft_convergence_limit",9,0,3,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)},
		{"trust_region",8,6,9,0,kw_938,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final},
		{"truth_surrogate_bypass",8,0,4,0,0,0.,0.,0,N_mdm(true,surrBasedLocalLayerBypass)}
		},
	kw_940[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_941[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_942[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_941,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_943[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_944[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_945[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_946[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_945,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_947[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_948[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_949[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_950[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_951[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_950}
		},
	kw_952[2] = {
		{"compute",8,4,2,0,kw_951},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_953[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_954[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_955[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_956[21] = {
		{"diagonal_covariance",8,0,12,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,10,0,kw_940},
		{"export_approx_points_file",11,3,14,0,kw_942,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,14,0,kw_942,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,5,0,kw_943},
		{"fixed_seed",8,0,16,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,12,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,9,0,kw_944,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,13,0,kw_946,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"model_pointer",11,0,17,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,7,0,kw_947,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,4,0,kw_948},
		{"reliability_levels",14,1,8,0,kw_949,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,6,0,kw_952,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,3,0,kw_953},
		{"sample_refinement",0,4,4,0,kw_948,0.,0.,-4},
		{"sample_type",8,2,2,0,kw_954},
		{"samples",1,0,1,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,1,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,15,0,0,0.,0.,0,N_mdm(int,randomSeed)},
		{"variance_based_decomp",8,2,11,0,kw_955,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_957[4] = {
		{"final_point",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,finalPoint)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"num_steps",9,0,2,2,0,0.,0.,0,N_mdm(int,numSteps)},
		{"step_vector",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,stepVector)}
		},
	kw_958[104] = {
		{"acv_sampling",0,23,4,1,kw_40,0.,0.,2,N_mdm(utype,methodName_APPROXIMATE_CONTROL_VARIATE)},
		{"adaptive_sampling",8,19,4,1,kw_54,0.,0.,0,N_mdm(utype,methodName_ADAPTIVE_SAMPLING)},
		{"approximate_control_variate",8,23,4,1,kw_40,0.,0.,0,N_mdm(utype,methodName_APPROXIMATE_CONTROL_VARIATE)},
		{"asynch_pattern_search",8,13,4,1,kw_57,0.,0.,0,N_mdm(utype,methodName_ASYNCH_PATTERN_SEARCH)},
		{"bayes_calibration",8,18,4,1,kw_390,0.,0.,0,N_mdm(utype,methodName_BAYES_CALIBRATION)},
		{"branch_and_bound",8,3,4,1,kw_392,0.,0.,0,N_mdm(utype,methodName_BRANCH_AND_BOUND)},
		{"centered_parameter_study",8,4,4,1,kw_393,0.,0.,0,N_mdm(utype,methodName_CENTERED_PARAMETER_STUDY)},
		{"coliny_apps",0,13,4,1,kw_57,0.,0.,-4,N_mdm(utype,methodName_ASYNCH_PATTERN_SEARCH)},
		{"coliny_beta",8,11,4,1,kw_394,0.,0.,0,N_mdm(utype,methodName_COLINY_BETA)},
		{"coliny_cobyla",8,12,4,1,kw_395,0.,0.,0,N_mdm(utype,methodName_COLINY_COBYLA)},
		{"coliny_direct",8,16,4,1,kw_397,0.,0.,0,N_mdm(utype,methodName_COLINY_DIRECT)},
		{"coliny_ea",8,19,4,1,kw_406,0.,0.,0,N_mdm(utype,methodName_COLINY_EA)},
		{"coliny_pattern_search",8,22,4,1,kw_410,0.,0.,0,N_mdm(utype,methodName_COLINY_PATTERN_SEARCH)},
		{"coliny_solis_wets",8,18,4,1,kw_411,0.,0.,0,N_mdm(utype,methodName_COLINY_SOLIS_WETS)},
		{"conmin_frcg",8,7,4,1,kw_412,0.,0.,0,N_mdm(utype,methodName_CONMIN_FRCG)},
		{"conmin_mfd",8,7,4,1,kw_413,0.,0.,0,N_mdm(utype,methodName_CONMIN_MFD)},
		{"dace",8,15,4,1,kw_415,0.,0.,0,N_mdm(utype,methodName_DACE)},
		{"demo_tpl",8,7,4,1,kw_416,0.,0.,0,N_mdm(utype,methodName_DEMO_TPL)},
		{"dl_solver",11,3,4,1,kw_417,0.,0.,0,N_mdm(utype_lit,TYPE_DATA_methodName_DL_SOLVER)},
		{"dot_bfgs",8,7,4,1,kw_418,0.,0.,0,N_mdm(utype,methodName_DOT_BFGS)},
		{"dot_frcg",8,7,4,1,kw_419,0.,0.,0,N_mdm(utype,methodName_DOT_FRCG)},
		{"dot_mmfd",8,7,4,1,kw_420,0.,0.,0,N_mdm(utype,methodName_DOT_MMFD)},
		{"dot_slp",8,7,4,1,kw_421,0.,0.,0,N_mdm(utype,methodName_DOT_SLP)},
		{"dot_sqp",8,7,4,1,kw_422,0.,0.,0,N_mdm(utype,methodName_DOT_SQP)},
		{"efficient_global",8,14,4,1,kw_436,0.,0.,0,N_mdm(utype,methodName_EFFICIENT_GLOBAL)},
		{"final_solutions",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,numFinalSolutions)},
		{"fsu_cvt",8,10,4,1,kw_439,0.,0.,0,N_mdm(utype,methodName_FSU_CVT)},
		{"fsu_quasi_mc",8,12,4,1,kw_441},
		{"function_train",8,47,4,1,kw_464,0.,0.,0,N_mdm(utype,methodName_C3_FUNCTION_TRAIN)},
		{"gaussian_process_adaptive_importance_sampling",0,15,4,1,kw_476,0.,0.,6,N_mdm(utype,methodName_GPAIS)},
		{"genie_direct",8,4,4,1,kw_477,0.,0.,0,N_mdm(utype,methodName_GENIE_DIRECT)},
		{"genie_opt_darts",8,4,4,1,kw_478,0.,0.,0,N_mdm(utype,methodName_GENIE_OPT_DARTS)},
		{"global_evidence",8,12,4,1,kw_498,0.,0.,0,N_mdm(utype,methodName_GLOBAL_EVIDENCE)},
		{"global_interval_est",8,11,4,1,kw_519,0.,0.,0,N_mdm(utype,methodName_GLOBAL_INTERVAL_EST)},
		{"global_reliability",8,22,4,1,kw_537,0.,0.,0,N_mdm(utype,methodName_GLOBAL_RELIABILITY)},
		{"gpais",8,15,4,1,kw_476,0.,0.,0,N_mdm(utype,methodName_GPAIS)},
		{"hybrid",8,5,4,1,kw_548,0.,0.,0,N_mdm(utype,methodName_HYBRID)},
		{"id_method",11,0,1,0,0,0.,0.,0,N_mdm(str,idMethod)},
		{"importance_sampling",8,15,4,1,kw_556,0.,0.,0,N_mdm(utype,methodName_IMPORTANCE_SAMPLING)},
		{"list_parameter_study",8,3,4,1,kw_559,0.,0.,0,N_mdm(utype,methodName_LIST_PARAMETER_STUDY)},
		{"local_evidence",8,7,4,1,kw_566,0.,0.,0,N_mdm(utype,methodName_LOCAL_EVIDENCE)},
		{"local_interval_est",8,4,4,1,kw_567,0.,0.,0,N_mdm(utype,methodName_LOCAL_INTERVAL_EST)},
		{"local_reliability",8,10,4,1,kw_579,0.,0.,0,N_mdm(utype,methodName_LOCAL_RELIABILITY)},
		{"mesh_adaptive_search",8,14,4,1,kw_581,0.,0.,0,N_mdm(utype,methodName_MESH_ADAPTIVE_SEARCH)},
		{"mfmc",0,14,4,1,kw_589,0.,0.,9,N_mdm(utype,methodName_MULTIFIDELITY_SAMPLING)},
		{"mlmc",0,17,4,1,kw_603,0.,0.,15,N_mdm(utype,methodName_MULTILEVEL_SAMPLING)},
		{"mlmfmc",0,13,4,1,kw_610,0.,0.,12,N_mdm(utype,methodName_MULTILEVEL_MULTIFIDELITY_SAMPLING)},
		{"moga",8,17,4,1,kw_625,0.,0.,0,N_mdm(utype,methodName_MOGA)},
		{"multi_start",8,7,4,1,kw_629,0.,0.,0,N_mdm(utype,methodName_MULTI_START)},
		{"multidim_parameter_study",8,2,4,1,kw_630,0.,0.,0,N_mdm(utype,methodName_MULTIDIM_PARAMETER_STUDY)},
		{"multifidelity_function_train",8,51,4,1,kw_656,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_FUNCTION_TRAIN)},
		{"multifidelity_mc",0,14,4,1,kw_589,0.,0.,2,N_mdm(utype,methodName_MULTIFIDELITY_SAMPLING)},
		{"multifidelity_polynomial_chaos",8,38,4,1,kw_697,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_POLYNOMIAL_CHAOS)},
		{"multifidelity_sampling",8,14,4,1,kw_589,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_SAMPLING)},
		{"multifidelity_stoch_collocation",8,35,4,1,kw_724,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_STOCH_COLLOCATION)},
		{"multilevel_function_train",8,49,4,1,kw_748,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_FUNCTION_TRAIN)},
		{"multilevel_mc",0,17,4,1,kw_603,0.,0.,4,N_mdm(utype,methodName_MULTILEVEL_SAMPLING)},
		{"multilevel_multifidelity_mc",0,13,4,1,kw_610,0.,0.,1,N_mdm(utype,methodName_MULTILEVEL_MULTIFIDELITY_SAMPLING)},
		{"multilevel_multifidelity_sampling",8,13,4,1,kw_610,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_MULTIFIDELITY_SAMPLING)},
		{"multilevel_polynomial_chaos",8,34,4,1,kw_785,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_POLYNOMIAL_CHAOS)},
		{"multilevel_sampling",8,17,4,1,kw_603,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_SAMPLING)},
		{"ncsu_direct",8,9,4,1,kw_786,0.,0.,0,N_mdm(utype,methodName_NCSU_DIRECT)},
		{"nl2sol",8,15,4,1,kw_787,0.,0.,0,N_mdm(utype,methodName_NL2SOL)},
		{"nlpql_sqp",8,5,4,1,kw_788,0.,0.,0,N_mdm(utype,methodName_NLPQL_SQP)},
		{"nlssol_sqp",8,10,4,1,kw_789,0.,0.,0,N_mdm(utype,methodName_NLSSOL_SQP)},
		{"nond_adaptive_sampling",0,19,4,1,kw_54,0.,0.,-64,N_mdm(utype,methodName_ADAPTIVE_SAMPLING)},
		{"nond_bayes_calibration",0,18,4,1,kw_390,0.,0.,-62,N_mdm(utype,methodName_BAYES_CALIBRATION)},
		{"nond_global_evidence",0,12,4,1,kw_498,0.,0.,-35,N_mdm(utype,methodName_GLOBAL_EVIDENCE)},
		{"nond_global_interval_est",0,11,4,1,kw_519,0.,0.,-35,N_mdm(utype,methodName_GLOBAL_INTERVAL_EST)},
		{"nond_global_reliability",0,22,4,1,kw_537,0.,0.,-35,N_mdm(utype,methodName_GLOBAL_RELIABILITY)},
		{"nond_importance_sampling",0,15,4,1,kw_556,0.,0.,-32,N_mdm(utype,methodName_IMPORTANCE_SAMPLING)},
		{"nond_local_evidence",0,7,4,1,kw_566,0.,0.,-31,N_mdm(utype,methodName_LOCAL_EVIDENCE)},
		{"nond_local_interval_est",0,4,4,1,kw_567,0.,0.,-31,N_mdm(utype,methodName_LOCAL_INTERVAL_EST)},
		{"nond_local_reliability",0,10,4,1,kw_579,0.,0.,-31,N_mdm(utype,methodName_LOCAL_RELIABILITY)},
		{"nond_pof_darts",0,11,4,1,kw_798,0.,0.,16,N_mdm(utype,methodName_POF_DARTS)},
		{"nond_polynomial_chaos",0,37,4,1,kw_843,0.,0.,16,N_mdm(utype,methodName_POLYNOMIAL_CHAOS)},
		{"nond_rkd_darts",0,11,4,1,kw_852,0.,0.,18,N_mdm(utype,methodName_RKD_DARTS)},
		{"nond_sampling",0,19,4,1,kw_867,0.,0.,19,N_mdm(utype,methodName_RANDOM_SAMPLING)},
		{"nond_stoch_collocation",0,32,4,1,kw_891,0.,0.,21,N_mdm(utype,methodName_STOCH_COLLOCATION)},
		{"nonlinear_cg",8,5,4,1,kw_892,0.,0.,0,N_mdm(utype,methodName_NONLINEAR_CG)},
		{"nowpac",8,5,4,1,kw_894,0.,0.,0,N_mdm(utype,methodName_MIT_NOWPAC)},
		{"npsol_sqp",8,10,4,1,kw_895,0.,0.,0,N_mdm(utype,methodName_NPSOL_SQP)},
		{"optpp_cg",8,8,4,1,kw_896,0.,0.,0,N_mdm(utype,methodName_OPTPP_CG)},
		{"optpp_fd_newton",8,12,4,1,kw_899,0.,0.,0,N_mdm(utype,methodName_OPTPP_FD_NEWTON)},
		{"optpp_g_newton",8,12,4,1,kw_902,0.,0.,0,N_mdm(utype,methodName_OPTPP_G_NEWTON)},
		{"optpp_newton",8,12,4,1,kw_905,0.,0.,0,N_mdm(utype,methodName_OPTPP_NEWTON)},
		{"optpp_pds",8,6,4,1,kw_906,0.,0.,0,N_mdm(utype,methodName_OPTPP_PDS)},
		{"optpp_q_newton",8,12,4,1,kw_909,0.,0.,0,N_mdm(utype,methodName_OPTPP_Q_NEWTON)},
		{"output",8,5,2,0,kw_910},
		{"pareto_set",8,10,4,1,kw_914,0.,0.,0,N_mdm(utype,methodName_PARETO_SET)},
		{"pof_darts",8,11,4,1,kw_798,0.,0.,0,N_mdm(utype,methodName_POF_DARTS)},
		{"polynomial_chaos",8,37,4,1,kw_843,0.,0.,0,N_mdm(utype,methodName_POLYNOMIAL_CHAOS)},
		{"psuade_moat",8,4,4,1,kw_915,0.,0.,0,N_mdm(utype,methodName_PSUADE_MOAT)},
		{"richardson_extrap",8,7,4,1,kw_916,0.,0.,0,N_mdm(utype,methodName_RICHARDSON_EXTRAP)},
		{"rkd_darts",8,11,4,1,kw_852,0.,0.,0,N_mdm(utype,methodName_RKD_DARTS)},
		{"rol",8,7,4,1,kw_917,0.,0.,0,N_mdm(utype,methodName_ROL)},
		{"sampling",8,19,4,1,kw_867,0.,0.,0,N_mdm(utype,methodName_RANDOM_SAMPLING)},
		{"snowpac",8,6,4,1,kw_919,0.,0.,0,N_mdm(utype,methodName_MIT_SNOWPAC)},
		{"soga",8,15,4,1,kw_932,0.,0.,0,N_mdm(utype,methodName_SOGA)},
		{"stoch_collocation",8,32,4,1,kw_891,0.,0.,0,N_mdm(utype,methodName_STOCH_COLLOCATION)},
		{"surrogate_based_global",8,8,4,1,kw_933,0.,0.,0,N_mdm(utype,methodName_SURROGATE_BASED_GLOBAL)},
		{"surrogate_based_local",8,16,4,1,kw_939,0.,0.,0,N_mdm(utype,methodName_SURROGATE_BASED_LOCAL)},
		{"surrogate_based_uq",8,21,4,1,kw_956,0.,0.,0,N_mdm(utype,methodName_SURROGATE_BASED_UQ)},
		{"vector_parameter_study",8,4,4,1,kw_957,0.,0.,0,N_mdm(utype,methodName_VECTOR_PARAMETER_STUDY)}
		},
	kw_959[1] = {
		{"refinement_samples",13,0,1,0,0,0.,0.,0,N_mom(ivec,refineSamples)}
		},
	kw_960[3] = {
		{"local_gradient",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_LOCAL_GRAD)},
		{"mean_gradient",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_MEAN_GRAD)},
		{"mean_value",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_MEAN_VALUE)}
		},
	kw_961[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceSampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceSampleType_SUBMETHOD_RANDOM)}
		},
	kw_962[7] = {
		{"decrease",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_DECREASE_TOLERANCE)},
		{"decrease_tolerance",10,0,3,0,0,0.,0.,0,N_mom(Real,decreaseTolerance)},
		{"exhaustive",8,0,5,0,0,0.,0.,0,N_mom(false,subspaceCVIncremental)},
		{"max_rank",9,0,4,0,0,0.,0.,0,N_mom(int,subspaceCVMaxRank)},
		{"minimum",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_MINIMUM_METRIC)},
		{"relative",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_RELATIVE_TOLERANCE)},
		{"relative_tolerance",10,0,2,0,0,0.,0.,0,N_mom(Real,relTolerance)}
		},
	kw_963[1] = {
		{"truncation_tolerance",10,0,1,0,0,0.,0.,0,N_mom(Real,truncationTolerance)}
		},
	kw_964[4] = {
		{"bing_li",8,0,1,0,0,0.,0.,0,N_mom(true,subspaceIdBingLi)},
		{"constantine",8,0,2,0,0,0.,0.,0,N_mom(true,subspaceIdConstantine)},
		{"cross_validation",8,7,4,0,kw_962,0.,0.,0,N_mom(true,subspaceIdCV)},
		{"energy",8,1,3,0,kw_963,0.,0.,0,N_mom(true,subspaceIdEnergy)}
		},
	kw_965[8] = {
		{"actual_model_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"bootstrap_samples",9,0,6,0,0,0.,0.,0,N_mom(int,numReplicates)},
		{"build_surrogate",8,1,7,0,kw_959,0.,0.,0,N_mom(true,subspaceBuildSurrogate)},
		{"dimension",9,0,5,0,0,0.,0.,0,N_mom(int,subspaceDimension)},
		{"initial_samples",9,0,2,0,0,0.,0.,0,N_mom(int,initialSamples)},
		{"normalization",8,3,8,0,kw_960},
		{"sample_type",8,2,3,0,kw_961},
		{"truncation_method",8,4,4,0,kw_964}
		},
	kw_966[1] = {
		{"collocation_ratio",10,0,1,1,0,0.,0.,0,N_mom(Real,adaptedBasisCollocRatio)}
		},
	kw_967[2] = {
		{"ranked",8,0,1,1,0,0.,0.,0,N_mom(type,method_rotation_ROTATION_METHOD_RANKED)},
		{"unranked",8,0,1,1,0,0.,0.,0,N_mom(type,method_rotation_ROTATION_METHOD_UNRANKED)}
		},
	kw_968[4] = {
		{"dimension",9,0,2,0,0,0.,0.,0,N_mom(int,subspaceDimension)},
		{"expansion_order",9,1,1,1,kw_966,0.,0.,0,N_mom(ushint,adaptedBasisExpOrder)},
		{"rotation_method",8,2,3,0,kw_967},
		{"sparse_grid_level",9,0,1,1,0,0.,0.,0,N_mom(ushint,adaptedBasisSparseGridLev)}
		},
	kw_969[2] = {
		{"actual_model_pointer",11,4,1,1,kw_968,0.,0.,0,N_mom(str,actualModelPointer)},
		{"truncation_tolerance",10,0,2,0,0,0.,0.,0,N_mom(Real,adaptedBasisTruncationTolerance)}
		},
	kw_970[1] = {
		{"optional_interface_responses_pointer",11,0,1,0,0,0.,0.,0,N_mom(str,optionalInterfRespPointer)}
		},
	kw_971[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mom(type,subMethodScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mom(type,subMethodScheduling_PEER_SCHEDULING)}
		},
	kw_972[8] = {
		{"identity_response_mapping",8,0,8,0,0,0.,0.,0,N_mom(true,identityRespMap)},
		{"iterator_scheduling",8,2,2,0,kw_971},
		{"iterator_servers",0x19,0,1,0,0,0.,0.,0,N_mom(int,subMethodServers)},
		{"primary_response_mapping",14,0,6,0,0,0.,0.,0,N_mom(RealDL,primaryRespCoeffs)},
		{"primary_variable_mapping",15,0,4,0,0,0.,0.,0,N_mom(strL,primaryVarMaps)},
		{"processors_per_iterator",0x19,0,3,0,0,0.,0.,0,N_mom(int,subMethodProcs)},
		{"secondary_response_mapping",14,0,7,0,0,0.,0.,0,N_mom(RealDL,secondaryRespCoeffs)},
		{"secondary_variable_mapping",15,0,5,0,0,0.,0.,0,N_mom(strL,secondaryVarMaps)}
		},
	kw_973[2] = {
		{"optional_interface_pointer",11,1,1,0,kw_970,0.,0.,0,N_mom(str,interfacePointer)},
		{"sub_method_pointer",11,8,2,1,kw_972,0.,0.,0,N_mom(str,subMethodPointer)}
		},
	kw_974[2] = {
		{"exponential",8,0,1,1,0,0.,0.,0,N_mom(utype,analyticCovIdForm_EXP_L1)},
		{"squared_exponential",8,0,1,1,0,0.,0.,0,N_mom(utype,analyticCovIdForm_EXP_L2)}
		},
	kw_975[3] = {
		{"analytic_covariance",8,2,1,1,kw_974},
		{"dace_method_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,subMethodPointer)},
		{"rf_data_file",11,0,1,1,0,0.,0.,0,N_mom(str,rfDataFileName)}
		},
	kw_976[2] = {
		{"karhunen_loeve",8,0,1,1,0,0.,0.,0,N_mom(utype,randomFieldIdForm_RF_KARHUNEN_LOEVE)},
		{"principal_components",8,0,1,1,0,0.,0.,0,N_mom(utype,randomFieldIdForm_RF_PCA_GP)}
		},
	kw_977[5] = {
		{"build_source",8,3,1,0,kw_975},
		{"expansion_bases",9,0,3,0,0,0.,0.,0,N_mom(int,subspaceDimension)},
		{"expansion_form",8,2,2,0,kw_976},
		{"propagation_model_pointer",11,0,5,1,0,0.,0.,0,N_mom(str,propagationModelPointer)},
		{"truncation_tolerance",10,0,4,0,0,0.,0.,0,N_mom(Real,truncationTolerance)}
		},
	kw_978[1] = {
		{"solution_level_control",11,0,1,0,0,0.,0.,0,N_mom(str,solutionLevelControl)}
		},
	kw_979[2] = {
		{"interface_pointer",11,0,1,0,0,0.,0.,0,N_mom(str,interfacePointer)},
		{"solution_level_cost",14,1,2,0,kw_978,0.,0.,0,N_mom(RealDL,solutionLevelCost)}
		},
	kw_980[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importChalUseVariableLabels)}
		},
	kw_981[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importChalUseVariableLabels)}
		},
	kw_982[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_EVAL_ID)},
		{"header",8,1,1,0,kw_981,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_IFACE_ID)}
		},
	kw_983[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mom(true,importChallengeActive)},
		{"annotated",8,1,1,0,kw_980,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_982,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_NONE)}
		},
	kw_984[6] = {
		{"additive",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_ADDITIVE_CORRECTION)},
		{"combined",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_COMBINED_CORRECTION)},
		{"first_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_1)},
		{"multiplicative",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_MULTIPLICATIVE_CORRECTION)},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_2)},
		{"zeroth_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_0)}
		},
	kw_985[1] = {
		{"folds",0x19,0,1,0,0,0.,0.,0,N_mom(int,refineCVFolds)}
		},
	kw_986[5] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mom(Real,convergenceTolerance)},
		{"cross_validation_metric",11,1,5,0,kw_985,0.,0.,0,N_mom(str,refineCVMetric)},
		{"max_function_evaluations",0x19,0,2,0,0,0.,0.,0,N_mom(sizet,maxFunctionEvals)},
		{"max_iterations",0x19,0,1,0,0,0.,0.,0,N_mom(sizet,maxIterations)},
		{"soft_convergence_limit",0x29,0,4,0,0,0.,0.,0,N_mom(int,softConvergenceLimit)}
		},
	kw_987[1] = {
		{"auto_refinement",8,5,1,0,kw_986,0.,0.,0,N_mom(true,autoRefine)}
		},
	kw_988[2] = {
		{"folds",9,0,1,0,0,0.,0.,0,N_mom(int,numFolds)},
		{"percent",10,0,1,0,0,0.,0.,0,N_mom(Real,percentFold)}
		},
	kw_989[2] = {
		{"cross_validation",8,2,1,0,kw_988,0.,0.,0,N_mom(true,crossValidateFlag)},
		{"press",8,0,2,0,0,0.,0.,0,N_mom(true,pressFlag)}
		},
	kw_990[2] = {
		{"gradient_threshold",10,0,1,1,0,0.,0.,0,N_mom(Real,discontGradThresh)},
		{"jump_threshold",10,0,1,1,0,0.,0.,0,N_mom(Real,discontJumpThresh)}
		},
	kw_991[3] = {
		{"cell_type",11,0,1,0,0,0.,0.,0,N_mom(str,decompCellType)},
		{"discontinuity_detection",8,2,3,0,kw_990,0.,0.,0,N_mom(true,decompDiscontDetect)},
		{"support_layers",9,0,2,0,0,0.,0.,0,N_mom(int,decompSupportLayers)}
		},
	kw_992[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_IFACE_ID)}
		},
	kw_993[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_992,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_NONE)}
		},
	kw_994[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_995[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_994}
		},
	kw_996[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_997[5] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_constant)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_linear)},
		{"none",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_none)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_quadratic)},
		{"reduced_quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_reduced_quadratic)}
		},
	kw_998[8] = {
		{"export_approx_variance_file",11,3,5,0,kw_993,0.,0.,0,N_mom(str,exportApproxVarianceFile)},
		{"export_model",8,2,6,0,kw_995,0.,0.,0,N_mom(true,exportSurrogate)},
		{"find_nugget",9,0,3,0,0,0.,0.,0,N_mom(shint,krigingFindNugget)},
		{"import_model",8,3,7,0,kw_996,0.,0.,0,N_mom(true,importSurrogate)},
		{"nugget",0x1a,0,3,0,0,0.,0.,0,N_mom(Real,krigingNugget)},
		{"num_restarts",0x19,0,2,0,0,1.,0.,0,N_mom(int,numRestarts)},
		{"options_file",11,0,4,0,0,0.,0.,0,N_mom(str,advancedOptionsFilename)},
		{"trend",8,5,1,0,kw_997}
		},
	kw_999[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1000[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_999}
		},
	kw_1001[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1002[4] = {
		{"basis_order",0x29,0,1,1,0,0.,0.,0,N_mom(shint,polynomialOrder)},
		{"export_model",8,2,3,0,kw_1000,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,4,0,kw_1001,0.,0.,0,N_mom(true,importSurrogate)},
		{"options_file",11,0,2,0,0,0.,0.,0,N_mom(str,advancedOptionsFilename)}
		},
	kw_1003[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_1004[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_1003,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_1005[1] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)}
		},
	kw_1006[1] = {
		{"l2_penalty",10,0,1,1,0,0.,0.,0,N_mom(Real,regressionL2Penalty)}
		},
	kw_1007[2] = {
		{"ls",8,0,1,1,0,0.,0.,0,N_mom(type,regressionType_FT_LS)},
		{"rls2",8,1,1,1,kw_1006,0.,0.,0,N_mom(type,regressionType_FT_RLS2)}
		},
	kw_1008[20] = {
		{"adapt_order",8,0,10,0,0,0.,0.,0,N_mom(true,adaptOrder)},
		{"adapt_rank",8,0,15,0,0,0.,0.,0,N_mom(true,adaptRank)},
		{"arithmetic_tolerance",10,0,8,0,0,0.,0.,0,N_mom(Real,statsRoundingTol)},
		{"kick_order",0x19,0,11,0,0,0.,0.,0,N_mom(ushint,kickOrder)},
		{"kick_rank",0x19,0,16,0,0,0.,0.,0,N_mom(sizet,kickRank)},
		{"max_cross_iterations",0x29,0,3,0,0,0.,0.,0,N_mom(int,maxCrossIterations)},
		{"max_cv_order_candidates",0x29,0,13,0,0,0.,0.,0,N_mom(ushint,maxCVOrderCandidates)},
		{"max_cv_rank_candidates",0x29,0,18,0,0,0.,0.,0,N_mom(sizet,maxCVRankCandidates)},
		{"max_order",0x29,0,12,0,0,0.,0.,0,N_mom(ushint,maxOrder)},
		{"max_rank",0x29,0,17,0,0,0.,0.,0,N_mom(sizet,maxRank)},
		{"max_solver_iterations",0x29,0,2,0,0,0.,0.,0,N_mom(sizet,maxSolverIterations)},
		{"order",0x21,1,9,0,kw_1005,0.,0.,6,N_mom(ushint,startOrder)},
		{"rank",0x21,0,14,0,0,0.,0.,6,N_mom(sizet,startRank)},
		{"regression_type",8,2,1,0,kw_1007},
		{"response_scaling",8,0,5,0,0,0.,0.,0,N_mom(true,respScalingFlag)},
		{"rounding_tolerance",10,0,7,0,0,0.,0.,0,N_mom(Real,solverRoundingTol)},
		{"solver_tolerance",10,0,4,0,0,0.,0.,0,N_mom(Real,solverTol)},
		{"start_order",0x29,1,9,0,kw_1005,0.,0.,0,N_mom(ushint,startOrder)},
		{"start_rank",0x29,0,14,0,0,0.,0.,0,N_mom(sizet,startRank)},
		{"tensor_grid",8,0,6,0,0,0.,0.,0,N_mom(true,tensorGridFlag)}
		},
	kw_1009[3] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_constant)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_linear)},
		{"reduced_quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_reduced_quadratic)}
		},
	kw_1010[2] = {
		{"point_selection",8,0,1,0,0,0.,0.,0,N_mom(true,pointSelection)},
		{"trend",8,3,2,0,kw_1009}
		},
	kw_1011[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,exportApproxVarianceFormat_TABULAR_IFACE_ID)}
		},
	kw_1012[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_1011,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxVarianceFormat_TABULAR_NONE)}
		},
	kw_1013[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1014[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_1013}
		},
	kw_1015[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1016[4] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_constant)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_linear)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_quadratic)},
		{"reduced_quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_reduced_quadratic)}
		},
	kw_1017[8] = {
		{"correlation_lengths",14,0,5,0,0,0.,0.,0,N_mom(RealDL,krigingCorrelations)},
		{"export_model",8,2,6,0,kw_1014,0.,0.,0,N_mom(true,exportSurrogate)},
		{"find_nugget",9,0,4,0,0,0.,0.,0,N_mom(shint,krigingFindNugget)},
		{"import_model",8,3,7,0,kw_1015,0.,0.,0,N_mom(true,importSurrogate)},
		{"max_trials",0x19,0,3,0,0,0.,0.,0,N_mom(shint,krigingMaxTrials)},
		{"nugget",0x1a,0,4,0,0,0.,0.,0,N_mom(Real,krigingNugget)},
		{"optimization_method",11,0,2,0,0,0.,0.,0,N_mom(str,krigingOptMethod)},
		{"trend",8,4,1,0,kw_1016}
		},
	kw_1018[3] = {
		{"dakota",8,2,1,1,kw_1010,0.,0.,0,N_mom(lit,surrogateType_global_gaussian)},
		{"export_approx_variance_file",11,3,2,0,kw_1012,0.,0.,0,N_mom(str,exportApproxVarianceFile)},
		{"surfpack",8,8,1,1,kw_1017,0.,0.,0,N_mom(lit,surrogateType_global_kriging)}
		},
	kw_1019[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importUseVariableLabels)}
		},
	kw_1020[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importUseVariableLabels)}
		},
	kw_1021[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,1,1,0,kw_1020,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_1022[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mom(true,importBuildActive)},
		{"annotated",8,1,1,0,kw_1019,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_1021,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_1023[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1024[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_1023}
		},
	kw_1025[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1026[2] = {
		{"cubic",8,0,1,1,0,0.,0.,0,N_mom(lit,marsInterpolation_cubic)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,marsInterpolation_linear)}
		},
	kw_1027[4] = {
		{"export_model",8,2,3,0,kw_1024,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,4,0,kw_1025,0.,0.,0,N_mom(true,importSurrogate)},
		{"interpolation",8,2,2,0,kw_1026},
		{"max_bases",9,0,1,0,0,0.,0.,0,N_mom(shint,marsMaxBases)}
		},
	kw_1028[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1029[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_1028}
		},
	kw_1030[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1031[5] = {
		{"basis_order",0x29,0,1,0,0,0.,0.,0,N_mom(shint,polynomialOrder)},
		{"export_model",8,2,3,0,kw_1029,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,4,0,kw_1030,0.,0.,0,N_mom(true,importSurrogate)},
		{"poly_order",0x21,0,1,0,0,0.,0.,-3,N_mom(shint,polynomialOrder)},
		{"weight_function",9,0,2,0,0,0.,0.,0,N_mom(shint,mlsWeightFunction)}
		},
	kw_1032[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1033[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_1032}
		},
	kw_1034[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1035[6] = {
		{"export_model",8,2,4,0,kw_1033,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,5,0,kw_1034,0.,0.,0,N_mom(true,importSurrogate)},
		{"max_nodes",9,0,1,0,0,0.,0.,0,N_mom(shint,annNodes)},
		{"nodes",1,0,1,0,0,0.,0.,-1,N_mom(shint,annNodes)},
		{"random_weight",9,0,3,0,0,0.,0.,0,N_mom(shint,annRandomWeight)},
		{"range",10,0,2,0,0,0.,0.,0,N_mom(Real,annRange)}
		},
	kw_1036[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1037[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_1036}
		},
	kw_1038[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1039[6] = {
		{"basis_order",0x29,0,1,1,0,0.,0.,0,N_mom(shint,polynomialOrder)},
		{"cubic",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_3)},
		{"export_model",8,2,2,0,kw_1037,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,3,0,kw_1038,0.,0.,0,N_mom(true,importSurrogate)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_1)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_2)}
		},
	kw_1040[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_1041[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_1040}
		},
	kw_1042[3] = {
		{"binary_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_BINARY_ARCHIVE)},
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelImportPrefix)},
		{"text_archive",8,0,2,1,0,0.,0.,0,N_mom(augment_utype,modelImportFormat_TEXT_ARCHIVE)}
		},
	kw_1043[6] = {
		{"bases",9,0,1,0,0,0.,0.,0,N_mom(shint,rbfBases)},
		{"export_model",8,2,5,0,kw_1041,0.,0.,0,N_mom(true,exportSurrogate)},
		{"import_model",8,3,6,0,kw_1042,0.,0.,0,N_mom(true,importSurrogate)},
		{"max_pts",9,0,2,0,0,0.,0.,0,N_mom(shint,rbfMaxPts)},
		{"max_subsets",9,0,4,0,0,0.,0.,0,N_mom(shint,rbfMaxSubsets)},
		{"min_partition",9,0,3,0,0,0.,0.,0,N_mom(shint,rbfMinPartition)}
		},
	kw_1044[3] = {
		{"all",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_all)},
		{"none",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_none)},
		{"region",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_region)}
		},
	kw_1045[30] = {
		{"actual_model_pointer",11,0,4,0,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"challenge_points_file",3,4,11,0,kw_983,0.,0.,12,N_mom(str,importChallengePtsFile)},
		{"correction",8,6,9,0,kw_984},
		{"dace_method_pointer",11,1,4,0,kw_987,0.,0.,0,N_mom(str,subMethodPointer)},
		{"diagnostics",7,2,10,0,kw_989,0.,0.,13,N_mom(strL,diagMetrics)},
		{"domain_decomposition",8,3,2,0,kw_991,0.,0.,0,N_mom(true,domainDecomp)},
		{"experimental_gaussian_process",8,8,1,1,kw_998,0.,0.,0,N_mom(lit,surrogateType_global_exp_gauss_proc)},
		{"experimental_polynomial",8,4,1,1,kw_1002,0.,0.,0,N_mom(lit,surrogateType_global_exp_poly)},
		{"export_approx_points_file",11,3,7,0,kw_1004,0.,0.,0,N_mom(str,exportApproxPtsFile)},
		{"export_points_file",3,3,7,0,kw_1004,0.,0.,-1,N_mom(str,exportApproxPtsFile)},
		{"function_train",8,20,1,1,kw_1008,0.,0.,0,N_mom(lit,surrogateType_global_function_train)},
		{"gaussian_process",8,3,1,1,kw_1018},
		{"import_build_points_file",11,4,6,0,kw_1022,0.,0.,0,N_mom(str,importBuildPtsFile)},
		{"import_challenge_points_file",11,4,11,0,kw_983,0.,0.,0,N_mom(str,importChallengePtsFile)},
		{"import_points_file",3,4,6,0,kw_1022,0.,0.,-2,N_mom(str,importBuildPtsFile)},
		{"kriging",0,3,1,1,kw_1018,0.,0.,-4},
		{"mars",8,4,1,1,kw_1027,0.,0.,0,N_mom(lit,surrogateType_global_mars)},
		{"metrics",15,2,10,0,kw_989,0.,0.,0,N_mom(strL,diagMetrics)},
		{"minimum_points",8,0,3,0,0,0.,0.,0,N_mom(type,pointsManagement_MINIMUM_POINTS)},
		{"moving_least_squares",8,5,1,1,kw_1031,0.,0.,0,N_mom(lit,surrogateType_global_moving_least_squares)},
		{"neural_network",8,6,1,1,kw_1035,0.,0.,0,N_mom(lit,surrogateType_global_neural_network)},
		{"polynomial",8,6,1,1,kw_1039,0.,0.,0,N_mom(lit,surrogateType_global_polynomial)},
		{"radial_basis",8,6,1,1,kw_1043,0.,0.,0,N_mom(lit,surrogateType_global_radial_basis)},
		{"recommended_points",8,0,3,0,0,0.,0.,0,N_mom(type,pointsManagement_RECOMMENDED_POINTS)},
		{"reuse_points",8,3,5,0,kw_1044},
		{"reuse_samples",0,3,5,0,kw_1044,0.,0.,-1},
		{"samples_file",3,4,6,0,kw_1022,0.,0.,-14,N_mom(str,importBuildPtsFile)},
		{"total_points",9,0,3,0,0,0.,0.,0,N_mom(int,pointsTotal)},
		{"truth_model_pointer",3,0,4,0,0,0.,0.,-28,N_mom(str,actualModelPointer)},
		{"use_derivatives",8,0,8,0,0,0.,0.,0,N_mom(true,modelUseDerivsFlag)}
		},
	kw_1046[6] = {
		{"additive",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_ADDITIVE_CORRECTION)},
		{"combined",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_COMBINED_CORRECTION)},
		{"first_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_1)},
		{"multiplicative",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_MULTIPLICATIVE_CORRECTION)},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_2)},
		{"zeroth_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_0)}
		},
	kw_1047[3] = {
		{"correction",8,6,2,0,kw_1046},
		{"model_fidelity_sequence",7,0,1,1,0,0.,0.,1,N_mom(strL,ensembleModelPointers)},
		{"ordered_model_fidelities",15,0,1,1,0,0.,0.,0,N_mom(strL,ensembleModelPointers)}
		},
	kw_1048[3] = {
		{"actual_model_pointer",11,0,2,2,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"taylor_series",8,0,1,1},
		{"truth_model_pointer",3,0,2,2,0,0.,0.,-2,N_mom(str,actualModelPointer)}
		},
	kw_1049[4] = {
		{"actual_model_pointer",3,0,1,1,0,0.,0.,2,N_mom(str,actualModelPointer)},
		{"approximation_models",7,0,2,2,0,0.,0.,2,N_mom(strL,ensembleModelPointers)},
		{"truth_model_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"unordered_model_fidelities",15,0,2,2,0,0.,0.,0,N_mom(strL,ensembleModelPointers)}
		},
	kw_1050[4] = {
		{"actual_model_pointer",11,0,2,2,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"qmea",8,0,1,1,0,0.,0.,0,N_mom(lit,surrogateType_multipoint_qmea)},
		{"tana",8,0,1,1,0,0.,0.,0,N_mom(lit,surrogateType_multipoint_tana)},
		{"truth_model_pointer",3,0,2,2,0,0.,0.,-3,N_mom(str,actualModelPointer)}
		},
	kw_1051[7] = {
		{"global",8,30,2,1,kw_1045},
		{"hierarchical",8,3,2,1,kw_1047,0.,0.,0,N_mom(lit,surrogateType_hierarchical)},
		{"id_surrogates",13,0,1,0,0,0.,0.,0,N_mom(id_to_index_set,surrogateFnIndices)},
		{"local",8,3,2,1,kw_1048,0.,0.,0,N_mom(lit,surrogateType_local_taylor)},
		{"model_ensemble",0,4,2,1,kw_1049,0.,0.,2,N_mom(lit,surrogateType_non_hierarchical)},
		{"multipoint",8,4,2,1,kw_1050},
		{"non_hierarchical",8,4,2,1,kw_1049,0.,0.,0,N_mom(lit,surrogateType_non_hierarchical)}
		},
	kw_1052[12] = {
		{"active_subspace",8,8,2,1,kw_965,0.,0.,0,N_mom(lit,modelType_active_subspace)},
		{"adapted_basis",8,2,2,1,kw_969,0.,0.,0,N_mom(lit,modelType_adapted_basis)},
		{"hierarchical_tagging",8,0,5,0,0,0.,0.,0,N_mom(true,hierarchicalTags)},
		{"id_model",11,0,1,0,0,0.,0.,0,N_mom(str,idModel)},
		{"nested",8,2,2,1,kw_973,0.,0.,0,N_mom(lit,modelType_nested)},
		{"random_field",8,5,2,1,kw_977,0.,0.,0,N_mom(lit,modelType_random_field)},
		{"responses_pointer",11,0,4,0,0,0.,0.,0,N_mom(str,responsesPointer)},
		{"simulation",0,2,2,1,kw_979,0.,0.,1,N_mom(lit,modelType_simulation)},
		{"single",8,2,2,1,kw_979,0.,0.,0,N_mom(lit,modelType_simulation)},
		{"subspace",0,8,2,1,kw_965,0.,0.,-9,N_mom(lit,modelType_active_subspace)},
		{"surrogate",8,7,2,1,kw_1051,0.,0.,0,N_mom(lit,modelType_surrogate)},
		{"variables_pointer",11,0,3,0,0,0.,0.,0,N_mom(str,variablesPointer)}
		},
	kw_1053[2] = {
		{"exp_id",8,0,2,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_HEADER)}
		},
	kw_1054[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_EXPER_ANNOT)},
		{"custom_annotated",8,2,1,0,kw_1053,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)}
		},
	kw_1055[6] = {
		{"experiment_variance_type",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,varianceType)},
		{"interpolate",8,0,5,0,0,0.,0.,0,N_rem(true,interpolateFlag)},
		{"num_config_variables",0x29,0,2,0,0,0.,0.,0,N_rem(sizet,numExpConfigVars)},
		{"num_experiments",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numExperiments)},
		{"scalar_data_file",11,3,4,0,kw_1054,0.,0.,0,N_rem(str,scalarDataFileName)},
		{"variance_type",0x807,0,3,0,0,0.,0.,-5,N_rem(strL,varianceType)}
		},
	kw_1056[2] = {
		{"exp_id",8,0,2,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_HEADER)}
		},
	kw_1057[7] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_EXPER_ANNOT)},
		{"custom_annotated",8,2,1,0,kw_1056,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"experiment_variance_type",0x80f,0,4,0,0,0.,0.,0,N_rem(strL,varianceType)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"num_config_variables",0x29,0,3,0,0,0.,0.,0,N_rem(sizet,numExpConfigVars)},
		{"num_experiments",0x29,0,2,0,0,0.,0.,0,N_rem(sizet,numExperiments)},
		{"variance_type",0x807,0,4,0,0,0.,0.,-4,N_rem(strL,varianceType)}
		},
	kw_1058[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_1059[6] = {
		{"nonlinear_equality_scale_types",0x807,0,2,0,0,0.,0.,3,N_rem(strL,nonlinearEqScaleTypes)},
		{"nonlinear_equality_scales",0x806,0,3,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqScales)},
		{"nonlinear_equality_targets",6,0,1,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqTargets)},
		{"scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,nonlinearEqScaleTypes)},
		{"scales",0x80e,0,3,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqScales)},
		{"targets",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqTargets)}
		},
	kw_1060[8] = {
		{"lower_bounds",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_lower_bounds",6,0,1,0,0,0.,0.,-1,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_scale_types",0x807,0,3,0,0,0.,0.,3,N_rem(strL,nonlinearIneqScaleTypes)},
		{"nonlinear_inequality_scales",0x806,0,4,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqScales)},
		{"nonlinear_inequality_upper_bounds",6,0,2,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqUpperBnds)},
		{"scale_types",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,nonlinearIneqScaleTypes)},
		{"scales",0x80e,0,4,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqScales)},
		{"upper_bounds",14,0,2,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqUpperBnds)}
		},
	kw_1061[16] = {
		{"calibration_data",8,6,5,0,kw_1055,0.,0.,0,N_rem(true,calibrationDataFlag)},
		{"calibration_data_file",11,7,5,0,kw_1057,0.,0.,0,N_rem(str,scalarDataFileName)},
		{"calibration_term_scales",6,0,3,0,0,0.,0.,10,N_rem(RealDL,primaryRespFnScales)},
		{"calibration_weights",6,0,4,0,0,0.,0.,12,N_rem(RealDL,primaryRespFnWeights)},
		{"field_calibration_terms",0x29,3,2,0,kw_1058,0.,0.,0,N_rem(sizet,numFieldLeastSqTerms)},
		{"least_squares_data_file",3,7,5,0,kw_1057,0.,0.,-4,N_rem(str,scalarDataFileName)},
		{"least_squares_term_scales",6,0,3,0,0,0.,0.,6,N_rem(RealDL,primaryRespFnScales)},
		{"least_squares_weights",6,0,4,0,0,0.,0.,8,N_rem(RealDL,primaryRespFnWeights)},
		{"nonlinear_equality_constraints",0x29,6,8,0,kw_1059,0.,0.,0,N_rem(sizet,numNonlinearEqConstraints)},
		{"nonlinear_inequality_constraints",0x29,8,7,0,kw_1060,0.,0.,0,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_nonlinear_equality_constraints",0x21,6,8,0,kw_1059,0.,0.,-2,N_rem(sizet,numNonlinearEqConstraints)},
		{"num_nonlinear_inequality_constraints",0x21,8,7,0,kw_1060,0.,0.,-2,N_rem(sizet,numNonlinearIneqConstraints)},
		{"primary_scales",14,0,3,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnScales)},
		{"scalar_calibration_terms",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numScalarLeastSqTerms)},
		{"simulation_variance",0x80e,0,6,0,0,0.,0.,0,N_rem(RealL,simVariance)},
		{"weights",14,0,4,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnWeights)}
		},
	kw_1062[4] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_bounds)},
		{"ignore_bounds",8,0,1,0,0,0.,0.,0,N_rem(true,ignoreBounds)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_relative)}
		},
	kw_1063[10] = {
		{"central",8,0,6,0,0,0.,0.,0,N_rem(lit,intervalType_central)},
		{"dakota",8,4,4,0,kw_1062,0.,0.,0,N_rem(lit,methodSource_dakota)},
		{"fd_gradient_step_size",6,0,7,0,0,0.,0.,1,N_rem(RealL,fdGradStepSize)},
		{"fd_step_size",14,0,7,0,0,0.,0.,0,N_rem(RealL,fdGradStepSize)},
		{"forward",8,0,6,0,0,0.,0.,0,N_rem(lit,intervalType_forward)},
		{"id_analytic_gradients",13,0,2,2,0,0.,0.,0,N_rem(intset,idAnalyticGrads)},
		{"id_numerical_gradients",13,0,1,1,0,0.,0.,0,N_rem(intset,idNumericalGrads)},
		{"interval_type",8,0,5},
		{"method_source",8,0,3},
		{"vendor",8,0,4,0,0,0.,0.,0,N_rem(lit,methodSource_vendor)}
		},
	kw_1064[2] = {
		{"fd_hessian_step_size",6,0,1,0,0,0.,0.,1,N_rem(RealL,fdHessStepSize)},
		{"fd_step_size",14,0,1,0,0,0.,0.,0,N_rem(RealL,fdHessStepSize)}
		},
	kw_1065[1] = {
		{"damped",8,0,1,0,0,0.,0.,0,N_rem(lit,quasiHessianType_damped_bfgs)}
		},
	kw_1066[2] = {
		{"bfgs",8,1,1,1,kw_1065,0.,0.,0,N_rem(lit,quasiHessianType_bfgs)},
		{"sr1",8,0,1,1,0,0.,0.,0,N_rem(lit,quasiHessianType_sr1)}
		},
	kw_1067[8] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_bounds)},
		{"central",8,0,3,0,0,0.,0.,0,N_rem(true,centralHess)},
		{"forward",8,0,3,0,0,0.,0.,0,N_rem(false,centralHess)},
		{"id_analytic_hessians",13,0,5,0,0,0.,0.,0,N_rem(intset,idAnalyticHessians)},
		{"id_numerical_hessians",13,2,1,0,kw_1064,0.,0.,0,N_rem(intset,idNumericalHessians)},
		{"id_quasi_hessians",13,2,4,0,kw_1066,0.,0.,0,N_rem(intset,idQuasiHessians)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_relative)}
		},
	kw_1068[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_1069[6] = {
		{"nonlinear_equality_scale_types",0x807,0,2,0,0,0.,0.,3,N_rem(strL,nonlinearEqScaleTypes)},
		{"nonlinear_equality_scales",0x806,0,3,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqScales)},
		{"nonlinear_equality_targets",6,0,1,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqTargets)},
		{"scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,nonlinearEqScaleTypes)},
		{"scales",0x80e,0,3,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqScales)},
		{"targets",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqTargets)}
		},
	kw_1070[8] = {
		{"lower_bounds",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_lower_bounds",6,0,1,0,0,0.,0.,-1,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_scale_types",0x807,0,3,0,0,0.,0.,3,N_rem(strL,nonlinearIneqScaleTypes)},
		{"nonlinear_inequality_scales",0x806,0,4,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqScales)},
		{"nonlinear_inequality_upper_bounds",6,0,2,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqUpperBnds)},
		{"scale_types",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,nonlinearIneqScaleTypes)},
		{"scales",0x80e,0,4,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqScales)},
		{"upper_bounds",14,0,2,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqUpperBnds)}
		},
	kw_1071[15] = {
		{"field_objectives",0x29,3,8,0,kw_1068,0.,0.,0,N_rem(sizet,numFieldObjectiveFunctions)},
		{"multi_objective_weights",6,0,4,0,0,0.,0.,13,N_rem(RealDL,primaryRespFnWeights)},
		{"nonlinear_equality_constraints",0x29,6,6,0,kw_1069,0.,0.,0,N_rem(sizet,numNonlinearEqConstraints)},
		{"nonlinear_inequality_constraints",0x29,8,5,0,kw_1070,0.,0.,0,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_field_objectives",0x21,3,8,0,kw_1068,0.,0.,-4,N_rem(sizet,numFieldObjectiveFunctions)},
		{"num_nonlinear_equality_constraints",0x21,6,6,0,kw_1069,0.,0.,-3,N_rem(sizet,numNonlinearEqConstraints)},
		{"num_nonlinear_inequality_constraints",0x21,8,5,0,kw_1070,0.,0.,-3,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_scalar_objectives",0x21,0,7,0,0,0.,0.,5,N_rem(sizet,numScalarObjectiveFunctions)},
		{"objective_function_scale_types",0x807,0,2,0,0,0.,0.,2,N_rem(strL,primaryRespFnScaleTypes)},
		{"objective_function_scales",6,0,3,0,0,0.,0.,2,N_rem(RealDL,primaryRespFnScales)},
		{"primary_scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,primaryRespFnScaleTypes)},
		{"primary_scales",14,0,3,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnScales)},
		{"scalar_objectives",0x29,0,7,0,0,0.,0.,0,N_rem(sizet,numScalarObjectiveFunctions)},
		{"sense",0x80f,0,1,0,0,0.,0.,0,N_rem(strL,primaryRespFnSense)},
		{"weights",14,0,4,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnWeights)}
		},
	kw_1072[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_1073[4] = {
		{"field_responses",0x29,3,2,0,kw_1072,0.,0.,0,N_rem(sizet,numFieldResponseFunctions)},
		{"num_field_responses",0x21,3,2,0,kw_1072,0.,0.,-1,N_rem(sizet,numFieldResponseFunctions)},
		{"num_scalar_responses",0x21,0,1,0,0,0.,0.,1,N_rem(sizet,numScalarResponseFunctions)},
		{"scalar_responses",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numScalarResponseFunctions)}
		},
	kw_1074[4] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_bounds)},
		{"ignore_bounds",8,0,1,0,0,0.,0.,0,N_rem(true,ignoreBounds)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_relative)}
		},
	kw_1075[8] = {
		{"central",8,0,4,0,0,0.,0.,0,N_rem(lit,intervalType_central)},
		{"dakota",8,4,2,0,kw_1074,0.,0.,0,N_rem(lit,methodSource_dakota)},
		{"fd_gradient_step_size",6,0,5,0,0,0.,0.,1,N_rem(RealL,fdGradStepSize)},
		{"fd_step_size",14,0,5,0,0,0.,0.,0,N_rem(RealL,fdGradStepSize)},
		{"forward",8,0,4,0,0,0.,0.,0,N_rem(lit,intervalType_forward)},
		{"interval_type",8,0,3},
		{"method_source",8,0,1},
		{"vendor",8,0,2,0,0,0.,0.,0,N_rem(lit,methodSource_vendor)}
		},
	kw_1076[7] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_bounds)},
		{"central",8,0,3,0,0,0.,0.,0,N_rem(true,centralHess)},
		{"fd_hessian_step_size",6,0,1,0,0,0.,0.,1,N_rem(RealL,fdHessStepSize)},
		{"fd_step_size",14,0,1,0,0,0.,0.,0,N_rem(RealL,fdHessStepSize)},
		{"forward",8,0,3,0,0,0.,0.,0,N_rem(false,centralHess)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_relative)}
		},
	kw_1077[1] = {
		{"damped",8,0,1,0,0,0.,0.,0,N_rem(lit,quasiHessianType_damped_bfgs)}
		},
	kw_1078[2] = {
		{"bfgs",8,1,1,1,kw_1077,0.,0.,0,N_rem(lit,quasiHessianType_bfgs)},
		{"sr1",8,0,1,1,0,0.,0.,0,N_rem(lit,quasiHessianType_sr1)}
		},
	kw_1079[19] = {
		{"analytic_gradients",8,0,4,2,0,0.,0.,0,N_rem(lit,gradientType_analytic)},
		{"analytic_hessians",8,0,5,3,0,0.,0.,0,N_rem(lit,hessianType_analytic)},
		{"calibration_terms",0x29,16,3,1,kw_1061,0.,0.,0,N_rem(sizet,numLeastSqTerms)},
		{"descriptors",15,0,2,0,0,0.,0.,0,N_rem(strL,responseLabels)},
		{"id_responses",11,0,1,0,0,0.,0.,0,N_rem(str,idResponses)},
		{"least_squares_terms",0x21,16,3,1,kw_1061,0.,0.,-3,N_rem(sizet,numLeastSqTerms)},
		{"mixed_gradients",8,10,4,2,kw_1063,0.,0.,0,N_rem(lit,gradientType_mixed)},
		{"mixed_hessians",8,8,5,3,kw_1067,0.,0.,0,N_rem(lit,hessianType_mixed)},
		{"no_gradients",8,0,4,2,0,0.,0.,0,N_rem(lit,gradientType_none)},
		{"no_hessians",8,0,5,3,0,0.,0.,0,N_rem(lit,hessianType_none)},
		{"num_least_squares_terms",0x21,16,3,1,kw_1061,0.,0.,-8,N_rem(sizet,numLeastSqTerms)},
		{"num_objective_functions",0x21,15,3,1,kw_1071,0.,0.,4,N_rem(sizet,numObjectiveFunctions)},
		{"num_response_functions",0x21,4,3,1,kw_1073,0.,0.,6,N_rem(sizet,numResponseFunctions)},
		{"numerical_gradients",8,8,4,2,kw_1075,0.,0.,0,N_rem(lit,gradientType_numerical)},
		{"numerical_hessians",8,7,5,3,kw_1076,0.,0.,0,N_rem(lit,hessianType_numerical)},
		{"objective_functions",0x29,15,3,1,kw_1071,0.,0.,0,N_rem(sizet,numObjectiveFunctions)},
		{"quasi_hessians",8,2,5,3,kw_1078,0.,0.,0,N_rem(lit,hessianType_quasi)},
		{"response_descriptors",7,0,2,0,0,0.,0.,-14,N_rem(strL,responseLabels)},
		{"response_functions",0x29,4,3,1,kw_1073,0.,0.,0,N_rem(sizet,numResponseFunctions)}
		},
	kw_1080[6] = {
		{"aleatory",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_ALEATORY_UNCERTAIN_VIEW)},
		{"all",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_ALL_VIEW)},
		{"design",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_DESIGN_VIEW)},
		{"epistemic",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_EPISTEMIC_UNCERTAIN_VIEW)},
		{"state",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_STATE_VIEW)},
		{"uncertain",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_UNCERTAIN_VIEW)}
		},
	kw_1081[11] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,betaUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(RealLb,betaUncBetas)},
		{"buv_alphas",6,0,1,1,0,0.,0.,-2,N_vam(RealLb,betaUncAlphas)},
		{"buv_betas",6,0,2,2,0,0.,0.,-2,N_vam(RealLb,betaUncBetas)},
		{"buv_descriptors",7,0,6,0,0,0.,0.,3,N_vae(caulbl,CAUVar_beta)},
		{"buv_lower_bounds",6,0,3,3,0,0.,0.,4,N_vam(rvec,betaUncLowerBnds)},
		{"buv_upper_bounds",6,0,4,4,0,0.,0.,4,N_vam(rvec,betaUncUpperBnds)},
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(caulbl,CAUVar_beta)},
		{"initial_point",14,0,5,0,0,0.,0.,0,N_vam(rvec,betaUncVars)},
		{"lower_bounds",14,0,3,3,0,0.,0.,0,N_vam(rvec,betaUncLowerBnds)},
		{"upper_bounds",14,0,4,4,0,0.,0.,0,N_vam(rvec,betaUncUpperBnds)}
		},
	kw_1082[5] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_binomial)},
		{"initial_point",13,0,3,0,0,0.,0.,0,N_vam(IntLb,binomialUncVars)},
		{"num_trials",13,0,2,2,0,0.,0.,0,N_vam(IntLb,binomialUncNumTrials)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,binomialUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,binomialUncProbPerTrial)}
		},
	kw_1083[12] = {
		{"cdv_descriptors",7,0,6,0,0,0.,0.,6,N_vam(strL,continuousDesignLabels)},
		{"cdv_initial_point",6,0,1,0,0,0.,0.,6,N_vam(rvec,continuousDesignVars)},
		{"cdv_lower_bounds",6,0,2,0,0,0.,0.,6,N_vam(rvec,continuousDesignLowerBnds)},
		{"cdv_scale_types",0x807,0,4,0,0,0.,0.,6,N_vam(strL,continuousDesignScaleTypes)},
		{"cdv_scales",0x806,0,5,0,0,0.,0.,6,N_vam(rvec,continuousDesignScales)},
		{"cdv_upper_bounds",6,0,3,0,0,0.,0.,6,N_vam(rvec,continuousDesignUpperBnds)},
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vam(strL,continuousDesignLabels)},
		{"initial_point",14,0,1,0,0,0.,0.,0,N_vam(rvec,continuousDesignVars)},
		{"lower_bounds",14,0,2,0,0,0.,0.,0,N_vam(rvec,continuousDesignLowerBnds)},
		{"scale_types",0x80f,0,4,0,0,0.,0.,0,N_vam(strL,continuousDesignScaleTypes)},
		{"scales",0x80e,0,5,0,0,0.,0.,0,N_vam(rvec,continuousDesignScales)},
		{"upper_bounds",14,0,3,0,0,0.,0.,0,N_vam(rvec,continuousDesignUpperBnds)}
		},
	kw_1084[10] = {
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(ceulbl,CEUVar_interval)},
		{"initial_point",14,0,5,0,0,0.,0.,0,N_vam(rvec,continuousIntervalUncVars)},
		{"interval_probabilities",14,0,2,0,0,0.,0.,0,N_vam(newrvec,Var_Info_CIp)},
		{"interval_probs",6,0,2,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_CIp)},
		{"iuv_descriptors",7,0,6,0,0,0.,0.,-4,N_vae(ceulbl,CEUVar_interval)},
		{"iuv_interval_probs",6,0,2,0,0,0.,0.,-3,N_vam(newrvec,Var_Info_CIp)},
		{"iuv_num_intervals",5,0,1,0,0,0.,0.,2,N_vam(newiarray,Var_Info_nCI)},
		{"lower_bounds",14,0,3,1,0,0.,0.,0,N_vam(newrvec,Var_Info_CIlb)},
		{"num_intervals",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nCI)},
		{"upper_bounds",14,0,4,2,0,0.,0.,0,N_vam(newrvec,Var_Info_CIub)}
		},
	kw_1085[8] = {
		{"csv_descriptors",7,0,4,0,0,0.,0.,4,N_vam(strL,continuousStateLabels)},
		{"csv_initial_state",6,0,1,0,0,0.,0.,4,N_vam(rvec,continuousStateVars)},
		{"csv_lower_bounds",6,0,2,0,0,0.,0.,4,N_vam(rvec,continuousStateLowerBnds)},
		{"csv_upper_bounds",6,0,3,0,0,0.,0.,4,N_vam(rvec,continuousStateUpperBnds)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,continuousStateLabels)},
		{"initial_state",14,0,1,0,0,0.,0.,0,N_vam(rvec,continuousStateVars)},
		{"lower_bounds",14,0,2,0,0,0.,0.,0,N_vam(rvec,continuousStateLowerBnds)},
		{"upper_bounds",14,0,3,0,0,0.,0.,0,N_vam(rvec,continuousStateUpperBnds)}
		},
	kw_1086[8] = {
		{"ddv_descriptors",7,0,4,0,0,0.,0.,4,N_vam(strL,discreteDesignRangeLabels)},
		{"ddv_initial_point",5,0,1,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeVars)},
		{"ddv_lower_bounds",5,0,2,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeLowerBnds)},
		{"ddv_upper_bounds",5,0,3,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeUpperBnds)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteDesignRangeLabels)},
		{"initial_point",13,0,1,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeVars)},
		{"lower_bounds",13,0,2,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeLowerBnds)},
		{"upper_bounds",13,0,3,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeUpperBnds)}
		},
	kw_1087[1] = {
		{"adjacency_matrix",13,0,1,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddsia)}
		},
	kw_1088[7] = {
		{"categorical",15,1,3,0,kw_1087,0.,0.,0,N_vam(categorical,discreteDesignSetIntCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetIntLabels)},
		{"elements",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_ddsi)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddsi)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(ivec,discreteDesignSetIntVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddsi)},
		{"set_values",5,0,2,1,0,0.,0.,-4,N_vam(newivec,Var_Info_ddsi)}
		},
	kw_1089[1] = {
		{"adjacency_matrix",13,0,1,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddsra)}
		},
	kw_1090[7] = {
		{"categorical",15,1,3,0,kw_1089,0.,0.,0,N_vam(categorical,discreteDesignSetRealCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetRealLabels)},
		{"elements",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_ddsr)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddsr)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,discreteDesignSetRealVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddsr)},
		{"set_values",6,0,2,1,0,0.,0.,-4,N_vam(newrvec,Var_Info_ddsr)}
		},
	kw_1091[7] = {
		{"adjacency_matrix",13,0,3,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddssa)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetStrLabels)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_ddss)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddss)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteDesignSetStrVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddss)},
		{"set_values",7,0,2,1,0,0.,0.,-4,N_vam(newsarray,Var_Info_ddss)}
		},
	kw_1092[3] = {
		{"integer",0x19,7,1,0,kw_1088,0.,0.,0,N_vam(sizet,numDiscreteDesSetIntVars)},
		{"real",0x19,7,3,0,kw_1090,0.,0.,0,N_vam(sizet,numDiscreteDesSetRealVars)},
		{"string",0x19,7,2,0,kw_1091,0.,0.,0,N_vam(sizet,numDiscreteDesSetStrVars)}
		},
	kw_1093[9] = {
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(deuilbl,DEUIVar_interval)},
		{"initial_point",13,0,5,0,0,0.,0.,0,N_vam(ivec,discreteIntervalUncVars)},
		{"interval_probabilities",14,0,2,0,0,0.,0.,0,N_vam(newrvec,Var_Info_DIp)},
		{"interval_probs",6,0,2,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_DIp)},
		{"lower_bounds",13,0,3,1,0,0.,0.,0,N_vam(newivec,Var_Info_DIlb)},
		{"num_intervals",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nDI)},
		{"range_probabilities",6,0,2,0,0,0.,0.,-4,N_vam(newrvec,Var_Info_DIp)},
		{"range_probs",6,0,2,0,0,0.,0.,-5,N_vam(newrvec,Var_Info_DIp)},
		{"upper_bounds",13,0,4,2,0,0.,0.,0,N_vam(newivec,Var_Info_DIub)}
		},
	kw_1094[8] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteStateRangeLabels)},
		{"dsv_descriptors",7,0,4,0,0,0.,0.,-1,N_vam(strL,discreteStateRangeLabels)},
		{"dsv_initial_state",5,0,1,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeVars)},
		{"dsv_lower_bounds",5,0,2,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeLowerBnds)},
		{"dsv_upper_bounds",5,0,3,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeUpperBnds)},
		{"initial_state",13,0,1,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeVars)},
		{"lower_bounds",13,0,2,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeLowerBnds)},
		{"upper_bounds",13,0,3,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeUpperBnds)}
		},
	kw_1095[7] = {
		{"categorical",15,0,3,0,0,0.,0.,0,N_vam(categorical,discreteStateSetIntCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteStateSetIntLabels)},
		{"elements",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_dssi)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndssi)},
		{"initial_state",13,0,4,0,0,0.,0.,0,N_vam(ivec,discreteStateSetIntVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndssi)},
		{"set_values",5,0,2,1,0,0.,0.,-4,N_vam(newivec,Var_Info_dssi)}
		},
	kw_1096[7] = {
		{"categorical",15,0,3,0,0,0.,0.,0,N_vam(categorical,discreteStateSetRealCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteStateSetRealLabels)},
		{"elements",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_dssr)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndssr)},
		{"initial_state",14,0,4,0,0,0.,0.,0,N_vam(rvec,discreteStateSetRealVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndssr)},
		{"set_values",6,0,2,1,0,0.,0.,-4,N_vam(newrvec,Var_Info_dssr)}
		},
	kw_1097[6] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteStateSetStrLabels)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_dsss)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndsss)},
		{"initial_state",15,0,3,0,0,0.,0.,0,N_vam(strL,discreteStateSetStrVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndsss)},
		{"set_values",7,0,2,1,0,0.,0.,-4,N_vam(newsarray,Var_Info_dsss)}
		},
	kw_1098[3] = {
		{"integer",0x19,7,1,0,kw_1095,0.,0.,0,N_vam(sizet,numDiscreteStateSetIntVars)},
		{"real",0x19,7,3,0,kw_1096,0.,0.,0,N_vam(sizet,numDiscreteStateSetRealVars)},
		{"string",0x19,6,2,0,kw_1097,0.,0.,0,N_vam(sizet,numDiscreteStateSetStrVars)}
		},
	kw_1099[9] = {
		{"categorical",15,0,4,0,0,0.,0.,0,N_vam(categorical,discreteUncSetIntCat)},
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(deuilbl,DEUIVar_set_int)},
		{"elements",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_dusi)},
		{"elements_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndusi)},
		{"initial_point",13,0,5,0,0,0.,0.,0,N_vam(ivec,discreteUncSetIntVars)},
		{"num_set_values",5,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndusi)},
		{"set_probabilities",14,0,3,0,0,0.,0.,0,N_vam(newrvec,Var_Info_DSIp)},
		{"set_probs",6,0,3,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_DSIp)},
		{"set_values",5,0,2,1,0,0.,0.,-6,N_vam(newivec,Var_Info_dusi)}
		},
	kw_1100[9] = {
		{"categorical",15,0,4,0,0,0.,0.,0,N_vam(categorical,discreteUncSetRealCat)},
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(deurlbl,DEURVar_set_real)},
		{"elements",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_dusr)},
		{"elements_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndusr)},
		{"initial_point",14,0,5,0,0,0.,0.,0,N_vam(rvec,discreteUncSetRealVars)},
		{"num_set_values",5,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndusr)},
		{"set_probabilities",14,0,3,0,0,0.,0.,0,N_vam(newrvec,Var_Info_DSRp)},
		{"set_probs",6,0,3,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_DSRp)},
		{"set_values",6,0,2,1,0,0.,0.,-6,N_vam(newrvec,Var_Info_dusr)}
		},
	kw_1101[8] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(deuslbl,DEUSVar_set_str)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_duss)},
		{"elements_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nduss)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteUncSetStrVars)},
		{"num_set_values",5,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nduss)},
		{"set_probabilities",14,0,3,0,0,0.,0.,0,N_vam(newrvec,Var_Info_DSSp)},
		{"set_probs",6,0,3,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_DSSp)},
		{"set_values",7,0,2,1,0,0.,0.,-6,N_vam(newsarray,Var_Info_duss)}
		},
	kw_1102[3] = {
		{"integer",0x19,9,1,0,kw_1099,0.,0.,0,N_vam(sizet,numDiscreteUncSetIntVars)},
		{"real",0x19,9,3,0,kw_1100,0.,0.,0,N_vam(sizet,numDiscreteUncSetRealVars)},
		{"string",0x19,8,2,0,kw_1101,0.,0.,0,N_vam(sizet,numDiscreteUncSetStrVars)}
		},
	kw_1103[5] = {
		{"betas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,exponentialUncBetas)},
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(caulbl,CAUVar_exponential)},
		{"euv_betas",6,0,1,1,0,0.,0.,-2,N_vam(RealLb,exponentialUncBetas)},
		{"euv_descriptors",7,0,3,0,0,0.,0.,-2,N_vae(caulbl,CAUVar_exponential)},
		{"initial_point",14,0,2,0,0,0.,0.,0,N_vam(RealLb,exponentialUncVars)}
		},
	kw_1104[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,frechetUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(rvec,frechetUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_frechet)},
		{"fuv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,frechetUncAlphas)},
		{"fuv_betas",6,0,2,2,0,0.,0.,-3,N_vam(rvec,frechetUncBetas)},
		{"fuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_frechet)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,frechetUncVars)}
		},
	kw_1105[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,gammaUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(RealLb,gammaUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_gamma)},
		{"gauv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,gammaUncAlphas)},
		{"gauv_betas",6,0,2,2,0,0.,0.,-3,N_vam(RealLb,gammaUncBetas)},
		{"gauv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_gamma)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,gammaUncVars)}
		},
	kw_1106[4] = {
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_geometric)},
		{"initial_point",13,0,2,0,0,0.,0.,0,N_vam(IntLb,geometricUncVars)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,geometricUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,geometricUncProbPerTrial)}
		},
	kw_1107[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,gumbelUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(rvec,gumbelUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_gumbel)},
		{"guuv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,gumbelUncAlphas)},
		{"guuv_betas",6,0,2,2,0,0.,0.,-3,N_vam(rvec,gumbelUncBetas)},
		{"guuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_gumbel)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,gumbelUncVars)}
		},
	kw_1108[11] = {
		{"abscissas",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_hba)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hbc)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(caulbl,CAUVar_histogram_bin)},
		{"huv_bin_abscissas",6,0,2,1,0,0.,0.,-3,N_vam(newrvec,Var_Info_hba)},
		{"huv_bin_counts",6,0,3,2,0,0.,0.,-3,N_vam(newrvec,Var_Info_hbc)},
		{"huv_bin_descriptors",7,0,5,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_histogram_bin)},
		{"huv_bin_ordinates",6,0,3,2,0,0.,0.,3,N_vam(newrvec,Var_Info_hbo)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,histogramBinUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,2,N_vam(newiarray,Var_Info_nhbp)},
		{"ordinates",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hbo)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhbp)}
		},
	kw_1109[6] = {
		{"abscissas",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_hpia)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hpic)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_histogram_point_int)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(ivec,histogramPointIntUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhpip)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhpip)}
		},
	kw_1110[6] = {
		{"abscissas",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_hpra)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hprc)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(daurlbl,DAURVar_histogram_point_real)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,histogramPointRealUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhprp)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhprp)}
		},
	kw_1111[6] = {
		{"abscissas",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_hpsa)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hpsc)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauslbl,DAUSVar_histogram_point_str)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,histogramPointStrUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhpsp)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhpsp)}
		},
	kw_1112[3] = {
		{"integer",0x19,6,1,0,kw_1109,0.,0.,0,N_vam(sizet,numHistogramPtIntUncVars)},
		{"real",0x19,6,3,0,kw_1110,0.,0.,0,N_vam(sizet,numHistogramPtRealUncVars)},
		{"string",0x19,6,2,0,kw_1111,0.,0.,0,N_vam(sizet,numHistogramPtStrUncVars)}
		},
	kw_1113[5] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_hypergeometric)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(IntLb,hyperGeomUncVars)},
		{"num_drawn",13,0,3,3,0,0.,0.,0,N_vam(IntLb,hyperGeomUncNumDrawn)},
		{"selected_population",13,0,2,2,0,0.,0.,0,N_vam(IntLb,hyperGeomUncSelectedPop)},
		{"total_population",13,0,1,1,0,0.,0.,0,N_vam(IntLb,hyperGeomUncTotalPop)}
		},
	kw_1114[2] = {
		{"lnuv_zetas",6,0,1,1,0,0.,0.,1,N_vam(RealLb,lognormalUncZetas)},
		{"zetas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncZetas)}
		},
	kw_1115[4] = {
		{"error_factors",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncErrFacts)},
		{"lnuv_error_factors",6,0,1,1,0,0.,0.,-1,N_vam(RealLb,lognormalUncErrFacts)},
		{"lnuv_std_deviations",6,0,1,1,0,0.,0.,1,N_vam(RealLb,lognormalUncStdDevs)},
		{"std_deviations",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncStdDevs)}
		},
	kw_1116[11] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(caulbl,CAUVar_lognormal)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(RealLb,lognormalUncVars)},
		{"lambdas",14,2,1,1,kw_1114,0.,0.,0,N_vam(rvec,lognormalUncLambdas)},
		{"lnuv_descriptors",7,0,5,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_lognormal)},
		{"lnuv_lambdas",6,2,1,1,kw_1114,0.,0.,-2,N_vam(rvec,lognormalUncLambdas)},
		{"lnuv_lower_bounds",6,0,2,0,0,0.,0.,3,N_vam(RealLb,lognormalUncLowerBnds)},
		{"lnuv_means",6,4,1,1,kw_1115,0.,0.,3,N_vam(RealLb,lognormalUncMeans)},
		{"lnuv_upper_bounds",6,0,3,0,0,0.,0.,3,N_vam(RealUb,lognormalUncUpperBnds)},
		{"lower_bounds",14,0,2,0,0,0.,0.,0,N_vam(RealLb,lognormalUncLowerBnds)},
		{"means",14,4,1,1,kw_1115,0.,0.,0,N_vam(RealLb,lognormalUncMeans)},
		{"upper_bounds",14,0,3,0,0,0.,0.,0,N_vam(RealUb,lognormalUncUpperBnds)}
		},
	kw_1117[7] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_loguniform)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,loguniformUncVars)},
		{"lower_bounds",14,0,1,1,0,0.,0.,0,N_vam(RealLb,loguniformUncLowerBnds)},
		{"luuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_loguniform)},
		{"luuv_lower_bounds",6,0,1,1,0,0.,0.,-2,N_vam(RealLb,loguniformUncLowerBnds)},
		{"luuv_upper_bounds",6,0,2,2,0,0.,0.,1,N_vam(RealUb,loguniformUncUpperBnds)},
		{"upper_bounds",14,0,2,2,0,0.,0.,0,N_vam(RealUb,loguniformUncUpperBnds)}
		},
	kw_1118[5] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_negative_binomial)},
		{"initial_point",13,0,3,0,0,0.,0.,0,N_vam(IntLb,negBinomialUncVars)},
		{"num_trials",13,0,2,2,0,0.,0.,0,N_vam(IntLb,negBinomialUncNumTrials)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,negBinomialUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,negBinomialUncProbPerTrial)}
		},
	kw_1119[11] = {
		{"descriptors",15,0,6,0,0,0.,0.,0,N_vae(caulbl,CAUVar_normal)},
		{"initial_point",14,0,5,0,0,0.,0.,0,N_vam(rvec,normalUncVars)},
		{"lower_bounds",14,0,3,0,0,0.,0.,0,N_vam(rvec,normalUncLowerBnds)},
		{"means",14,0,1,1,0,0.,0.,0,N_vam(rvec,normalUncMeans)},
		{"nuv_descriptors",7,0,6,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_normal)},
		{"nuv_lower_bounds",6,0,3,0,0,0.,0.,-3,N_vam(rvec,normalUncLowerBnds)},
		{"nuv_means",6,0,1,1,0,0.,0.,-3,N_vam(rvec,normalUncMeans)},
		{"nuv_std_deviations",6,0,2,2,0,0.,0.,2,N_vam(RealLb,normalUncStdDevs)},
		{"nuv_upper_bounds",6,0,4,0,0,0.,0.,2,N_vam(rvec,normalUncUpperBnds)},
		{"std_deviations",14,0,2,2,0,0.,0.,0,N_vam(RealLb,normalUncStdDevs)},
		{"upper_bounds",14,0,4,0,0,0.,0.,0,N_vam(rvec,normalUncUpperBnds)}
		},
	kw_1120[3] = {
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_poisson)},
		{"initial_point",13,0,2,0,0,0.,0.,0,N_vam(IntLb,poissonUncVars)},
		{"lambdas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,poissonUncLambdas)}
		},
	kw_1121[9] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(caulbl,CAUVar_triangular)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,triangularUncVars)},
		{"lower_bounds",14,0,2,2,0,0.,0.,0,N_vam(RealLb,triangularUncLowerBnds)},
		{"modes",14,0,1,1,0,0.,0.,0,N_vam(rvec,triangularUncModes)},
		{"tuv_descriptors",7,0,5,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_triangular)},
		{"tuv_lower_bounds",6,0,2,2,0,0.,0.,-3,N_vam(RealLb,triangularUncLowerBnds)},
		{"tuv_modes",6,0,1,1,0,0.,0.,-3,N_vam(rvec,triangularUncModes)},
		{"tuv_upper_bounds",6,0,3,3,0,0.,0.,1,N_vam(RealUb,triangularUncUpperBnds)},
		{"upper_bounds",14,0,3,3,0,0.,0.,0,N_vam(RealUb,triangularUncUpperBnds)}
		},
	kw_1122[7] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_uniform)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,uniformUncVars)},
		{"lower_bounds",14,0,1,1,0,0.,0.,0,N_vam(RealLb,uniformUncLowerBnds)},
		{"upper_bounds",14,0,2,2,0,0.,0.,0,N_vam(RealUb,uniformUncUpperBnds)},
		{"uuv_descriptors",7,0,4,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_uniform)},
		{"uuv_lower_bounds",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,uniformUncLowerBnds)},
		{"uuv_upper_bounds",6,0,2,2,0,0.,0.,-3,N_vam(RealUb,uniformUncUpperBnds)}
		},
	kw_1123[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,weibullUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(RealLb,weibullUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_weibull)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,weibullUncVars)},
		{"wuv_alphas",6,0,1,1,0,0.,0.,-4,N_vam(RealLb,weibullUncAlphas)},
		{"wuv_betas",6,0,2,2,0,0.,0.,-4,N_vam(RealLb,weibullUncBetas)},
		{"wuv_descriptors",7,0,4,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_weibull)}
		},
	kw_1124[42] = {
		{"active",8,6,2,0,kw_1080},
		{"beta_uncertain",0x19,11,13,0,kw_1081,0.,0.,0,N_vam(sizet,numBetaUncVars)},
		{"binomial_uncertain",0x19,5,20,0,kw_1082,0.,0.,0,N_vam(sizet,numBinomialUncVars)},
		{"continuous_design",0x19,12,4,0,kw_1083,0.,0.,0,N_vam(sizet,numContinuousDesVars)},
		{"continuous_interval_uncertain",0x19,10,26,0,kw_1084,0.,0.,0,N_vam(sizet,numContinuousIntervalUncVars)},
		{"continuous_state",0x19,8,29,0,kw_1085,0.,0.,0,N_vam(sizet,numContinuousStateVars)},
		{"discrete_design_range",0x19,8,5,0,kw_1086,0.,0.,0,N_vam(sizet,numDiscreteDesRangeVars)},
		{"discrete_design_set",8,3,6,0,kw_1092},
		{"discrete_interval_uncertain",0x19,9,27,0,kw_1093,0.,0.,0,N_vam(sizet,numDiscreteIntervalUncVars)},
		{"discrete_state_range",0x19,8,30,0,kw_1094,0.,0.,0,N_vam(sizet,numDiscreteStateRangeVars)},
		{"discrete_state_set",8,3,31,0,kw_1098},
		{"discrete_uncertain_set",8,3,28,0,kw_1102},
		{"exponential_uncertain",0x19,5,12,0,kw_1103,0.,0.,0,N_vam(sizet,numExponentialUncVars)},
		{"frechet_uncertain",0x19,7,16,0,kw_1104,0.,0.,0,N_vam(sizet,numFrechetUncVars)},
		{"gamma_uncertain",0x19,7,14,0,kw_1105,0.,0.,0,N_vam(sizet,numGammaUncVars)},
		{"geometric_uncertain",0x19,4,22,0,kw_1106,0.,0.,0,N_vam(sizet,numGeometricUncVars)},
		{"gumbel_uncertain",0x19,7,15,0,kw_1107,0.,0.,0,N_vam(sizet,numGumbelUncVars)},
		{"histogram_bin_uncertain",0x19,11,18,0,kw_1108,0.,0.,0,N_vam(sizet,numHistogramBinUncVars)},
		{"histogram_point_uncertain",8,3,24,0,kw_1112},
		{"hypergeometric_uncertain",0x19,5,23,0,kw_1113,0.,0.,0,N_vam(sizet,numHyperGeomUncVars)},
		{"id_variables",11,0,1,0,0,0.,0.,0,N_vam(str,idVariables)},
		{"interval_uncertain",0x11,10,26,0,kw_1084,0.,0.,-17,N_vam(sizet,numContinuousIntervalUncVars)},
		{"linear_equality_constraint_matrix",14,0,37,0,0,0.,0.,0,N_vam(rvec,linearEqConstraintCoeffs)},
		{"linear_equality_scale_types",15,0,39,0,0,0.,0.,0,N_vam(strL,linearEqScaleTypes)},
		{"linear_equality_scales",14,0,40,0,0,0.,0.,0,N_vam(rvec,linearEqScales)},
		{"linear_equality_targets",14,0,38,0,0,0.,0.,0,N_vam(rvec,linearEqTargets)},
		{"linear_inequality_constraint_matrix",14,0,32,0,0,0.,0.,0,N_vam(rvec,linearIneqConstraintCoeffs)},
		{"linear_inequality_lower_bounds",14,0,33,0,0,0.,0.,0,N_vam(rvec,linearIneqLowerBnds)},
		{"linear_inequality_scale_types",15,0,35,0,0,0.,0.,0,N_vam(strL,linearIneqScaleTypes)},
		{"linear_inequality_scales",14,0,36,0,0,0.,0.,0,N_vam(rvec,linearIneqScales)},
		{"linear_inequality_upper_bounds",14,0,34,0,0,0.,0.,0,N_vam(rvec,linearIneqUpperBnds)},
		{"lognormal_uncertain",0x19,11,8,0,kw_1116,0.,0.,0,N_vam(sizet,numLognormalUncVars)},
		{"loguniform_uncertain",0x19,7,10,0,kw_1117,0.,0.,0,N_vam(sizet,numLoguniformUncVars)},
		{"mixed",8,0,3,0,0,0.,0.,0,N_vam(type,varsDomain_MIXED_DOMAIN)},
		{"negative_binomial_uncertain",0x19,5,21,0,kw_1118,0.,0.,0,N_vam(sizet,numNegBinomialUncVars)},
		{"normal_uncertain",0x19,11,7,0,kw_1119,0.,0.,0,N_vam(sizet,numNormalUncVars)},
		{"poisson_uncertain",0x19,3,19,0,kw_1120,0.,0.,0,N_vam(sizet,numPoissonUncVars)},
		{"relaxed",8,0,3,0,0,0.,0.,0,N_vam(type,varsDomain_RELAXED_DOMAIN)},
		{"triangular_uncertain",0x19,9,11,0,kw_1121,0.,0.,0,N_vam(sizet,numTriangularUncVars)},
		{"uncertain_correlation_matrix",14,0,25,0,0,0.,0.,0,N_vam(newrvec,Var_Info_ucm)},
		{"uniform_uncertain",0x19,7,9,0,kw_1122,0.,0.,0,N_vam(sizet,numUniformUncVars)},
		{"weibull_uncertain",0x19,7,17,0,kw_1123,0.,0.,0,N_vam(sizet,numWeibullUncVars)}
		},
	kw_1125[6] = {
		{"environment",0x108,15,1,1,kw_15,0.,0.,0,NIDRProblemDescDB::env_start},
		{"interface",8,12,5,5,kw_32,0.,0.,0,N_ifm3(start,0,stop)},
		{"method",8,104,2,2,kw_958,0.,0.,0,N_mdm3(start,0,stop)},
		{"model",8,12,3,3,kw_1052,0.,0.,0,N_mom3(start,0,stop)},
		{"responses",8,19,6,6,kw_1079,0.,0.,0,N_rem3(start,0,stop)},
		{"variables",8,42,4,4,kw_1124,0.,0.,0,N_vam3(start,0,stop)}
		};

} // namespace Dakota

#ifdef __cplusplus
extern "C" {
#endif
KeyWord Dakota_Keyword_Top = {"KeywordTop",0,6,0,0,Dakota::kw_1125};
#ifdef __cplusplus
}
#endif
