
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
		{"processors_per_analysis",0x19,0,1,0,0,0.,0.,0,N_ifm(pint,procsPerAnalysis)}
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
	kw_22[10] = {
		{"analysis_components",15,0,4,0,0,0.,0.,0,N_ifm(str2D,analysisComponents)},
		{"direct",8,1,3,1,kw_16,0.,0.,0,N_ifm(type,interfaceType_TEST_INTERFACE)},
		{"fork",8,10,3,1,kw_18,0.,0.,0,N_ifm(type,interfaceType_FORK_INTERFACE)},
		{"grid",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_GRID_INTERFACE)},
		{"input_filter",11,0,1,0,0,0.,0.,0,N_ifm(str,inputFilter)},
		{"matlab",8,0,3,1,0,0.,0.,0,N_ifm(type,interfaceType_MATLAB_INTERFACE)},
		{"output_filter",11,0,2,0,0,0.,0.,0,N_ifm(str,outputFilter)},
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
		{"analysis_concurrency",0x19,0,3,0,0,0.,0.,0,N_ifm(pint,asynchLocalAnalysisConcurrency)},
		{"evaluation_concurrency",0x19,0,1,0,0,0.,0.,0,N_ifm(pint,asynchLocalEvalConcurrency)},
		{"local_evaluation_scheduling",8,2,2,0,kw_24}
		},
	kw_26[1] = {
		{"size",0x19,0,1,0,0,0.,0.,0,N_ifm(pint,asynchLocalEvalConcurrency)}
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
		{"analysis_drivers",15,10,2,0,kw_22,0.,0.,0,N_ifm(strL,analysisDrivers)},
		{"analysis_scheduling",8,2,11,0,kw_23},
		{"analysis_servers",0x19,0,10,0,0,0.,0.,0,N_ifm(pint,analysisServers)},
		{"asynchronous",8,3,6,0,kw_25,0.,0.,0,N_ifm(true, asynchFlag)},
		{"batch",8,1,6,0,kw_26,0.,0.,0,N_ifm(true,batchEvalFlag)},
		{"deactivate",8,4,5,0,kw_28},
		{"evaluation_scheduling",8,2,8,0,kw_30},
		{"evaluation_servers",0x19,0,7,0,0,0.,0.,0,N_ifm(pint,evalServers)},
		{"failure_capture",8,4,4,0,kw_31},
		{"id_interface",11,0,1,0,0,0.,0.,0,N_ifm(str,idInterface)},
		{"processors_per_evaluation",0x19,0,9,0,0,0.,0.,0,N_ifm(pint,procsPerEval)}
		},
	kw_33[4] = {
		{"constant_liar",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_cl)},
		{"distance_penalty",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_distance)},
		{"naive",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_naive)},
		{"topology",8,0,1,1,0,0.,0.,0,N_mdm(lit,batchSelectionType_topology)}
		},
	kw_34[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_35[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_36[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_35,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_37[3] = {
		{"distance",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_distance)},
		{"gradient",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_gradient)},
		{"predicted_variance",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessMetricType_predicted_variance)}
		},
	kw_38[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_39[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_40[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_39,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_41[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_42[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_43[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_42}
		},
	kw_44[2] = {
		{"compute",8,3,2,0,kw_43},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_45[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_46[19] = {
		{"batch_selection",8,4,5,0,kw_33},
		{"distribution",8,2,14,0,kw_34},
		{"export_approx_points_file",11,3,8,0,kw_36,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,8,0,kw_36,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"fitness_metric",8,3,4,0,kw_37},
		{"gen_reliability_levels",14,1,13,0,kw_38,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,7,0,kw_40,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,7,0,kw_40,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,9,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,16,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,12,0,kw_41,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,6,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"response_levels",14,2,11,0,kw_44,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,15,0,kw_45},
		{"samples",1,0,1,0,0,0.,0.,-8,N_mdm(int,numSamples)},
		{"samples_on_emulator",9,0,3,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_47[7] = {
		{"merit1",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit1)},
		{"merit1_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit1_smooth)},
		{"merit2",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2)},
		{"merit2_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2_smooth)},
		{"merit2_squared",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit2_squared)},
		{"merit_max",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit_max)},
		{"merit_max_smooth",8,0,1,1,0,0.,0.,0,N_mdm(lit,meritFunction_merit_max_smooth)}
		},
	kw_48[2] = {
		{"blocking",8,0,1,1,0,0.,0.,0,N_mdm(lit,evalSynchronize_blocking)},
		{"nonblocking",8,0,1,1,0,0.,0.,0,N_mdm(lit,evalSynchronize_nonblocking)}
		},
	kw_49[13] = {
		{"constraint_penalty",10,0,7,0,0,0.,0.,0,N_mdm(Real,constrPenalty)},
		{"constraint_tolerance",10,0,9,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"contraction_factor",10,0,2,0,0,0.,0.,0,N_mdm(Real,contractStepLength)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initStepLength)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"merit_function",8,7,6,0,kw_47},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"smoothing_factor",10,0,8,0,0,0.,0.,0,N_mdm(Real,smoothFactor)},
		{"solution_accuracy",2,0,4,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,4,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"synchronization",8,2,5,0,kw_48},
		{"variable_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_50[1] = {
		{"hyperprior_betas",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,hyperPriorBetas)}
		},
	kw_51[5] = {
		{"both",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_BOTH)},
		{"hyperprior_alphas",14,1,2,0,kw_50,0.,0.,0,N_mdm(RealDL,hyperPriorAlphas)},
		{"one",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_ONE)},
		{"per_experiment",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_PER_EXPER)},
		{"per_response",8,0,1,1,0,0.,0.,0,N_mdm(utype,calibrateErrorMode_CALIBRATE_PER_RESP)}
		},
	kw_52[1] = {
		{"confidence_intervals",8,0,1,0,0,0.,0.,0,N_mdm(true,chainDiagnosticsCI)}
		},
	kw_53[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_54[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_53,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_55[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_54,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_54,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_56[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_57[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_58[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_59[3] = {
		{"adapted",8,2,1,1,kw_58,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_60[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_61[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_62[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_63[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_64[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_65[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_66[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_60,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_60,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_61,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_62,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_63,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_63,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_62,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_64,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_65,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_65,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_67[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_68[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_69[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_70[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_71[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_72[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_73[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_67,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_67,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_68,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_69,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_70,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_70,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_69,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_71,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_72,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_72,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_74[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_75[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_76[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_75,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_77[7] = {
		{"basis_type",8,3,2,0,kw_59},
		{"collocation_points_sequence",13,18,3,1,kw_66,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_73,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_74,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_76,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_76,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_78[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_79[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_78,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_80[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_79,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_79,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_81[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_82[2] = {
		{"dimension_adaptive",8,3,1,1,kw_81},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_83[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_84[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_85[16] = {
		{"allocation_control",8,1,3,0,kw_56},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_57},
		{"expansion_order_sequence",13,7,5,1,kw_77,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_80,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_80,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_80,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_82,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_83,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_84,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_86[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_87[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_88[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_89[3] = {
		{"dimension_adaptive",8,2,1,1,kw_88},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_90[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_91[2] = {
		{"dimension_adaptive",8,2,1,1,kw_90},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_92[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_93[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_94[11] = {
		{"allocation_control",8,1,3,0,kw_86},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"discrepancy_emulation",8,2,4,0,kw_87},
		{"h_refinement",8,3,1,0,kw_89,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_91,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_92,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_93,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_95[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_96[3] = {
		{"estimator_variance",8,1,1,1,kw_95,0.,0.,0,N_mdm(type,mlmfAllocControl_ESTIMATOR_VARIANCE)},
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_RIP_SAMPLING)}
		},
	kw_97[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_98[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_99[3] = {
		{"adapted",8,2,1,1,kw_98,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_100[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_101[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_102[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_103[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_104[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_105[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_106[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_100,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_100,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_101,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_102,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_103,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_103,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_102,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_104,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_105,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_105,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_107[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_108[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_109[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_110[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_111[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_112[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_113[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_107,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_107,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_108,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_109,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_110,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_110,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_109,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_111,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_112,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_112,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_114[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_115[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_116[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_115,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_117[7] = {
		{"basis_type",8,3,2,0,kw_99},
		{"collocation_points_sequence",13,18,3,1,kw_106,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_113,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_114,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_116,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_116,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_118[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_119[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_118,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_120[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_119,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_119,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_121[14] = {
		{"allocation_control",8,3,2,0,kw_96},
		{"askey",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,3,0,kw_97},
		{"expansion_order_sequence",13,7,4,1,kw_117,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,7,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"initial_samples",5,0,1,0,0,0.,0.,5,N_mdm(szarray,pilotSamples)},
		{"least_interpolation",0,6,4,1,kw_120,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,6,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,4,1,kw_120,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,4,1,kw_120,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"pilot_samples",13,0,1,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"wiener",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_122[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_123[3] = {
		{"adapted",8,2,1,1,kw_122,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_124[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_125[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_126[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_127[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_128[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_129[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_130[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_124,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_124,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_125,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_126,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_127,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_127,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_126,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_128,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_129,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_129,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_131[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_132[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_133[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_134[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_135[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_136[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_137[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_131,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_131,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_132,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_133,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_134,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_134,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_133,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_135,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_136,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_136,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_138[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_139[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_140[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_139,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_141[8] = {
		{"basis_type",8,3,2,0,kw_123},
		{"collocation_points",9,18,3,1,kw_130,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,18,3,1,kw_137,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,3,3,1,kw_138,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_140,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_140,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_142[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_143[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_142,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_144[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_143,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_143,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_145[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_146[2] = {
		{"dimension_adaptive",8,3,1,1,kw_145},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_147[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_148[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_149[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_141,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_144,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_144,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_144,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_146,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_147,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_148,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_150[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_151[3] = {
		{"dimension_adaptive",8,2,1,1,kw_150},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_152[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_153[2] = {
		{"dimension_adaptive",8,2,1,1,kw_152},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_154[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_155[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_156[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_151,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_153,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_154,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_155,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_157[7] = {
		{"gaussian_process",8,6,1,1,kw_55},
		{"kriging",0,6,1,1,kw_55,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_85,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,11,1,1,kw_94,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,14,1,1,kw_121,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_149,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_156,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
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
	kw_160[11] = {
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"chains",0x29,0,3,0,0,3.,0.,0,N_mdm(int,numChains)},
		{"crossover_chain_pairs",0x29,0,5,0,0,0.,0.,0,N_mdm(int,crossoverChainPairs)},
		{"emulator",8,7,8,0,kw_157},
		{"export_chain_points_file",11,3,10,0,kw_159,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"gr_threshold",0x1a,0,6,0,0,0.,0.,0,N_mdm(Real,grThreshold)},
		{"jump_step",0x29,0,7,0,0,0.,0.,0,N_mdm(int,jumpStep)},
		{"num_cr",0x29,0,4,0,0,1.,0.,0,N_mdm(int,numCR)},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"standardized_space",8,0,9,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_161[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importCandFormat_TABULAR_IFACE_ID)}
		},
	kw_162[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_161,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importCandFormat_TABULAR_NONE)}
		},
	kw_163[7] = {
		{"batch_size",0x29,0,4,0,0,1.,0.,0,N_mdm(int,batchSize)},
		{"import_candidate_points_file",11,3,5,0,kw_162,0.,0.,0,N_mdm(str,importCandPtsFile)},
		{"initial_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,numSamples)},
		{"ksg2",8,0,6,0,0,0.,0.,0,N_mdm(true,mutualInfoKSG2)},
		{"max_hifi_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(int,maxHifiEvals)},
		{"num_candidates",0x19,0,2,2,0,0.,0.,0,N_mdm(sizet,numCandidates)},
		{"samples",1,0,1,1,0,0.,0.,-4,N_mdm(int,numSamples)}
		},
	kw_164[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_165[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_164,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_166[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_167[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_166,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_168[1] = {
		{"update_period",9,0,1,0,0,0.,0.,0,N_mdm(int,proposalCovUpdatePeriod)}
		},
	kw_169[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_170[1] = {
		{"multiplier",0x1a,0,1,0,0,0.,0.,0,N_mdm(Real,priorPropCovMult)}
		},
	kw_171[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_172[4] = {
		{"derivatives",8,1,1,1,kw_168,0.,0.,0,N_mdm(lit,proposalCovType_derivatives)},
		{"filename",11,2,1,1,kw_169,0.,0.,0,N_mdm(str,proposalCovFile)},
		{"prior",8,1,1,1,kw_170,0.,0.,0,N_mdm(lit,proposalCovType_prior)},
		{"values",14,2,1,1,kw_171,0.,0.,0,N_mdm(RealDL,proposalCovData)}
		},
	kw_173[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_174[17] = {
		{"adaptive_metropolis",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_adaptive_metropolis)},
		{"build_samples",9,0,4,2,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"delayed_rejection",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_delayed_rejection)},
		{"dram",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_dram)},
		{"export_chain_points_file",11,3,9,0,kw_165,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"gpmsa_normalize",8,0,8,0,0,0.,0.,0,N_mdm(true,gpmsaNormalize)},
		{"import_build_points_file",11,3,5,0,kw_167,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,3,5,0,kw_167,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"logit_transform",8,0,7,0,0,0.,0.,0,N_mdm(true,logitTransform)},
		{"metropolis_hastings",8,0,10,0,0,0.,0.,0,N_mdm(lit,mcmcType_metropolis_hastings)},
		{"options_file",11,0,12,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"proposal_covariance",8,4,11,0,kw_172,0.,0.,0,N_mdm(lit,proposalCovType_user)},
		{"rng",8,2,3,0,kw_173},
		{"samples",1,0,1,1,0,0.,0.,-12,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"standardized_space",8,0,6,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_175[3] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mdm(order,approxCorrectionOrder_0)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mdm(order,approxCorrectionOrder_1)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mdm(order,approxCorrectionOrder_2)}
		},
	kw_176[4] = {
		{"correction_order",8,3,2,0,kw_175},
		{"gaussian_process",8,0,1,1,0,0.,0.,0,N_mdm(lit,modelDiscrepancyType_global_kriging)},
		{"kriging",0,0,1,1,0,0.,0.,-1,N_mdm(lit,modelDiscrepancyType_global_kriging)},
		{"polynomial",8,0,1,1,0,0.,0.,0,N_mdm(lit,modelDiscrepancyType_global_polynomial)}
		},
	kw_177[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrModelFormat_TABULAR_IFACE_ID)}
		},
	kw_178[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_177,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrModelFormat_TABULAR_NONE)}
		},
	kw_179[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportCorrVarFormat_TABULAR_IFACE_ID)}
		},
	kw_180[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_179,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportCorrVarFormat_TABULAR_NONE)}
		},
	kw_181[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportDiscrepFormat_TABULAR_IFACE_ID)}
		},
	kw_182[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_181,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportDiscrepFormat_TABULAR_NONE)}
		},
	kw_183[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importPredConfigFormat_TABULAR_IFACE_ID)}
		},
	kw_184[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_183,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importPredConfigFormat_TABULAR_NONE)}
		},
	kw_185[7] = {
		{"discrepancy_type",8,4,1,0,kw_176},
		{"export_corrected_model_file",11,3,6,0,kw_178,0.,0.,0,N_mdm(str,exportCorrModelFile)},
		{"export_corrected_variance_file",11,3,7,0,kw_180,0.,0.,0,N_mdm(str,exportCorrVarFile)},
		{"export_discrepancy_file",11,3,5,0,kw_182,0.,0.,0,N_mdm(str,exportDiscrepFile)},
		{"import_prediction_configs",11,3,4,0,kw_184,0.,0.,0,N_mdm(str,importPredConfigs)},
		{"num_prediction_configs",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numPredConfigs)},
		{"prediction_configs",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,predictionConfigList)}
		},
	kw_186[3] = {
		{"evidence_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,evidenceSamples)},
		{"laplace_approx",8,0,3,0,0,0.,0.,0,N_mdm(true,modelEvidLaplace)},
		{"mc_approx",8,0,1,0,0,0.,0.,0,N_mdm(true,modelEvidMC)}
		},
	kw_187[1] = {
		{"ksg2",8,0,1,0,0,0.,0.,0,N_mdm(true,mutualInfoKSG2)}
		},
	kw_188[3] = {
		{"kde",8,0,3,0,0,0.,0.,0,N_mdm(true,posteriorStatsKDE)},
		{"kl_divergence",8,0,1,0,0,0.,0.,0,N_mdm(true,posteriorStatsKL)},
		{"mutual_info",8,1,2,0,kw_187,0.,0.,0,N_mdm(true,posteriorStatsMutual)}
		},
	kw_189[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_190[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_191[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_190,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_192[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_191,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_191,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_193[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_194[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_195[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_196[3] = {
		{"adapted",8,2,1,1,kw_195,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_197[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_198[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_199[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_200[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_201[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_202[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_203[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_197,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_197,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_198,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_199,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_200,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_200,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_199,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_201,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_202,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_202,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_204[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_205[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_206[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_207[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_208[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_209[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_210[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_204,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_204,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_205,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_206,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_207,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_207,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_206,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_208,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_209,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_209,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_211[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_212[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_213[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_212,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_214[7] = {
		{"basis_type",8,3,2,0,kw_196},
		{"collocation_points_sequence",13,18,3,1,kw_203,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_210,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_211,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_213,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_213,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_215[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_216[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_215,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_217[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_216,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_216,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_218[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_219[2] = {
		{"dimension_adaptive",8,3,1,1,kw_218},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_220[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_221[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_222[16] = {
		{"allocation_control",8,1,3,0,kw_193},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_194},
		{"expansion_order_sequence",13,7,5,1,kw_214,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_217,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_217,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_217,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_219,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_220,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_221,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_223[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_224[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_225[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_226[3] = {
		{"dimension_adaptive",8,2,1,1,kw_225},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_227[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_228[2] = {
		{"dimension_adaptive",8,2,1,1,kw_227},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_229[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_230[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_231[11] = {
		{"allocation_control",8,1,3,0,kw_223},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"discrepancy_emulation",8,2,4,0,kw_224},
		{"h_refinement",8,3,1,0,kw_226,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_228,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_229,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_230,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_232[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_233[3] = {
		{"estimator_variance",8,1,1,1,kw_232,0.,0.,0,N_mdm(type,mlmfAllocControl_ESTIMATOR_VARIANCE)},
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_RIP_SAMPLING)}
		},
	kw_234[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_235[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_236[3] = {
		{"adapted",8,2,1,1,kw_235,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_237[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_238[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_239[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_240[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_241[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_242[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_243[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_237,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_237,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_238,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_239,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_240,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_240,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_239,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_241,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_242,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_242,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_244[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_245[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_246[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_247[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_248[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_249[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_250[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_244,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_244,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_245,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_246,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_247,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_247,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_246,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_248,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_249,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_249,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_251[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_252[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_253[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_252,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_254[7] = {
		{"basis_type",8,3,2,0,kw_236},
		{"collocation_points_sequence",13,18,3,1,kw_243,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_250,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_251,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_253,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_253,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_255[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_256[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_255,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_257[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_256,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_256,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_258[14] = {
		{"allocation_control",8,3,2,0,kw_233},
		{"askey",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,3,0,kw_234},
		{"expansion_order_sequence",13,7,4,1,kw_254,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,7,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"initial_samples",5,0,1,0,0,0.,0.,5,N_mdm(szarray,pilotSamples)},
		{"least_interpolation",0,6,4,1,kw_257,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,6,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,4,1,kw_257,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,4,1,kw_257,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"pilot_samples",13,0,1,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"wiener",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_259[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_260[3] = {
		{"adapted",8,2,1,1,kw_259,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_261[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_262[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_263[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_264[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_265[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_266[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_267[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_261,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_261,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_262,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_263,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_264,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_264,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_263,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_265,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_266,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_266,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_268[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_269[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_270[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_271[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_272[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_273[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_274[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_268,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_268,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_269,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_270,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_271,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_271,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_270,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_272,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_273,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_273,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_275[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_276[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_277[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_276,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_278[8] = {
		{"basis_type",8,3,2,0,kw_260},
		{"collocation_points",9,18,3,1,kw_267,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,18,3,1,kw_274,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,3,3,1,kw_275,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_277,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_277,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_279[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_280[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_279,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_281[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_280,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_280,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_282[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_283[2] = {
		{"dimension_adaptive",8,3,1,1,kw_282},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_284[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_285[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_286[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_278,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_281,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_281,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_281,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_283,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_284,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_285,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_287[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_288[3] = {
		{"dimension_adaptive",8,2,1,1,kw_287},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_289[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_290[2] = {
		{"dimension_adaptive",8,2,1,1,kw_289},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_291[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_292[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_293[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_288,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_290,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_291,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_292,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_294[7] = {
		{"gaussian_process",8,6,1,1,kw_192},
		{"kriging",0,6,1,1,kw_192,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_222,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,11,1,1,kw_231,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,14,1,1,kw_258,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_286,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_293,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
		},
	kw_295[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_296[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_295,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_297[3] = {
		{"nip",8,0,1,1,0,0.,0.,0,N_mdm(utype,preSolveMethod_SUBMETHOD_NIP)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(utype,preSolveMethod_SUBMETHOD_NONE)},
		{"sqp",8,0,1,1,0,0.,0.,0,N_mdm(utype,preSolveMethod_SUBMETHOD_SQP)}
		},
	kw_298[1] = {
		{"update_period",9,0,1,0,0,0.,0.,0,N_mdm(int,proposalCovUpdatePeriod)}
		},
	kw_299[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_300[1] = {
		{"multiplier",0x1a,0,1,0,0,0.,0.,0,N_mdm(Real,priorPropCovMult)}
		},
	kw_301[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,proposalCovInputType_matrix)}
		},
	kw_302[4] = {
		{"derivatives",8,1,1,1,kw_298,0.,0.,0,N_mdm(lit,proposalCovType_derivatives)},
		{"filename",11,2,1,1,kw_299,0.,0.,0,N_mdm(str,proposalCovFile)},
		{"prior",8,1,1,1,kw_300,0.,0.,0,N_mdm(lit,proposalCovType_prior)},
		{"values",14,2,1,1,kw_301,0.,0.,0,N_mdm(RealDL,proposalCovData)}
		},
	kw_303[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_304[16] = {
		{"adaptive_metropolis",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_adaptive_metropolis)},
		{"chain_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,chainSamples)},
		{"delayed_rejection",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_delayed_rejection)},
		{"dram",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_dram)},
		{"emulator",8,7,4,0,kw_294},
		{"export_chain_points_file",11,3,7,0,kw_296,0.,0.,0,N_mdm(str,exportMCMCPtsFile)},
		{"logit_transform",8,0,6,0,0,0.,0.,0,N_mdm(true,logitTransform)},
		{"metropolis_hastings",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_metropolis_hastings)},
		{"multilevel",8,0,8,0,0,0.,0.,0,N_mdm(lit,mcmcType_multilevel)},
		{"options_file",11,0,11,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"pre_solve",8,3,9,0,kw_297},
		{"proposal_covariance",8,4,10,0,kw_302,0.,0.,0,N_mdm(lit,proposalCovType_user)},
		{"rng",8,2,3,0,kw_303},
		{"samples",1,0,1,1,0,0.,0.,-12,N_mdm(int,chainSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"standardized_space",8,0,5,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_305[2] = {
		{"diagonal",8,0,1,1,0,0.,0.,0,N_mdm(lit,dataDistCovInputType_diagonal)},
		{"matrix",8,0,1,1,0,0.,0.,0,N_mdm(lit,dataDistCovInputType_matrix)}
		},
	kw_306[2] = {
		{"covariance",14,2,2,2,kw_305,0.,0.,0,N_mdm(RealDL,dataDistCovariance)},
		{"means",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,dataDistMeans)}
		},
	kw_307[2] = {
		{"gaussian",8,2,1,1,kw_306},
		{"obs_data_filename",11,0,1,1,0,0.,0.,0,N_mdm(str,dataDistFile)}
		},
	kw_308[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_309[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_308,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_310[6] = {
		{"build_samples",9,0,2,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"import_build_points_file",11,4,4,0,kw_309,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_309,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,3,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_311[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_312[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_313[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_314[3] = {
		{"adapted",8,2,1,1,kw_313,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_315[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_316[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_317[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_318[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_319[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_320[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_321[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_315,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_315,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_316,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_317,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_318,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_318,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_317,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_319,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_320,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_320,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_322[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_323[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_324[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_325[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_326[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_327[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_328[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_322,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_322,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_323,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_324,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_325,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_325,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_324,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_326,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_327,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_327,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_329[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_330[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_331[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_330,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_332[7] = {
		{"basis_type",8,3,2,0,kw_314},
		{"collocation_points_sequence",13,18,3,1,kw_321,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_328,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_329,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_331,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_331,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_333[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_334[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_333,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_335[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_334,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_334,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_336[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_337[2] = {
		{"dimension_adaptive",8,3,1,1,kw_336},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_338[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_339[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_340[16] = {
		{"allocation_control",8,1,3,0,kw_311},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_312},
		{"expansion_order_sequence",13,7,5,1,kw_332,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,9,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,6,5,1,kw_335,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_335,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_335,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_337,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order_sequence",13,3,5,1,kw_338,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_339,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_341[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_342[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_343[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_344[3] = {
		{"dimension_adaptive",8,2,1,1,kw_343},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_345[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_346[2] = {
		{"dimension_adaptive",8,2,1,1,kw_345},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_347[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_348[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_349[11] = {
		{"allocation_control",8,1,3,0,kw_341},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"discrepancy_emulation",8,2,4,0,kw_342},
		{"h_refinement",8,3,1,0,kw_344,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_346,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order_sequence",13,3,5,1,kw_347,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_348,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_350[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_351[3] = {
		{"estimator_variance",8,1,1,1,kw_350,0.,0.,0,N_mdm(type,mlmfAllocControl_ESTIMATOR_VARIANCE)},
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_RIP_SAMPLING)}
		},
	kw_352[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_353[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_354[3] = {
		{"adapted",8,2,1,1,kw_353,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_355[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_356[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_357[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_358[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_359[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_360[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_361[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_355,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_355,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_356,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_357,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_358,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_358,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_357,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_359,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_360,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_360,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_362[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_363[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_364[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_365[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_366[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_367[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_368[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_362,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_362,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_363,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_364,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_365,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_365,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_364,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_366,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_367,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_367,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_369[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_370[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_371[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_370,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_372[7] = {
		{"basis_type",8,3,2,0,kw_354},
		{"collocation_points_sequence",13,18,3,1,kw_361,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_368,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_369,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_371,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_371,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_373[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_374[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_373,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_375[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_374,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_374,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_376[14] = {
		{"allocation_control",8,3,2,0,kw_351},
		{"askey",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,3,0,kw_352},
		{"expansion_order_sequence",13,7,4,1,kw_372,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_expansion_file",11,0,7,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,8,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"initial_samples",5,0,1,0,0,0.,0.,5,N_mdm(szarray,pilotSamples)},
		{"least_interpolation",0,6,4,1,kw_375,0.,0.,3,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"normalized",8,0,6,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,4,1,kw_375,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,4,1,kw_375,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"pilot_samples",13,0,1,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"wiener",8,0,5,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_377[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_378[3] = {
		{"adapted",8,2,1,1,kw_377,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_379[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_380[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_381[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_382[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_383[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_384[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_385[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_379,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_379,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_380,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_381,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_382,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_382,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_381,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_383,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_384,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_384,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_386[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_387[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_388[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_389[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_390[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_391[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_392[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_386,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_386,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_387,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_388,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_389,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_389,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_388,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_390,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_391,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_391,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_393[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_394[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_395[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_394,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_396[8] = {
		{"basis_type",8,3,2,0,kw_378},
		{"collocation_points",9,18,3,1,kw_385,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,18,3,1,kw_392,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,3,3,1,kw_393,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_395,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_395,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)}
		},
	kw_397[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_398[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_397,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_399[7] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_398,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_398,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"posterior_adaptive",8,0,5,0,0,0.,0.,0,N_mdm(true,adaptPosteriorRefine)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_400[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_401[2] = {
		{"dimension_adaptive",8,3,1,1,kw_400},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_402[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_403[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_404[15] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"expansion_order",9,8,3,1,kw_396,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"full_covariance",8,0,7,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"least_interpolation",0,7,3,1,kw_399,0.,0.,4,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,7,3,1,kw_399,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,7,3,1,kw_399,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_401,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"quadrature_order",9,3,3,1,kw_402,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,5,3,1,kw_403,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_405[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_406[3] = {
		{"dimension_adaptive",8,2,1,1,kw_405},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_407[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_408[2] = {
		{"dimension_adaptive",8,2,1,1,kw_407},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_409[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_410[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_411[11] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"diagonal_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"full_covariance",8,0,6,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"h_refinement",8,3,1,0,kw_406,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"p_refinement",8,2,1,0,kw_408,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"quadrature_order",9,3,3,1,kw_409,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"sparse_grid_level",9,7,3,1,kw_410,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_412[7] = {
		{"gaussian_process",8,6,1,1,kw_310},
		{"kriging",0,6,1,1,kw_310,0.,0.,-1},
		{"mf_pce",8,16,1,1,kw_340,0.,0.,0,N_mdm(type,emulatorType_MF_PCE_EMULATOR)},
		{"mf_sc",8,11,1,1,kw_349,0.,0.,0,N_mdm(type,emulatorType_MF_SC_EMULATOR)},
		{"ml_pce",8,14,1,1,kw_376,0.,0.,0,N_mdm(type,emulatorType_ML_PCE_EMULATOR)},
		{"pce",8,15,1,1,kw_404,0.,0.,0,N_mdm(type,emulatorType_PCE_EMULATOR)},
		{"sc",8,11,1,1,kw_411,0.,0.,0,N_mdm(type,emulatorType_SC_EMULATOR)}
		},
	kw_413[1] = {
		{"posterior_density_export_filename",11,0,1,0,0,0.,0.,0,N_mdm(str,posteriorDensityExportFilename)}
		},
	kw_414[1] = {
		{"posterior_samples_export_filename",11,0,1,0,0,0.,0.,0,N_mdm(str,posteriorSamplesExportFilename)}
		},
	kw_415[8] = {
		{"data_distribution",8,2,5,2,kw_307},
		{"emulator",8,7,3,0,kw_412},
		{"evaluate_posterior_density",8,1,8,0,kw_413,0.,0.,0,N_mdm(true,evaluatePosteriorDensity)},
		{"generate_posterior_samples",8,1,7,0,kw_414,0.,0.,0,N_mdm(true,generatePosteriorSamples)},
		{"posterior_samples_import_filename",11,0,6,0,0,0.,0.,0,N_mdm(str,posteriorSamplesImportFilename)},
		{"pushforward_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,numPushforwardSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"standardized_space",8,0,4,0,0,0.,0.,0,N_mdm(true,standardizedSpace)}
		},
	kw_416[17] = {
		{"burn_in_samples",9,0,4,0,0,0.,0.,0,N_mdm(int,burnInSamples)},
		{"calibrate_error_multipliers",8,5,3,0,kw_51},
		{"chain_diagnostics",8,1,6,0,kw_52,0.,0.,0,N_mdm(true,chainDiagnostics)},
		{"convergence_tolerance",10,0,11,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"dream",8,11,1,1,kw_160,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_DREAM)},
		{"experimental_design",8,7,2,0,kw_163,0.,0.,0,N_mdm(true,adaptExpDesign)},
		{"gpmsa",8,17,1,1,kw_174,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_GPMSA)},
		{"max_iterations",0x29,0,12,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_discrepancy",8,7,8,0,kw_185,0.,0.,0,N_mdm(true,calModelDiscrepancy)},
		{"model_evidence",8,3,7,0,kw_186,0.,0.,0,N_mdm(true,modelEvidence)},
		{"model_pointer",11,0,13,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"posterior_stats",8,3,5,0,kw_188},
		{"probability_levels",14,1,10,0,kw_189,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"queso",8,16,1,1,kw_304,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_QUESO)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"sub_sampling_period",9,0,9,0,0,0.,0.,0,N_mdm(int,subSamplingPeriod)},
		{"wasabi",8,8,1,1,kw_415,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_WASABI)}
		},
	kw_417[1] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,modelPointer)}
		},
	kw_418[3] = {
		{"method_name",11,1,1,1,kw_417,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"scaling",8,0,2,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_419[4] = {
		{"deltas_per_variable",5,0,2,2,0,0.,0.,3,N_mdm(ivec,stepsPerVariable)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"step_vector",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,stepVector)},
		{"steps_per_variable",13,0,2,2,0,0.,0.,0,N_mdm(ivec,stepsPerVariable)}
		},
	kw_420[11] = {
		{"beta_solver_name",11,0,1,1,0,0.,0.,0,N_mdm(str,betaSolverName)},
		{"convergence_tolerance",10,0,7,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,5,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,4,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,2,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,2,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_421[12] = {
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,9,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,6,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,10,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,5,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,3,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,3,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_422[2] = {
		{"all_dimensions",8,0,1,1,0,0.,0.,0,N_mdm(lit,boxDivision_all_dimensions)},
		{"major_dimension",8,0,1,1,0,0.,0.,0,N_mdm(lit,boxDivision_major_dimension)}
		},
	kw_423[16] = {
		{"constraint_penalty",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"convergence_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"division",8,2,1,0,kw_422},
		{"global_balance_parameter",10,0,2,0,0,0.,0.,0,N_mdm(Real,globalBalanceParam)},
		{"local_balance_parameter",10,0,3,0,0,0.,0.,0,N_mdm(Real,localBalanceParam)},
		{"max_boxsize_limit",10,0,4,0,0,0.,0.,0,N_mdm(Real,maxBoxSize)},
		{"max_function_evaluations",0x29,0,13,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"min_boxsize_limit",10,0,5,0,0,0.,0.,0,N_mdm(Real,minBoxSize)},
		{"misc_options",15,0,10,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,8,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,9,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,7,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,7,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_424[3] = {
		{"blend",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_blend)},
		{"two_point",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_two_point)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,crossoverType_uniform)}
		},
	kw_425[2] = {
		{"linear_rank",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_linear_rank)},
		{"merit_function",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_proportional)}
		},
	kw_426[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_427[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_428[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_429[2] = {
		{"mutation_range",9,0,2,0,0,0.,0.,0,N_mdm(int,mutationRange)},
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real,mutationScale)}
		},
	kw_430[5] = {
		{"non_adaptive",8,0,2,0,0,0.,0.,0,N_mdm(false,mutationAdaptive)},
		{"offset_cauchy",8,2,1,1,kw_427,0.,0.,0,N_mdm(lit,mutationType_offset_cauchy)},
		{"offset_normal",8,2,1,1,kw_428,0.,0.,0,N_mdm(lit,mutationType_offset_normal)},
		{"offset_uniform",8,2,1,1,kw_429,0.,0.,0,N_mdm(lit,mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_431[4] = {
		{"chc",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_chc)},
		{"elitist",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_elitist)},
		{"new_solutions_generated",9,0,2,0,0,0.,0.,0,N_mdm(int,newSolnsGenerated)},
		{"random",9,0,1,1,0,0.,0.,0,N_mdm(ilit2,TYPE_DATA_replacementType_random)}
		},
	kw_432[19] = {
		{"constraint_penalty",10,0,9,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"convergence_tolerance",10,0,15,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"crossover_rate",10,0,5,0,0,0.,0.,0,N_mdm(Real,crossoverRate)},
		{"crossover_type",8,3,6,0,kw_424},
		{"fitness_type",8,2,3,0,kw_425},
		{"initialization_type",8,3,2,0,kw_426},
		{"max_function_evaluations",0x29,0,16,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,14,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,13,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,18,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_rate",10,0,7,0,0,0.,0.,0,N_mdm(Real,mutationRate)},
		{"mutation_type",8,5,8,0,kw_430},
		{"population_size",0x19,0,1,0,0,0.,0.,0,N_mdm(pint,populationSize)},
		{"replacement_type",8,4,4,0,kw_431},
		{"scaling",8,0,17,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,11,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,12,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,10,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,10,0,0,0.,0.,0,N_mdm(Real,solnTarget)}
		},
	kw_433[3] = {
		{"adaptive_pattern",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_adaptive)},
		{"basic_pattern",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_simple)},
		{"multi_step",8,0,1,1,0,0.,0.,0,N_mdm(lit,exploratoryMoves_multi_step)}
		},
	kw_434[2] = {
		{"coordinate",8,0,1,1,0,0.,0.,0,N_mdm(lit,patternBasis_coordinate)},
		{"simplex",8,0,1,1,0,0.,0.,0,N_mdm(lit,patternBasis_simplex)}
		},
	kw_435[2] = {
		{"blocking",8,0,1,1,0,0.,0.,0,N_mdm(lit,evalSynchronize_blocking)},
		{"nonblocking",8,0,1,1,0,0.,0.,0,N_mdm(lit,evalSynchronize_nonblocking)}
		},
	kw_436[22] = {
		{"constant_penalty",8,0,1,0,0,0.,0.,0,N_mdm(true,constantPenalty)},
		{"constraint_penalty",10,0,10,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"contraction_factor",10,0,9,0,0,0.,0.,0,N_mdm(Real,contractFactor)},
		{"convergence_tolerance",10,0,18,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"expand_after_success",9,0,3,0,0,0.,0.,0,N_mdm(int,expandAfterSuccess)},
		{"exploratory_moves",8,3,7,0,kw_433},
		{"initial_delta",10,0,11,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,19,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,17,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,16,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,21,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"no_expansion",8,0,2,0,0,0.,0.,0,N_mdm(false,expansionFlag)},
		{"pattern_basis",8,2,4,0,kw_434},
		{"scaling",8,0,20,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,14,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,15,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,13,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,13,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"stochastic",8,0,5,0,0,0.,0.,0,N_mdm(true,randomizeOrderFlag)},
		{"synchronization",8,2,8,0,kw_435},
		{"total_pattern_size",9,0,6,0,0,0.,0.,0,N_mdm(int,totalPatternSize)},
		{"variable_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_437[18] = {
		{"constant_penalty",8,0,4,0,0,0.,0.,0,N_mdm(true,constantPenalty)},
		{"constraint_penalty",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintPenalty)},
		{"contract_after_failure",9,0,1,0,0,0.,0.,0,N_mdm(int,contractAfterFail)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,contractFactor)},
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"expand_after_success",9,0,3,0,0,0.,0.,0,N_mdm(int,expandAfterSuccess)},
		{"initial_delta",10,0,7,0,0,0.,0.,0,N_mdm(Real,initDelta)},
		{"max_function_evaluations",0x29,0,15,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,13,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,12,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,17,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"no_expansion",8,0,2,0,0,0.,0.,0,N_mdm(false,expansionFlag)},
		{"scaling",8,0,16,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,10,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"show_misc_options",8,0,11,0,0,0.,0.,0,N_mdm(true,showMiscOptions)},
		{"solution_accuracy",2,0,9,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,9,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,threshDelta)}
		},
	kw_438[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_439[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_440[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_441[15] = {
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
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"symbols",9,0,8,0,0,0.,0.,0,N_mdm(int,numSymbols)},
		{"variance_based_decomp",8,1,7,0,kw_440,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_442[7] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"options_file",11,0,6,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"solution_accuracy",2,0,5,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,5,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"variable_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_443[3] = {
		{"max_function_evaluations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,2,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_444[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_445[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_446[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_447[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_448[7] = {
		{"constraint_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,4,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_449[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_450[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_449,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_451[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_452[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_453[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_452,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_454[13] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"export_approx_points_file",11,3,9,0,kw_450,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,9,0,kw_450,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,6,0,kw_451},
		{"import_build_points_file",11,4,8,0,kw_453,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,8,0,kw_453,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"kriging",0,2,6,0,kw_451,0.,0.,-4},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"x_conv_tol",10,0,5,0,0,0.,0.,0,N_mdm(Real,xConvTol)}
		},
	kw_455[3] = {
		{"grid",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_grid)},
		{"halton",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_halton)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(lit,trialType_random)}
		},
	kw_456[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_457[10] = {
		{"fixed_seed",8,0,3,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"latinize",8,0,4,0,0,0.,0.,0,N_mdm(true,latinizeFlag)},
		{"max_iterations",0x29,0,9,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"num_trials",9,0,8,0,0,0.,0.,0,N_mdm(int,numTrials)},
		{"quality_metrics",8,0,5,0,0,0.,0.,0,N_mdm(true,volQualityFlag)},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"trial_type",8,3,7,0,kw_455},
		{"variance_based_decomp",8,1,6,0,kw_456,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_458[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_459[12] = {
		{"fixed_sequence",8,0,6,0,0,0.,0.,0,N_mdm(true,fixedSequenceFlag)},
		{"halton",8,0,1,1,0,0.,0.,0,N_mdm(utype,methodName_FSU_HALTON)},
		{"hammersley",8,0,1,1,0,0.,0.,0,N_mdm(utype,methodName_FSU_HAMMERSLEY)},
		{"latinize",8,0,2,0,0,0.,0.,0,N_mdm(true,latinizeFlag)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"prime_base",13,0,9,0,0,0.,0.,0,N_mdm(ivec,primeBase)},
		{"quality_metrics",8,0,3,0,0,0.,0.,0,N_mdm(true,volQualityFlag)},
		{"samples",9,0,5,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sequence_leap",13,0,8,0,0,0.,0.,0,N_mdm(ivec,sequenceLeap)},
		{"sequence_start",13,0,7,0,0,0.,0.,0,N_mdm(ivec,sequenceStart)},
		{"variance_based_decomp",8,1,4,0,kw_458,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_460[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_461[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_462[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_461,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_463[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_464[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_465[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_466[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_465,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_467[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_468[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_469[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_470[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_471[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_470}
		},
	kw_472[2] = {
		{"compute",8,4,2,0,kw_471},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_473[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_474[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_475[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_476[22] = {
		{"convergence_tolerance",10,0,15,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,14,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,12,0,kw_460},
		{"export_approx_points_file",11,3,17,0,kw_462,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,17,0,kw_462,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,7,0,kw_463},
		{"fixed_seed",8,0,4,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,14,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,11,0,kw_464,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,16,0,kw_466,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"model_pointer",11,0,18,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,9,0,kw_467,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,6,0,kw_468},
		{"reliability_levels",14,1,10,0,kw_469,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,8,0,kw_472,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,5,0,kw_473},
		{"sample_refinement",0,4,6,0,kw_468,0.,0.,-4},
		{"sample_type",8,2,2,0,kw_474},
		{"samples",1,0,1,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,1,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"variance_based_decomp",8,2,13,0,kw_475,0.,0.,0,N_mdm(true,vbdFlag)}
		},
	kw_477[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_478[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_479[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_478,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_480[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_481[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_482[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_481,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_483[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_484[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_485[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_484}
		},
	kw_486[2] = {
		{"compute",8,3,2,0,kw_485},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_487[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_488[15] = {
		{"build_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,10,0,kw_477},
		{"export_approx_points_file",11,3,5,0,kw_479,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,5,0,kw_479,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gen_reliability_levels",14,1,9,0,kw_480,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,4,0,kw_482,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_482,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,8,0,kw_483,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,7,0,kw_486,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_487},
		{"samples",1,0,1,0,0,0.,0.,-12,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,3,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_489[4] = {
		{"max_function_evaluations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,3,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_490[4] = {
		{"max_function_evaluations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,3,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_491[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
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
	kw_498[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_499[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_500[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_501[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_500}
		},
	kw_502[2] = {
		{"compute",8,3,2,0,kw_501},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_503[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_504[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_505[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_504,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_506[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_507[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_508[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_507,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_509[7] = {
		{"export_approx_points_file",11,3,4,0,kw_505,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_505,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,1,0,kw_506},
		{"import_build_points_file",11,4,3,0,kw_508,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_508,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,2,1,0,kw_506,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_510[12] = {
		{"distribution",8,2,7,0,kw_491},
		{"ea",8,0,3,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EA)},
		{"ego",8,7,3,0,kw_497,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGO)},
		{"gen_reliability_levels",14,1,6,0,kw_498,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lhs",8,0,3,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_LHS)},
		{"model_pointer",11,0,9,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,5,0,kw_499,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,4,0,kw_502,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,8,0,kw_503},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sbo",8,7,3,0,kw_509,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SBO)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_511[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_512[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_511,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_513[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_514[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_515[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_514,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_516[7] = {
		{"export_approx_points_file",11,3,4,0,kw_512,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_512,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,1,0,kw_513},
		{"import_build_points_file",11,4,3,0,kw_515,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_515,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,2,1,0,kw_513,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_517[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_518[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_519[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_518,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_520[2] = {
		{"dakota",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"surfpack",8,0,1,1,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)}
		},
	kw_521[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_522[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_521,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_523[7] = {
		{"export_approx_points_file",11,3,4,0,kw_519,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,4,0,kw_519,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gaussian_process",8,2,1,0,kw_520},
		{"import_build_points_file",11,4,3,0,kw_522,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,3,0,kw_522,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"kriging",0,2,1,0,kw_520,0.,0.,-3},
		{"use_derivatives",8,0,2,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_524[11] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"ea",8,0,6,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EA)},
		{"ego",8,7,6,0,kw_516,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EGO)},
		{"lhs",8,0,6,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_LHS)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"rng",8,2,7,0,kw_517},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"sbo",8,7,6,0,kw_523,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SBO)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_525[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_526[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_527[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_526,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_528[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_529[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_530[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_529,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_531[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_532[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_533[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_532}
		},
	kw_534[2] = {
		{"compute",8,3,2,0,kw_533},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_535[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_536[21] = {
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"dakota",8,0,3,0,0,0.,0.,0,N_mdm(type,emulatorType_GP_EMULATOR)},
		{"distribution",8,2,12,0,kw_525},
		{"export_approx_points_file",11,3,5,0,kw_527,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,5,0,kw_527,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"gen_reliability_levels",14,1,11,0,kw_528,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_build_points_file",11,4,4,0,kw_530,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_530,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"initial_samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,13,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,10,0,kw_531,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,9,0,kw_534,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,8,0,kw_535},
		{"seed",0x19,0,7,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"surfpack",8,0,3,0,0,0.,0.,0,N_mdm(type,emulatorType_KRIGING_EMULATOR)},
		{"u_gaussian_process",8,0,2,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_EGRA_U)},
		{"u_kriging",0,0,2,1,0,0.,0.,-1,N_mdm(utype,reliabilitySearchType_EGRA_U)},
		{"use_derivatives",8,0,6,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"x_gaussian_process",8,0,2,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_EGRA_X)},
		{"x_kriging",0,0,2,1,0,0.,0.,-1,N_mdm(utype,reliabilitySearchType_EGRA_X)}
		},
	kw_537[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_538[1] = {
		{"model_pointer_list",11,0,1,0,0,0.,0.,0,N_mdm(strL,hybridModelPointers)}
		},
	kw_539[5] = {
		{"iterator_scheduling",8,2,3,0,kw_537},
		{"iterator_servers",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,iteratorServers)},
		{"method_name_list",15,1,1,1,kw_538,0.,0.,0,N_mdm(strL,hybridMethodNames)},
		{"method_pointer_list",15,0,1,1,0,0.,0.,0,N_mdm(strL,hybridMethodPointers)},
		{"processors_per_iterator",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,procsPerIterator)}
		},
	kw_540[1] = {
		{"global_model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,hybridGlobalModelPointer)}
		},
	kw_541[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_542[1] = {
		{"local_model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,hybridLocalModelPointer)}
		},
	kw_543[8] = {
		{"global_method_name",11,1,1,1,kw_540,0.,0.,0,N_mdm(str,hybridGlobalMethodName)},
		{"global_method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,hybridGlobalMethodPointer)},
		{"iterator_scheduling",8,2,5,0,kw_541},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,iteratorServers)},
		{"local_method_name",11,1,2,2,kw_542,0.,0.,0,N_mdm(str,hybridLocalMethodName)},
		{"local_method_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,hybridLocalMethodPointer)},
		{"local_search_probability",10,0,3,0,0,0.,0.,0,N_mdm(Real,hybridLSProb)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(pint,procsPerIterator)}
		},
	kw_544[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_545[1] = {
		{"model_pointer_list",11,0,1,0,0,0.,0.,0,N_mdm(strL,hybridModelPointers)}
		},
	kw_546[5] = {
		{"iterator_scheduling",8,2,3,0,kw_544},
		{"iterator_servers",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,iteratorServers)},
		{"method_name_list",15,1,1,1,kw_545,0.,0.,0,N_mdm(strL,hybridMethodNames)},
		{"method_pointer_list",15,0,1,1,0,0.,0.,0,N_mdm(strL,hybridMethodPointers)},
		{"processors_per_iterator",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,procsPerIterator)}
		},
	kw_547[5] = {
		{"collaborative",8,5,1,1,kw_539,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_COLLABORATIVE)},
		{"coupled",0,8,1,1,kw_543,0.,0.,1,N_mdm(utype,subMethod_SUBMETHOD_EMBEDDED)},
		{"embedded",8,8,1,1,kw_543,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_EMBEDDED)},
		{"sequential",8,5,1,1,kw_546,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SEQUENTIAL)},
		{"uncoupled",0,5,1,1,kw_546,0.,0.,-1,N_mdm(utype,subMethod_SUBMETHOD_SEQUENTIAL)}
		},
	kw_548[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_549[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_550[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_551[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_552[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_551}
		},
	kw_553[2] = {
		{"compute",8,3,2,0,kw_552},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_554[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_555[15] = {
		{"adapt_import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"convergence_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,10,0,kw_548},
		{"gen_reliability_levels",14,1,9,0,kw_549,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"initial_samples",1,0,1,0,0,0.,0.,8,N_mdm(int,numSamples)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"mm_adapt_import",8,0,3,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,8,0,kw_550,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,4,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"response_levels",14,2,7,0,kw_553,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_554},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_556[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,pstudyFileFormat_TABULAR_IFACE_ID)}
		},
	kw_557[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,pstudyFileActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_556,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,pstudyFileFormat_TABULAR_NONE)}
		},
	kw_558[3] = {
		{"import_points_file",11,4,1,1,kw_557,0.,0.,0,N_mdm(str,pstudyFilename)},
		{"list_of_points",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,listOfPoints)},
		{"model_pointer",11,0,2,0,0,0.,0.,0,N_mdm(str,modelPointer)}
		},
	kw_559[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_560[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_561[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_562[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_563[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_562}
		},
	kw_564[2] = {
		{"compute",8,3,2,0,kw_563},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_565[7] = {
		{"distribution",8,2,5,0,kw_559},
		{"gen_reliability_levels",14,1,4,0,kw_560,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"nip",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NIP)},
		{"probability_levels",14,1,3,0,kw_561,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,2,0,kw_564,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"sqp",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SQP)}
		},
	kw_566[4] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"nip",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NIP)},
		{"sqp",8,0,1,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SQP)}
		},
	kw_567[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_568[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_569[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_570[5] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_571[4] = {
		{"first_order",8,0,1,1,0,0.,0.,0,N_mdm(lit,reliabilityIntegration_first_order)},
		{"probability_refinement",8,5,2,0,kw_570},
		{"sample_refinement",0,5,2,0,kw_570,0.,0.,-1},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mdm(lit,reliabilityIntegration_second_order)}
		},
	kw_572[12] = {
		{"integration",8,4,3,0,kw_571},
		{"nip",8,0,2,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_NIP)},
		{"no_approx",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_NO_APPROX)},
		{"sqp",8,0,2,0,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_SQP)},
		{"u_multi_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_QMEA_U)},
		{"u_taylor_mean",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_AMV_U)},
		{"u_taylor_mpp",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_AMV_PLUS_U)},
		{"u_two_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_TANA_U)},
		{"x_multi_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_QMEA_X)},
		{"x_taylor_mean",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_AMV_X)},
		{"x_taylor_mpp",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_AMV_PLUS_X)},
		{"x_two_point",8,0,1,1,0,0.,0.,0,N_mdm(utype,reliabilitySearchType_TANA_X)}
		},
	kw_573[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_574[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_575[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_576[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_575}
		},
	kw_577[2] = {
		{"compute",8,4,2,0,kw_576},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_578[10] = {
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,6,0,kw_567},
		{"final_moments",8,3,9,0,kw_568},
		{"gen_reliability_levels",14,1,5,0,kw_569,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mpp_search",8,12,1,0,kw_572},
		{"probability_levels",14,1,3,0,kw_573,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"reliability_levels",14,1,4,0,kw_574,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,2,0,kw_577,0.,0.,0,N_mdm(resplevs,responseLevels)}
		},
	kw_579[2] = {
		{"inform_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,useSurrogate_inform_search)},
		{"optimize",8,0,1,1,0,0.,0.,0,N_mdm(lit,useSurrogate_optimize)}
		},
	kw_580[14] = {
		{"display_all_evaluations",8,0,9,0,0,0.,0.,0,N_mdm(true,showAllEval)},
		{"display_format",11,0,6,0,0,0.,0.,0,N_mdm(str,displayFormat)},
		{"function_precision",10,0,3,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"history_file",11,0,5,0,0,0.,0.,0,N_mdm(str,historyFile)},
		{"initial_delta",10,0,1,0,0,0.,0.,0,N_mdm(Real,initMeshSize)},
		{"max_function_evaluations",0x29,0,12,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,14,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"neighbor_order",0x19,0,8,0,0,0.,0.,0,N_mdm(pint,neighborOrder)},
		{"scaling",8,0,13,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"use_surrogate",8,2,10,0,kw_579},
		{"variable_neighborhood_search",10,0,7,0,0,0.,0.,0,N_mdm(Real,vns)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,minMeshSize)}
		},
	kw_581[3] = {
		{"metric_tracker",8,0,1,1,0,0.,0.,0,N_mdm(lit,convergenceType_metric_tracker)},
		{"num_generations",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,2,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_582[2] = {
		{"num_offspring",0x19,0,2,0,0,0.,0.,0,N_mdm(pintz,numOffspring)},
		{"num_parents",0x19,0,1,0,0,0.,0.,0,N_mdm(pintz,numParents)}
		},
	kw_583[5] = {
		{"crossover_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_crossoverType_null_crossover)},
		{"multi_point_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_binary)},
		{"multi_point_parameterized_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_parameterized_binary)},
		{"multi_point_real",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_real)},
		{"shuffle_random",8,2,1,1,kw_582,0.,0.,0,N_mdm(litc,TYPE_DATA_crossoverType_shuffle_random)}
		},
	kw_584[2] = {
		{"domination_count",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_domination_count)},
		{"layer_rank",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_layer_rank)}
		},
	kw_585[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_586[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_587[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_588[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_589[6] = {
		{"bit_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_bit_random)},
		{"mutation_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_mutationType_null_mutation)},
		{"offset_cauchy",8,1,1,1,kw_586,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_cauchy)},
		{"offset_normal",8,1,1,1,kw_587,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_normal)},
		{"offset_uniform",8,1,1,1,kw_588,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_590[1] = {
		{"num_designs",0x29,0,1,0,0,2.,0.,0,N_mdm(pintz,numDesigns)}
		},
	kw_591[3] = {
		{"distance",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_distance)},
		{"max_designs",14,1,1,1,kw_590,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_max_designs)},
		{"radial",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_nichingType_radial)}
		},
	kw_592[1] = {
		{"orthogonal_distance",14,0,1,1,0,0.,0.,0,N_mdm(RealLlit,TYPE_DATA_postProcessorType_distance_postprocessor)}
		},
	kw_593[2] = {
		{"shrinkage_fraction",10,0,1,0,0,0.,0.,0,N_mdm(Real01,shrinkagePercent)},
		{"shrinkage_percentage",2,0,1,0,0,0.,0.,-1,N_mdm(Real01,shrinkagePercent)}
		},
	kw_594[4] = {
		{"below_limit",10,2,1,1,kw_593,0.,0.,0,N_mdm(litp,TYPE_DATA_replacementType_below_limit)},
		{"elitist",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_elitist)},
		{"roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_roulette_wheel)},
		{"unique_roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_unique_roulette_wheel)}
		},
	kw_595[17] = {
		{"convergence_tolerance",10,0,16,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"convergence_type",8,3,4,0,kw_581},
		{"crossover_type",8,5,13,0,kw_583},
		{"fitness_type",8,2,1,0,kw_584},
		{"initialization_type",8,3,12,0,kw_585},
		{"log_file",11,0,10,0,0,0.,0.,0,N_mdm(str,logFile)},
		{"max_function_evaluations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,17,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_type",8,6,14,0,kw_589},
		{"niching_type",8,3,3,0,kw_591},
		{"population_size",0x29,0,9,0,0,0.,0.,0,N_mdm(nnint,populationSize)},
		{"postprocessor_type",8,1,5,0,kw_592},
		{"print_each_pop",8,0,11,0,0,0.,0.,0,N_mdm(true,printPopFlag)},
		{"replacement_type",8,4,2,0,kw_594},
		{"scaling",8,0,8,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,15,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_596[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_597[1] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,subModelPointer)}
		},
	kw_598[1] = {
		{"seed",9,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_599[7] = {
		{"iterator_scheduling",8,2,5,0,kw_596},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,iteratorServers)},
		{"method_name",11,1,1,1,kw_597,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(pint,procsPerIterator)},
		{"random_starts",9,1,2,0,kw_598,0.,0.,0,N_mdm(int,concurrentRandomJobs)},
		{"starting_points",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,concurrentParameterSets)}
		},
	kw_600[2] = {
		{"model_pointer",11,0,2,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"partitions",13,0,1,1,0,0.,0.,0,N_mdm(usharray,varPartitions)}
		},
	kw_601[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_602[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_603[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_604[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_605[3] = {
		{"adapted",8,2,1,1,kw_604,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_606[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_607[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_608[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_609[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_610[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_611[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_612[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_606,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_606,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_607,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_608,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_609,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_609,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_608,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_610,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_611,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_611,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_613[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_614[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_615[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_616[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_617[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_618[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_619[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_613,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_613,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_614,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_615,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_616,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_616,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_615,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_617,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_618,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_618,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_620[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_621[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_622[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_621,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_623[7] = {
		{"basis_type",8,3,2,0,kw_605},
		{"collocation_points_sequence",13,18,3,1,kw_612,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_619,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_620,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_622,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_622,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_624[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_625[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_624,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_626[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_627[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_628[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_629[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_628,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_630[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_631[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_630,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_632[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_631,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_631,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_633[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_634[2] = {
		{"dimension_adaptive",8,3,1,1,kw_633},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_635[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_636[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_637[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_638[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_639[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_640[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_639}
		},
	kw_641[2] = {
		{"compute",8,4,2,0,kw_640},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_642[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_643[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_644[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_645[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_646[36] = {
		{"allocation_control",8,1,3,0,kw_601},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,23,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_602},
		{"distribution",8,2,20,0,kw_603},
		{"expansion_order_sequence",13,7,5,1,kw_623,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_approx_points_file",11,3,25,0,kw_625,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,25,0,kw_625,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,15,0,kw_626},
		{"fixed_seed",8,0,12,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,19,0,kw_627,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,24,0,kw_629,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"least_interpolation",0,6,5,1,kw_632,0.,0.,5,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"model_pointer",11,0,26,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_632,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_632,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_634,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"probability_levels",14,1,17,0,kw_635,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,14,0,kw_636},
		{"quadrature_order_sequence",13,3,5,1,kw_637,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"reliability_levels",14,1,18,0,kw_638,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,16,0,kw_641,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,13,0,kw_642},
		{"sample_refinement",0,4,14,0,kw_636,0.,0.,-5},
		{"sample_type",8,2,10,0,kw_643},
		{"samples",1,0,9,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,9,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,11,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"sparse_grid_level_sequence",13,5,5,1,kw_644,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"variance_based_decomp",8,2,21,0,kw_645,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_647[1] = {
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)}
		},
	kw_648[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_649[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_650[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_651[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_650,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_652[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_653[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_654[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_655[3] = {
		{"dimension_adaptive",8,2,1,1,kw_654},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_656[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_657[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_656,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_658[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_659[2] = {
		{"dimension_adaptive",8,2,1,1,kw_658},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_660[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_661[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_662[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_663[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_664[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_665[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_664}
		},
	kw_666[2] = {
		{"compute",8,4,2,0,kw_665},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_667[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_668[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_669[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_670[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_671[33] = {
		{"allocation_control",8,1,3,0,kw_647},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,22,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,21,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_648},
		{"distribution",8,2,19,0,kw_649},
		{"export_approx_points_file",11,3,24,0,kw_651,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,24,0,kw_651,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,14,0,kw_652},
		{"fixed_seed",8,0,11,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,21,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,18,0,kw_653,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"h_refinement",8,3,1,0,kw_655,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"import_approx_points_file",11,4,23,0,kw_657,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"model_pointer",11,0,25,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"p_refinement",8,2,1,0,kw_659,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,6,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"probability_levels",14,1,16,0,kw_660,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,13,0,kw_661},
		{"quadrature_order_sequence",13,3,5,1,kw_662,0.,0.,0,N_mdm(usharray,quadratureOrderSeq)},
		{"reliability_levels",14,1,17,0,kw_663,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,15,0,kw_666,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,12,0,kw_667},
		{"sample_refinement",0,4,13,0,kw_661,0.,0.,-5},
		{"sample_type",8,2,9,0,kw_668},
		{"samples",1,0,8,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,8,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,10,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"sparse_grid_level_sequence",13,7,5,1,kw_669,0.,0.,0,N_mdm(usharray,sparseGridLevelSeq)},
		{"use_derivatives",8,0,7,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"variance_based_decomp",8,2,20,0,kw_670,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_672[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_673[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportSamplesFormat_TABULAR_IFACE_ID)}
		},
	kw_674[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_673,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportSamplesFormat_TABULAR_NONE)}
		},
	kw_675[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_676[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_677[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_678[12] = {
		{"convergence_tolerance",10,0,7,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"distribution",8,2,9,0,kw_672},
		{"export_sample_sequence",8,3,5,0,kw_674,0.,0.,0,N_mdm(true,exportSampleSeqFlag)},
		{"final_moments",8,3,8,0,kw_675},
		{"fixed_seed",8,0,2,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"initial_samples",5,0,3,0,0,0.,0.,3,N_mdm(szarray,pilotSamples)},
		{"max_iterations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,11,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"pilot_samples",13,0,3,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"rng",8,2,10,0,kw_676},
		{"sample_type",8,2,4,0,kw_677},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_679[1] = {
		{"estimator_rate",10,0,1,0,0,0.,0.,0,N_mdm(Real,multilevEstimatorRate)}
		},
	kw_680[3] = {
		{"estimator_variance",8,1,1,1,kw_679,0.,0.,0,N_mdm(type,mlmfAllocControl_ESTIMATOR_VARIANCE)},
		{"greedy",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_GREEDY_REFINEMENT)},
		{"rip_sampling",8,0,1,1,0,0.,0.,0,N_mdm(type,mlmfAllocControl_RIP_SAMPLING)}
		},
	kw_681[2] = {
		{"distinct",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_DISTINCT_EMULATION)},
		{"recursive",8,0,1,1,0,0.,0.,0,N_mdm(type,multilevDiscrepEmulation_RECURSIVE_EMULATION)}
		},
	kw_682[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_683[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_684[3] = {
		{"adapted",8,2,1,1,kw_683,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_685[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_686[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_687[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_688[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_689[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_690[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_691[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_685,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_685,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_686,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_687,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_688,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_688,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_687,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_689,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_690,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_690,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_692[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_693[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_694[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_695[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_696[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_697[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_698[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_692,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_692,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_693,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_694,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_695,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_695,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_694,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_696,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_697,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_697,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_699[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_700[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_701[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_700,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_702[7] = {
		{"basis_type",8,3,2,0,kw_684},
		{"collocation_points_sequence",13,18,3,1,kw_691,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"collocation_ratio",10,18,3,1,kw_698,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples_sequence",13,3,3,1,kw_699,0.,0.,0,N_mdm(szarray,expansionSamplesSeq)},
		{"import_build_points_file",11,4,4,0,kw_701,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_701,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_703[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_704[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_703,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_705[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_706[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
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
	kw_709[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_710[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_709,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_711[6] = {
		{"collocation_points_sequence",13,0,1,1,0,0.,0.,0,N_mdm(szarray,collocationPointsSeq)},
		{"import_build_points_file",11,4,4,0,kw_710,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_710,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
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
	kw_714[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_715[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_716[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_715}
		},
	kw_717[2] = {
		{"compute",8,4,2,0,kw_716},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_718[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_719[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_720[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_721[35] = {
		{"allocation_control",8,3,3,0,kw_680},
		{"askey",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,23,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"discrepancy_emulation",8,2,4,0,kw_681},
		{"distribution",8,2,20,0,kw_682},
		{"expansion_order_sequence",13,7,5,1,kw_702,0.,0.,0,N_mdm(usharray,expansionOrderSeq)},
		{"export_approx_points_file",11,3,25,0,kw_704,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,8,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,25,0,kw_704,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,15,0,kw_705},
		{"fixed_seed",8,0,12,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,22,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,19,0,kw_706,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,24,0,kw_708,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"initial_samples",5,0,2,0,0,0.,0.,7,N_mdm(szarray,pilotSamples)},
		{"least_interpolation",0,6,5,1,kw_711,0.,0.,5,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,26,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,7,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,5,1,kw_711,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,5,1,kw_711,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"pilot_samples",13,0,2,0,0,0.,0.,0,N_mdm(szarray,pilotSamples)},
		{"probability_levels",14,1,17,0,kw_712,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,14,0,kw_713},
		{"reliability_levels",14,1,18,0,kw_714,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,16,0,kw_717,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,13,0,kw_718},
		{"sample_refinement",0,4,14,0,kw_713,0.,0.,-4},
		{"sample_type",8,2,10,0,kw_719},
		{"samples",1,0,9,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,9,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,11,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"variance_based_decomp",8,2,21,0,kw_720,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,6,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_722[9] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"min_boxsize_limit",10,0,2,0,0,0.,0.,0,N_mdm(Real,minBoxSize)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,7,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"solution_accuracy",2,0,1,0,0,0.,0.,1,N_mdm(Real,solnTarget)},
		{"solution_target",10,0,1,0,0,0.,0.,0,N_mdm(Real,solnTarget)},
		{"volume_boxsize_limit",10,0,3,0,0,0.,0.,0,N_mdm(Real,volBoxSize)}
		},
	kw_723[15] = {
		{"absolute_conv_tol",10,0,2,0,0,0.,0.,0,N_mdm(Real,absConvTol)},
		{"convergence_tolerance",10,0,10,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"covariance",9,0,8,0,0,0.,0.,0,N_mdm(int,covarianceType)},
		{"false_conv_tol",10,0,6,0,0,0.,0.,0,N_mdm(Real,falseConvTol)},
		{"function_precision",10,0,1,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"initial_trust_radius",10,0,7,0,0,0.,0.,0,N_mdm(Real,initTRRadius)},
		{"max_function_evaluations",0x29,0,13,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,11,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"regression_diagnostics",8,0,9,0,0,0.,0.,0,N_mdm(true,regressDiag)},
		{"scaling",8,0,14,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"singular_conv_tol",10,0,4,0,0,0.,0.,0,N_mdm(Real,singConvTol)},
		{"singular_radius",10,0,5,0,0,0.,0.,0,N_mdm(Real,singRadius)},
		{"speculative",8,0,12,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"x_conv_tol",10,0,3,0,0,0.,0.,0,N_mdm(Real,xConvTol)}
		},
	kw_724[5] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_725[10] = {
		{"constraint_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"function_precision",10,0,2,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"linesearch_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,lineSearchTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,7,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"verify_level",9,0,1,0,0,0.,0.,0,N_mdm(int,verifyLevel)}
		},
	kw_726[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_727[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_728[2] = {
		{"global",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_global)},
		{"local",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_local)}
		},
	kw_729[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_730[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_731[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_730}
		},
	kw_732[2] = {
		{"compute",8,3,2,0,kw_731},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_733[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_734[11] = {
		{"build_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,8,0,kw_726},
		{"gen_reliability_levels",14,1,7,0,kw_727,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lipschitz",8,2,3,0,kw_728},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,6,0,kw_729,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,5,0,kw_732,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,9,0,kw_733},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,4,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_735[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_736[2] = {
		{"advancements",9,0,1,0,0,0.,0.,0,N_mdm(ushint,adaptedBasisAdvancements)},
		{"soft_convergence_limit",9,0,2,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)}
		},
	kw_737[3] = {
		{"adapted",8,2,1,1,kw_736,0.,0.,0,N_mdm(type,expansionBasisType_ADAPTED_BASIS_EXPANDING_FRONT)},
		{"tensor_product",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TENSOR_PRODUCT_BASIS)},
		{"total_order",8,0,1,1,0,0.,0.,0,N_mdm(type,expansionBasisType_TOTAL_ORDER_BASIS)}
		},
	kw_738[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_739[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_740[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_741[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_742[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_743[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_744[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_738,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_738,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_739,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_740,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_741,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_741,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_740,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_742,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_743,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_743,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_745[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_746[1] = {
		{"noise_only",8,0,1,0,0,0.,0.,0,N_mdm(true,crossValidNoiseOnly)}
		},
	kw_747[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_748[2] = {
		{"l2_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Real,regressionL2Penalty)},
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_749[2] = {
		{"equality_constrained",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_EQ_CON_LS)},
		{"svd",8,0,1,0,0,0.,0.,0,N_mdm(type,lsRegressionType_SVD_LS)}
		},
	kw_750[1] = {
		{"noise_tolerance",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,regressionNoiseTol)}
		},
	kw_751[18] = {
		{"basis_pursuit",8,0,1,0,0,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"basis_pursuit_denoising",8,1,1,0,kw_745,0.,0.,0,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"bp",0,0,1,0,0,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT)},
		{"bpdn",0,1,1,0,kw_745,0.,0.,-2,N_mdm(type,regressionType_BASIS_PURSUIT_DENOISING)},
		{"cross_validation",8,1,2,0,kw_746,0.,0.,0,N_mdm(true,crossValidation)},
		{"lars",0,1,1,0,kw_747,0.,0.,3,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"lasso",0,2,1,0,kw_748,0.,0.,1,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_absolute_shrinkage",8,2,1,0,kw_748,0.,0.,0,N_mdm(type,regressionType_LASSO_REGRESSION)},
		{"least_angle_regression",8,1,1,0,kw_747,0.,0.,0,N_mdm(type,regressionType_LEAST_ANGLE_REGRESSION)},
		{"least_squares",8,2,1,0,kw_749,0.,0.,0,N_mdm(type,regressionType_DEFAULT_LEAST_SQ_REGRESSION)},
		{"max_solver_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"omp",0,1,1,0,kw_750,0.,0.,1,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"orthogonal_matching_pursuit",8,1,1,0,kw_750,0.,0.,0,N_mdm(type,regressionType_ORTHOG_MATCH_PURSUIT)},
		{"ratio_order",10,0,3,0,0,0.,0.,0,N_mdm(Realp,collocRatioTermsOrder)},
		{"reuse_points",8,0,6,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,6,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",8,0,5,0,0,0.,0.,0,N_mdm(true,tensorGridFlag)},
		{"use_derivatives",8,0,4,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)}
		},
	kw_752[3] = {
		{"incremental_lhs",8,0,2,0,0,0.,0.,0,N_mdm(lit,expansionSampleType_incremental_lhs)},
		{"reuse_points",8,0,1,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,1,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)}
		},
	kw_753[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_754[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_753,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_755[7] = {
		{"basis_type",8,3,2,0,kw_737},
		{"collocation_points",9,18,3,1,kw_744,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"collocation_ratio",10,18,3,1,kw_751,0.,0.,0,N_mdm(Realp,collocationRatio)},
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"expansion_samples",9,3,3,1,kw_752,0.,0.,0,N_mdm(sizet,expansionSamples)},
		{"import_build_points_file",11,4,4,0,kw_754,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_754,0.,0.,-1,N_mdm(str,importBuildPtsFile)}
		},
	kw_756[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_757[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_756,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_758[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_759[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_760[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_761[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_760,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_762[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_763[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importBuildActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_762,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_764[6] = {
		{"collocation_points",9,0,1,1,0,0.,0.,0,N_mdm(sizet,collocationPoints)},
		{"import_build_points_file",11,4,4,0,kw_763,0.,0.,0,N_mdm(str,importBuildPtsFile)},
		{"import_points_file",3,4,4,0,kw_763,0.,0.,-1,N_mdm(str,importBuildPtsFile)},
		{"reuse_points",8,0,3,0,0,0.,0.,0,N_mdm(lit,pointReuse_all)},
		{"reuse_samples",0,0,3,0,0,0.,0.,-1,N_mdm(lit,pointReuse_all)},
		{"tensor_grid",13,0,2,0,0,0.,0.,0,N_mdm(usharray,tensorGridOrder)}
		},
	kw_765[3] = {
		{"decay",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_DECAY)},
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_766[2] = {
		{"dimension_adaptive",8,3,1,1,kw_765},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_767[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_768[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_769[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_770[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_771[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_772[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_771}
		},
	kw_773[2] = {
		{"compute",8,4,2,0,kw_772},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_774[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_775[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_776[5] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,3,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,2,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_777[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_778[36] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,21,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"cubature_integrand",9,0,3,1,0,0.,0.,0,N_mdm(ushint,cubIntOrder)},
		{"diagonal_covariance",8,0,20,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,18,0,kw_735},
		{"expansion_order",9,7,3,1,kw_755,0.,0.,0,N_mdm(ushint,expansionOrder)},
		{"export_approx_points_file",11,3,23,0,kw_757,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_expansion_file",11,0,6,0,0,0.,0.,0,N_mdm(str,exportExpansionFile)},
		{"export_points_file",3,3,23,0,kw_757,0.,0.,-2,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,13,0,kw_758},
		{"fixed_seed",8,0,10,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,20,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,17,0,kw_759,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"import_approx_points_file",11,4,22,0,kw_761,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"import_expansion_file",11,0,3,1,0,0.,0.,0,N_mdm(str,importExpansionFile)},
		{"least_interpolation",0,6,3,1,kw_764,0.,0.,5,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"model_pointer",11,0,24,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"normalized",8,0,5,0,0,0.,0.,0,N_mdm(true,normalizedCoeffs)},
		{"oli",0,6,3,1,kw_764,0.,0.,1,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"orthogonal_least_interpolation",8,6,3,1,kw_764,0.,0.,0,N_mdm(type,regressionType_ORTHOG_LEAST_INTERPOLATION)},
		{"p_refinement",8,2,1,0,kw_766,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"probability_levels",14,1,15,0,kw_767,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,12,0,kw_768},
		{"quadrature_order",9,3,3,1,kw_769,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"reliability_levels",14,1,16,0,kw_770,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,14,0,kw_773,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,11,0,kw_774},
		{"sample_refinement",0,4,12,0,kw_768,0.,0.,-5},
		{"sample_type",8,2,8,0,kw_775},
		{"samples",1,0,7,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,7,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,9,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"sparse_grid_level",9,5,3,1,kw_776,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"variance_based_decomp",8,2,19,0,kw_777,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_779[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_780[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_781[2] = {
		{"global",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_global)},
		{"local",8,0,1,1,0,0.,0.,0,N_mdm(lit,lipschitzType_local)}
		},
	kw_782[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_783[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_784[3] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"system",8,2,2,0,kw_783}
		},
	kw_785[2] = {
		{"compute",8,3,2,0,kw_784},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_786[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_787[11] = {
		{"build_samples",9,0,1,1,0,0.,0.,0,N_mdm(int,buildSamples)},
		{"distribution",8,2,8,0,kw_779},
		{"gen_reliability_levels",14,1,7,0,kw_780,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"lipschitz",8,2,3,0,kw_781},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"probability_levels",14,1,6,0,kw_782,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"response_levels",14,2,5,0,kw_785,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,9,0,kw_786},
		{"samples",1,0,1,1,0,0.,0.,-8,N_mdm(int,buildSamples)},
		{"samples_on_emulator",9,0,4,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_788[2] = {
		{"candidate_designs",0x19,0,1,0,0,0.,0.,0,N_mdm(sizet,numCandidateDesigns)},
		{"leja_oversample_ratio",10,0,1,0,0,0.,0.,0,N_mdm(Real,collocationRatio)}
		},
	kw_789[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_790[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_791[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_792[1] = {
		{"percent_variance_explained",10,0,1,0,0,0.,0.,0,N_mdm(Real,percentVarianceExplained)}
		},
	kw_793[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_794[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_795[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_796[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_795}
		},
	kw_797[2] = {
		{"compute",8,4,2,0,kw_796},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_798[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_799[4] = {
		{"incremental_lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"incremental_random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)},
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_800[1] = {
		{"drop_tolerance",10,0,1,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)}
		},
	kw_801[5] = {
		{"confidence_level",10,0,2,0,0,0.,0.,0,N_mdm(Real,wilksConfidenceLevel)},
		{"one_sided_lower",8,0,3,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_ONE_SIDED_LOWER)},
		{"one_sided_upper",8,0,4,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_ONE_SIDED_UPPER)},
		{"order",9,0,1,0,0,0.,0.,0,N_mdm(ushint,wilksOrder)},
		{"two_sided",8,0,5,0,0,0.,0.,0,N_mdm(type,wilksSidedInterval_TWO_SIDED)}
		},
	kw_802[19] = {
		{"backfill",8,0,8,0,0,0.,0.,0,N_mdm(true,backfillFlag)},
		{"d_optimal",8,2,6,0,kw_788,0.,0.,0,N_mdm(true,dOptimal)},
		{"distribution",8,2,16,0,kw_789},
		{"final_moments",8,3,11,0,kw_790},
		{"fixed_seed",8,0,3,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"gen_reliability_levels",14,1,15,0,kw_791,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"initial_samples",1,0,1,0,0,0.,0.,9,N_mdm(int,numSamples)},
		{"model_pointer",11,0,18,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"principal_components",8,1,9,0,kw_792,0.,0.,0,N_mdm(true,pcaFlag)},
		{"probability_levels",14,1,13,0,kw_793,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"refinement_samples",13,0,5,0,0,0.,0.,0,N_mdm(ivec,refineSamples)},
		{"reliability_levels",14,1,14,0,kw_794,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,12,0,kw_797,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,17,0,kw_798},
		{"sample_type",8,4,4,0,kw_799},
		{"samples",9,0,1,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,2,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"variance_based_decomp",8,1,7,0,kw_800,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wilks",8,5,10,0,kw_801,0.,0.,0,N_mdm(true,wilksFlag)}
		},
	kw_803[2] = {
		{"complementary",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_COMPLEMENTARY)},
		{"cumulative",8,0,1,1,0,0.,0.,0,N_mdm(type,distributionType_CUMULATIVE)}
		},
	kw_804[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_805[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_804,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_806[3] = {
		{"central",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_CENTRAL_MOMENTS)},
		{"none",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_NO_MOMENTS)},
		{"standard",8,0,1,1,0,0.,0.,0,N_mdm(type,finalMomentsType_STANDARD_MOMENTS)}
		},
	kw_807[1] = {
		{"num_gen_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,genReliabilityLevels)}
		},
	kw_808[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_809[3] = {
		{"dimension_adaptive",8,2,1,1,kw_808},
		{"local_adaptive",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_LOCAL_ADAPTIVE_CONTROL)},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_810[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mdm(augment_utype,importApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_811[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mdm(true,importApproxActive)},
		{"annotated",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_810,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mdm(utype,importApproxFormat_TABULAR_NONE)}
		},
	kw_812[2] = {
		{"generalized",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_GENERALIZED)},
		{"sobol",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_DIMENSION_ADAPTIVE_CONTROL_SOBOL)}
		},
	kw_813[2] = {
		{"dimension_adaptive",8,2,1,1,kw_812},
		{"uniform",8,0,1,1,0,0.,0.,0,N_mdm(type,refinementControl_UNIFORM_CONTROL)}
		},
	kw_814[1] = {
		{"num_probability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,probabilityLevels)}
		},
	kw_815[4] = {
		{"adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_AIS)},
		{"import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_IS)},
		{"mm_adapt_import",8,0,1,1,0,0.,0.,0,N_mdm(utype,integrationRefine_MMAIS)},
		{"refinement_samples",13,0,2,0,0,0.,0.,0,N_mdm(ivec,refineSamples)}
		},
	kw_816[3] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"non_nested",8,0,2,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)}
		},
	kw_817[1] = {
		{"num_reliability_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,reliabilityLevels)}
		},
	kw_818[2] = {
		{"parallel",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_PARALLEL)},
		{"series",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTargetReduce_SYSTEM_SERIES)}
		},
	kw_819[4] = {
		{"gen_reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_GEN_RELIABILITIES)},
		{"probabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_PROBABILITIES)},
		{"reliabilities",8,0,1,1,0,0.,0.,0,N_mdm(type,responseLevelTarget_RELIABILITIES)},
		{"system",8,2,2,0,kw_818}
		},
	kw_820[2] = {
		{"compute",8,4,2,0,kw_819},
		{"num_response_levels",13,0,1,0,0,0.,0.,0,N_mdm(num_resplevs,responseLevels)}
		},
	kw_821[2] = {
		{"mt19937",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_mt19937)},
		{"rnum2",8,0,1,1,0,0.,0.,0,N_mdm(lit,rngName_rnum2)}
		},
	kw_822[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mdm(utype,sampleType_SUBMETHOD_RANDOM)}
		},
	kw_823[7] = {
		{"dimension_preference",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,anisoDimPref)},
		{"hierarchical",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_HIERARCHICAL_INTERPOLANT)},
		{"nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NESTED)},
		{"nodal",8,0,2,0,0,0.,0.,0,N_mdm(type,expansionBasisType_NODAL_INTERPOLANT)},
		{"non_nested",8,0,4,0,0,0.,0.,0,N_mdm(type,nestingOverride_NON_NESTED)},
		{"restricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_RESTRICTED)},
		{"unrestricted",8,0,3,0,0,0.,0.,0,N_mdm(type,growthOverride_UNRESTRICTED)}
		},
	kw_824[2] = {
		{"drop_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,vbdDropTolerance)},
		{"interaction_order",0x19,0,1,0,0,0.,0.,0,N_mdm(ushint,vbdOrder)}
		},
	kw_825[31] = {
		{"askey",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_ASKEY_U)},
		{"convergence_tolerance",10,0,20,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"diagonal_covariance",8,0,19,0,0,0.,0.,0,N_mdm(type,covarianceControl_DIAGONAL_COVARIANCE)},
		{"distribution",8,2,17,0,kw_803},
		{"export_approx_points_file",11,3,22,0,kw_805,0.,0.,0,N_mdm(str,exportApproxPtsFile)},
		{"export_points_file",3,3,22,0,kw_805,0.,0.,-1,N_mdm(str,exportApproxPtsFile)},
		{"final_moments",8,3,12,0,kw_806},
		{"fixed_seed",8,0,9,0,0,0.,0.,0,N_mdm(true,fixedSeedFlag)},
		{"full_covariance",8,0,19,0,0,0.,0.,0,N_mdm(type,covarianceControl_FULL_COVARIANCE)},
		{"gen_reliability_levels",14,1,16,0,kw_807,0.,0.,0,N_mdm(resplevs,genReliabilityLevels)},
		{"h_refinement",8,3,1,0,kw_809,0.,0.,0,N_mdm(type,refinementType_H_REFINEMENT)},
		{"import_approx_points_file",11,4,21,0,kw_811,0.,0.,0,N_mdm(str,importApproxPtsFile)},
		{"max_refinement_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxRefineIterations)},
		{"model_pointer",11,0,23,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"p_refinement",8,2,1,0,kw_813,0.,0.,0,N_mdm(type,refinementType_P_REFINEMENT)},
		{"piecewise",8,0,4,0,0,0.,0.,0,NIDRProblemDescDB::method_piecewise},
		{"probability_levels",14,1,14,0,kw_814,0.,0.,0,N_mdm(resplevs01,probabilityLevels)},
		{"probability_refinement",8,4,11,0,kw_815},
		{"quadrature_order",9,3,3,1,kw_816,0.,0.,0,N_mdm(ushint,quadratureOrder)},
		{"reliability_levels",14,1,15,0,kw_817,0.,0.,0,N_mdm(resplevs,reliabilityLevels)},
		{"response_levels",14,2,13,0,kw_820,0.,0.,0,N_mdm(resplevs,responseLevels)},
		{"rng",8,2,10,0,kw_821},
		{"sample_refinement",0,4,11,0,kw_815,0.,0.,-5},
		{"sample_type",8,2,7,0,kw_822},
		{"samples",1,0,6,0,0,0.,0.,1,N_mdm(int,samplesOnEmulator)},
		{"samples_on_emulator",9,0,6,0,0,0.,0.,0,N_mdm(int,samplesOnEmulator)},
		{"seed",0x19,0,8,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"sparse_grid_level",9,7,3,1,kw_823,0.,0.,0,N_mdm(ushint,sparseGridLevel)},
		{"use_derivatives",8,0,5,0,0,0.,0.,0,N_mdm(true,methodUseDerivsFlag)},
		{"variance_based_decomp",8,2,18,0,kw_824,0.,0.,0,N_mdm(true,vbdFlag)},
		{"wiener",8,0,4,0,0,0.,0.,0,N_mdm(type,expansionType_STD_NORMAL_U)}
		},
	kw_826[5] = {
		{"convergence_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"misc_options",15,0,1,0,0,0.,0.,0,N_mdm(strL,miscOptions)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)}
		},
	kw_827[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_828[5] = {
		{"max_function_evaluations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,4,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"trust_region",8,6,1,0,kw_827,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final}
		},
	kw_829[10] = {
		{"constraint_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"function_precision",10,0,2,0,0,0.,0.,0,N_mdm(Real,functionPrecision)},
		{"linesearch_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,lineSearchTolerance)},
		{"max_function_evaluations",0x29,0,8,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,10,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,9,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,7,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"verify_level",9,0,1,0,0,0.,0.,0,N_mdm(int,verifyLevel)}
		},
	kw_830[8] = {
		{"convergence_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,6,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"max_step",10,0,1,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"model_pointer",11,0,8,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,7,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"speculative",8,0,5,0,0,0.,0.,0,N_mdm(true,speculativeFlag)}
		},
	kw_831[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_832[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_833[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_831},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_832},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_834[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_835[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_836[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_834},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_835},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_837[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_838[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_839[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_837},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_838},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_840[6] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_function_evaluations",0x29,0,4,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,2,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_scheme_size",9,0,1,0,0,0.,0.,0,N_mdm(int,searchSchemeSize)}
		},
	kw_841[3] = {
		{"argaez_tapia",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_ArgaezTapia)},
		{"el_bakry",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_NormFmu)},
		{"van_shanno",8,0,1,1,0,0.,0.,0,N_mdm(type,meritFn_VanShanno)}
		},
	kw_842[4] = {
		{"gradient_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_gradient_based_line_search)},
		{"tr_pds",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_tr_pds)},
		{"trust_region",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_trust_region)},
		{"value_based_line_search",8,0,1,1,0,0.,0.,0,N_mdm(lit,searchMethod_value_based_line_search)}
		},
	kw_843[12] = {
		{"centering_parameter",10,0,4,0,0,0.,0.,0,N_mdm(Real,centeringParam)},
		{"convergence_tolerance",10,0,8,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"gradient_tolerance",10,0,6,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_function_evaluations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"max_step",10,0,5,0,0,0.,0.,0,N_mdm(Real,maxStep)},
		{"merit_function",8,3,2,0,kw_841},
		{"model_pointer",11,0,12,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,11,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"search_method",8,4,1,0,kw_842},
		{"speculative",8,0,9,0,0,0.,0.,0,N_mdm(true,speculativeFlag)},
		{"steplength_to_boundary",10,0,3,0,0,0.,0.,0,N_mdm(Real,stepLenToBoundary)}
		},
	kw_844[5] = {
		{"debug",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_DEBUG_OUTPUT)},
		{"normal",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_NORMAL_OUTPUT)},
		{"quiet",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_QUIET_OUTPUT)},
		{"silent",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_SILENT_OUTPUT)},
		{"verbose",8,0,1,1,0,0.,0.,0,N_mdm(type,methodOutput_VERBOSE_OUTPUT)}
		},
	kw_845[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mdm(type,iteratorScheduling_PEER_SCHEDULING)}
		},
	kw_846[2] = {
		{"model_pointer",11,0,1,0,0,0.,0.,0,N_mdm(str,subModelPointer)},
		{"opt_model_pointer",3,0,1,0,0,0.,0.,-1,N_mdm(str,subModelPointer)}
		},
	kw_847[1] = {
		{"seed",9,0,1,0,0,0.,0.,0,N_mdm(int,randomSeed)}
		},
	kw_848[10] = {
		{"iterator_scheduling",8,2,5,0,kw_845},
		{"iterator_servers",0x19,0,4,0,0,0.,0.,0,N_mdm(pint,iteratorServers)},
		{"method_name",11,2,1,1,kw_846,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"multi_objective_weight_sets",6,0,3,0,0,0.,0.,5,N_mdm(RealDL,concurrentParameterSets)},
		{"opt_method_name",3,2,1,1,kw_846,0.,0.,-3,N_mdm(str,subMethodName)},
		{"opt_method_pointer",3,0,1,1,0,0.,0.,-3,N_mdm(str,subMethodPointer)},
		{"processors_per_iterator",0x19,0,6,0,0,0.,0.,0,N_mdm(pint,procsPerIterator)},
		{"random_weight_sets",9,1,2,0,kw_847,0.,0.,0,N_mdm(int,concurrentRandomJobs)},
		{"weight_sets",14,0,3,0,0,0.,0.,0,N_mdm(RealDL,concurrentParameterSets)}
		},
	kw_849[4] = {
		{"model_pointer",11,0,4,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"partitions",13,0,1,0,0,0.,0.,0,N_mdm(usharray,varPartitions)},
		{"samples",9,0,2,0,0,0.,0.,0,N_mdm(int,numSamples)},
		{"seed",0x19,0,3,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_850[7] = {
		{"converge_order",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_CONVERGE_ORDER)},
		{"converge_qoi",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_CONVERGE_QOI)},
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"estimate_order",8,0,1,1,0,0.,0.,0,N_mdm(utype,subMethod_SUBMETHOD_ESTIMATE_ORDER)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,5,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"refinement_rate",10,0,2,0,0,0.,0.,0,N_mdm(Real,refinementRate)}
		},
	kw_851[7] = {
		{"constraint_tolerance",10,0,4,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"gradient_tolerance",10,0,3,0,0,0.,0.,0,N_mdm(Real,gradientTolerance)},
		{"max_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,7,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"options_file",11,0,5,0,0,0.,0.,0,N_mdm(str,advancedOptionsFilename)},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"variable_tolerance",10,0,2,0,0,0.,0.,0,N_mdm(Real,threshStepLength)}
		},
	kw_852[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_853[6] = {
		{"max_function_evaluations",0x29,0,4,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,3,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,6,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"scaling",8,0,5,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,1,0,0,0.,0.,0,N_mdm(pint,randomSeed)},
		{"trust_region",8,6,2,0,kw_852,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final}
		},
	kw_854[2] = {
		{"num_generations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,1,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_855[2] = {
		{"num_generations",0x29,0,2,0,0,0.,0.,0,N_mdm(sizet,numGenerations)},
		{"percent_change",10,0,1,0,0,0.,0.,0,N_mdm(Realz,convergenceTolerance)}
		},
	kw_856[2] = {
		{"average_fitness_tracker",8,2,1,1,kw_854,0.,0.,0,N_mdm(lit,convergenceType_average_fitness_tracker)},
		{"best_fitness_tracker",8,2,1,1,kw_855,0.,0.,0,N_mdm(lit,convergenceType_best_fitness_tracker)}
		},
	kw_857[2] = {
		{"num_offspring",0x19,0,2,0,0,0.,0.,0,N_mdm(pintz,numOffspring)},
		{"num_parents",0x19,0,1,0,0,0.,0.,0,N_mdm(pintz,numParents)}
		},
	kw_858[5] = {
		{"crossover_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_crossoverType_null_crossover)},
		{"multi_point_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_binary)},
		{"multi_point_parameterized_binary",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_parameterized_binary)},
		{"multi_point_real",9,0,1,1,0,0.,0.,0,N_mdm(ilit2p,TYPE_DATA_crossoverType_multi_point_real)},
		{"shuffle_random",8,2,1,1,kw_857,0.,0.,0,N_mdm(litc,TYPE_DATA_crossoverType_shuffle_random)}
		},
	kw_859[2] = {
		{"constraint_penalty",10,0,2,0,0,0.,0.,0,N_mdm(Realp,constraintTolerance)},
		{"merit_function",8,0,1,1,0,0.,0.,0,N_mdm(lit,fitnessType_merit_function)}
		},
	kw_860[3] = {
		{"flat_file",11,0,1,1,0,0.,0.,0,N_mdm(slit2,TYPE_DATA_initializationType_flat_file)},
		{"simple_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_random)},
		{"unique_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,initializationType_unique_random)}
		},
	kw_861[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_862[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_863[1] = {
		{"mutation_scale",10,0,1,0,0,0.,0.,0,N_mdm(Real01,mutationScale)}
		},
	kw_864[6] = {
		{"bit_random",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_bit_random)},
		{"mutation_rate",10,0,2,0,0,0.,0.,0,N_mdm(litz,TYPE_DATA_mutationType_null_mutation)},
		{"offset_cauchy",8,1,1,1,kw_861,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_cauchy)},
		{"offset_normal",8,1,1,1,kw_862,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_normal)},
		{"offset_uniform",8,1,1,1,kw_863,0.,0.,0,N_mdm(litc,TYPE_DATA_mutationType_offset_uniform)},
		{"replace_uniform",8,0,1,1,0,0.,0.,0,N_mdm(lit,mutationType_replace_uniform)}
		},
	kw_865[4] = {
		{"elitist",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_elitist)},
		{"favor_feasible",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_favor_feasible)},
		{"roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_roulette_wheel)},
		{"unique_roulette_wheel",8,0,1,1,0,0.,0.,0,N_mdm(lit,replacementType_unique_roulette_wheel)}
		},
	kw_866[15] = {
		{"convergence_tolerance",10,0,14,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"convergence_type",8,2,3,0,kw_856},
		{"crossover_type",8,5,11,0,kw_858},
		{"fitness_type",8,2,1,0,kw_859},
		{"initialization_type",8,3,10,0,kw_860},
		{"log_file",11,0,8,0,0,0.,0.,0,N_mdm(str,logFile)},
		{"max_function_evaluations",0x29,0,5,0,0,0.,0.,0,N_mdm(nnint,maxFunctionEvaluations)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"model_pointer",11,0,15,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"mutation_type",8,6,12,0,kw_864},
		{"population_size",0x29,0,7,0,0,0.,0.,0,N_mdm(nnint,populationSize)},
		{"print_each_pop",8,0,9,0,0,0.,0.,0,N_mdm(true,printPopFlag)},
		{"replacement_type",8,4,2,0,kw_865},
		{"scaling",8,0,6,0,0,0.,0.,0,N_mdm(true,methodScaling)},
		{"seed",0x19,0,13,0,0,0.,0.,0,N_mdm(pint,randomSeed)}
		},
	kw_867[8] = {
		{"approx_method_name",3,0,1,1,0,0.,0.,4,N_mdm(str,subMethodName)},
		{"approx_method_pointer",3,0,1,1,0,0.,0.,4,N_mdm(str,subMethodPointer)},
		{"approx_model_pointer",3,0,2,2,0,0.,0.,4,N_mdm(str,modelPointer)},
		{"max_iterations",0x29,0,4,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"method_name",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"model_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"replace_points",8,0,3,0,0,0.,0.,0,N_mdm(true,surrBasedGlobalReplacePts)}
		},
	kw_868[2] = {
		{"filter",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalAcceptLogic_FILTER)},
		{"tr_ratio",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalAcceptLogic_TR_RATIO)}
		},
	kw_869[7] = {
		{"augmented_lagrangian_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_AUGMENTED_LAGRANGIAN_OBJECTIVE)},
		{"lagrangian_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_LAGRANGIAN_OBJECTIVE)},
		{"linearized_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_LINEARIZED_CONSTRAINTS)},
		{"no_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_NO_CONSTRAINTS)},
		{"original_constraints",8,0,2,2,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbCon_ORIGINAL_CONSTRAINTS)},
		{"original_primary",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_ORIGINAL_PRIMARY)},
		{"single_objective",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalSubProbObj_SINGLE_OBJECTIVE)}
		},
	kw_870[1] = {
		{"homotopy",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalConstrRelax_HOMOTOPY)}
		},
	kw_871[4] = {
		{"adaptive_penalty_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_ADAPTIVE_PENALTY_MERIT)},
		{"augmented_lagrangian_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_AUGMENTED_LAGRANGIAN_MERIT)},
		{"lagrangian_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_LAGRANGIAN_MERIT)},
		{"penalty_merit",8,0,1,1,0,0.,0.,0,N_mdm(type,surrBasedLocalMeritFn_PENALTY_MERIT)}
		},
	kw_872[6] = {
		{"contract_threshold",10,0,3,0,0,0.,0.,0,N_mdm(Real,trustRegionContractTrigger)},
		{"contraction_factor",10,0,5,0,0,0.,0.,0,N_mdm(Real,trustRegionContract)},
		{"expand_threshold",10,0,4,0,0,0.,0.,0,N_mdm(Real,trustRegionExpandTrigger)},
		{"expansion_factor",10,0,6,0,0,0.,0.,0,N_mdm(Real,trustRegionExpand)},
		{"initial_size",14,0,1,0,0,0.,0.,0,N_mdm(RealDL,trustRegionInitSize)},
		{"minimum_size",10,0,2,0,0,0.,0.,0,N_mdm(Real,trustRegionMinSize)}
		},
	kw_873[16] = {
		{"acceptance_logic",8,2,7,0,kw_868},
		{"approx_method_name",3,0,1,1,0,0.,0.,9,N_mdm(str,subMethodName)},
		{"approx_method_pointer",3,0,1,1,0,0.,0.,9,N_mdm(str,subMethodPointer)},
		{"approx_model_pointer",3,0,2,2,0,0.,0.,9,N_mdm(str,modelPointer)},
		{"approx_subproblem",8,7,5,0,kw_869},
		{"constraint_relax",8,1,8,0,kw_870},
		{"constraint_tolerance",10,0,12,0,0,0.,0.,0,N_mdm(Real,constraintTolerance)},
		{"convergence_tolerance",10,0,11,0,0,0.,0.,0,N_mdm(Real,convergenceTolerance)},
		{"max_iterations",0x29,0,10,0,0,0.,0.,0,N_mdm(nnint,maxIterations)},
		{"merit_function",8,4,6,0,kw_871},
		{"method_name",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodName)},
		{"method_pointer",11,0,1,1,0,0.,0.,0,N_mdm(str,subMethodPointer)},
		{"model_pointer",11,0,2,2,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"soft_convergence_limit",9,0,3,0,0,0.,0.,0,N_mdm(ushint,softConvLimit)},
		{"trust_region",8,6,9,0,kw_872,0.,0.,0,0,0,NIDRProblemDescDB::method_tr_final},
		{"truth_surrogate_bypass",8,0,4,0,0,0.,0.,0,N_mdm(true,surrBasedLocalLayerBypass)}
		},
	kw_874[4] = {
		{"final_point",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,finalPoint)},
		{"model_pointer",11,0,3,0,0,0.,0.,0,N_mdm(str,modelPointer)},
		{"num_steps",9,0,2,2,0,0.,0.,0,N_mdm(int,numSteps)},
		{"step_vector",14,0,1,1,0,0.,0.,0,N_mdm(RealDL,stepVector)}
		},
	kw_875[92] = {
		{"adaptive_sampling",8,19,4,1,kw_46,0.,0.,0,N_mdm(utype,methodName_ADAPTIVE_SAMPLING)},
		{"asynch_pattern_search",8,13,4,1,kw_49,0.,0.,0,N_mdm(utype,methodName_ASYNCH_PATTERN_SEARCH)},
		{"bayes_calibration",8,17,4,1,kw_416,0.,0.,0,N_mdm(utype,methodName_BAYES_CALIBRATION)},
		{"branch_and_bound",8,3,4,1,kw_418,0.,0.,0,N_mdm(utype,methodName_BRANCH_AND_BOUND)},
		{"centered_parameter_study",8,4,4,1,kw_419,0.,0.,0,N_mdm(utype,methodName_CENTERED_PARAMETER_STUDY)},
		{"coliny_apps",0,13,4,1,kw_49,0.,0.,-4,N_mdm(utype,methodName_ASYNCH_PATTERN_SEARCH)},
		{"coliny_beta",8,11,4,1,kw_420,0.,0.,0,N_mdm(utype,methodName_COLINY_BETA)},
		{"coliny_cobyla",8,12,4,1,kw_421,0.,0.,0,N_mdm(utype,methodName_COLINY_COBYLA)},
		{"coliny_direct",8,16,4,1,kw_423,0.,0.,0,N_mdm(utype,methodName_COLINY_DIRECT)},
		{"coliny_ea",8,19,4,1,kw_432,0.,0.,0,N_mdm(utype,methodName_COLINY_EA)},
		{"coliny_pattern_search",8,22,4,1,kw_436,0.,0.,0,N_mdm(utype,methodName_COLINY_PATTERN_SEARCH)},
		{"coliny_solis_wets",8,18,4,1,kw_437,0.,0.,0,N_mdm(utype,methodName_COLINY_SOLIS_WETS)},
		{"conmin_frcg",8,7,4,1,kw_438,0.,0.,0,N_mdm(utype,methodName_CONMIN_FRCG)},
		{"conmin_mfd",8,7,4,1,kw_439,0.,0.,0,N_mdm(utype,methodName_CONMIN_MFD)},
		{"dace",8,15,4,1,kw_441,0.,0.,0,N_mdm(utype,methodName_DACE)},
		{"demo_tpl",8,7,4,1,kw_442,0.,0.,0,N_mdm(utype,methodName_DEMO_TPL)},
		{"dl_solver",11,3,4,1,kw_443,0.,0.,0,N_mdm(utype_lit,TYPE_DATA_methodName_DL_SOLVER)},
		{"dot_bfgs",8,7,4,1,kw_444,0.,0.,0,N_mdm(utype,methodName_DOT_BFGS)},
		{"dot_frcg",8,7,4,1,kw_445,0.,0.,0,N_mdm(utype,methodName_DOT_FRCG)},
		{"dot_mmfd",8,7,4,1,kw_446,0.,0.,0,N_mdm(utype,methodName_DOT_MMFD)},
		{"dot_slp",8,7,4,1,kw_447,0.,0.,0,N_mdm(utype,methodName_DOT_SLP)},
		{"dot_sqp",8,7,4,1,kw_448,0.,0.,0,N_mdm(utype,methodName_DOT_SQP)},
		{"efficient_global",8,13,4,1,kw_454,0.,0.,0,N_mdm(utype,methodName_EFFICIENT_GLOBAL)},
		{"final_solutions",0x29,0,3,0,0,0.,0.,0,N_mdm(sizet,numFinalSolutions)},
		{"fsu_cvt",8,10,4,1,kw_457,0.,0.,0,N_mdm(utype,methodName_FSU_CVT)},
		{"fsu_quasi_mc",8,12,4,1,kw_459},
		{"function_train_uq",8,22,4,1,kw_476,0.,0.,0,N_mdm(utype,methodName_C3_FUNCTION_TRAIN)},
		{"gaussian_process_adaptive_importance_sampling",0,15,4,1,kw_488,0.,0.,6,N_mdm(utype,methodName_GPAIS)},
		{"genie_direct",8,4,4,1,kw_489,0.,0.,0,N_mdm(utype,methodName_GENIE_DIRECT)},
		{"genie_opt_darts",8,4,4,1,kw_490,0.,0.,0,N_mdm(utype,methodName_GENIE_OPT_DARTS)},
		{"global_evidence",8,12,4,1,kw_510,0.,0.,0,N_mdm(utype,methodName_GLOBAL_EVIDENCE)},
		{"global_interval_est",8,11,4,1,kw_524,0.,0.,0,N_mdm(utype,methodName_GLOBAL_INTERVAL_EST)},
		{"global_reliability",8,21,4,1,kw_536,0.,0.,0,N_mdm(utype,methodName_GLOBAL_RELIABILITY)},
		{"gpais",8,15,4,1,kw_488,0.,0.,0,N_mdm(utype,methodName_GPAIS)},
		{"hybrid",8,5,4,1,kw_547,0.,0.,0,N_mdm(utype,methodName_HYBRID)},
		{"id_method",11,0,1,0,0,0.,0.,0,N_mdm(str,idMethod)},
		{"importance_sampling",8,15,4,1,kw_555,0.,0.,0,N_mdm(utype,methodName_IMPORTANCE_SAMPLING)},
		{"list_parameter_study",8,3,4,1,kw_558,0.,0.,0,N_mdm(utype,methodName_LIST_PARAMETER_STUDY)},
		{"local_evidence",8,7,4,1,kw_565,0.,0.,0,N_mdm(utype,methodName_LOCAL_EVIDENCE)},
		{"local_interval_est",8,4,4,1,kw_566,0.,0.,0,N_mdm(utype,methodName_LOCAL_INTERVAL_EST)},
		{"local_reliability",8,10,4,1,kw_578,0.,0.,0,N_mdm(utype,methodName_LOCAL_RELIABILITY)},
		{"mesh_adaptive_search",8,14,4,1,kw_580,0.,0.,0,N_mdm(utype,methodName_MESH_ADAPTIVE_SEARCH)},
		{"moga",8,17,4,1,kw_595,0.,0.,0,N_mdm(utype,methodName_MOGA)},
		{"multi_start",8,7,4,1,kw_599,0.,0.,0,N_mdm(utype,methodName_MULTI_START)},
		{"multidim_parameter_study",8,2,4,1,kw_600,0.,0.,0,N_mdm(utype,methodName_MULTIDIM_PARAMETER_STUDY)},
		{"multifidelity_polynomial_chaos",8,36,4,1,kw_646,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_POLYNOMIAL_CHAOS)},
		{"multifidelity_stoch_collocation",8,33,4,1,kw_671,0.,0.,0,N_mdm(utype,methodName_MULTIFIDELITY_STOCH_COLLOCATION)},
		{"multilevel_mc",0,12,4,1,kw_678,0.,0.,2,N_mdm(utype,methodName_MULTILEVEL_SAMPLING)},
		{"multilevel_polynomial_chaos",8,35,4,1,kw_721,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_POLYNOMIAL_CHAOS)},
		{"multilevel_sampling",8,12,4,1,kw_678,0.,0.,0,N_mdm(utype,methodName_MULTILEVEL_SAMPLING)},
		{"ncsu_direct",8,9,4,1,kw_722,0.,0.,0,N_mdm(utype,methodName_NCSU_DIRECT)},
		{"nl2sol",8,15,4,1,kw_723,0.,0.,0,N_mdm(utype,methodName_NL2SOL)},
		{"nlpql_sqp",8,5,4,1,kw_724,0.,0.,0,N_mdm(utype,methodName_NLPQL_SQP)},
		{"nlssol_sqp",8,10,4,1,kw_725,0.,0.,0,N_mdm(utype,methodName_NLSSOL_SQP)},
		{"nond_adaptive_sampling",0,19,4,1,kw_46,0.,0.,-54,N_mdm(utype,methodName_ADAPTIVE_SAMPLING)},
		{"nond_bayes_calibration",0,17,4,1,kw_416,0.,0.,-53,N_mdm(utype,methodName_BAYES_CALIBRATION)},
		{"nond_global_evidence",0,12,4,1,kw_510,0.,0.,-26,N_mdm(utype,methodName_GLOBAL_EVIDENCE)},
		{"nond_global_interval_est",0,11,4,1,kw_524,0.,0.,-26,N_mdm(utype,methodName_GLOBAL_INTERVAL_EST)},
		{"nond_global_reliability",0,21,4,1,kw_536,0.,0.,-26,N_mdm(utype,methodName_GLOBAL_RELIABILITY)},
		{"nond_importance_sampling",0,15,4,1,kw_555,0.,0.,-23,N_mdm(utype,methodName_IMPORTANCE_SAMPLING)},
		{"nond_local_evidence",0,7,4,1,kw_565,0.,0.,-22,N_mdm(utype,methodName_LOCAL_EVIDENCE)},
		{"nond_local_interval_est",0,4,4,1,kw_566,0.,0.,-22,N_mdm(utype,methodName_LOCAL_INTERVAL_EST)},
		{"nond_local_reliability",0,10,4,1,kw_578,0.,0.,-22,N_mdm(utype,methodName_LOCAL_RELIABILITY)},
		{"nond_pof_darts",0,11,4,1,kw_734,0.,0.,16,N_mdm(utype,methodName_POF_DARTS)},
		{"nond_polynomial_chaos",0,36,4,1,kw_778,0.,0.,16,N_mdm(utype,methodName_POLYNOMIAL_CHAOS)},
		{"nond_rkd_darts",0,11,4,1,kw_787,0.,0.,18,N_mdm(utype,methodName_RKD_DARTS)},
		{"nond_sampling",0,19,4,1,kw_802,0.,0.,19,N_mdm(utype,methodName_RANDOM_SAMPLING)},
		{"nond_stoch_collocation",0,31,4,1,kw_825,0.,0.,21,N_mdm(utype,methodName_STOCH_COLLOCATION)},
		{"nonlinear_cg",8,5,4,1,kw_826,0.,0.,0,N_mdm(utype,methodName_NONLINEAR_CG)},
		{"nowpac",8,5,4,1,kw_828,0.,0.,0,N_mdm(utype,methodName_MIT_NOWPAC)},
		{"npsol_sqp",8,10,4,1,kw_829,0.,0.,0,N_mdm(utype,methodName_NPSOL_SQP)},
		{"optpp_cg",8,8,4,1,kw_830,0.,0.,0,N_mdm(utype,methodName_OPTPP_CG)},
		{"optpp_fd_newton",8,12,4,1,kw_833,0.,0.,0,N_mdm(utype,methodName_OPTPP_FD_NEWTON)},
		{"optpp_g_newton",8,12,4,1,kw_836,0.,0.,0,N_mdm(utype,methodName_OPTPP_G_NEWTON)},
		{"optpp_newton",8,12,4,1,kw_839,0.,0.,0,N_mdm(utype,methodName_OPTPP_NEWTON)},
		{"optpp_pds",8,6,4,1,kw_840,0.,0.,0,N_mdm(utype,methodName_OPTPP_PDS)},
		{"optpp_q_newton",8,12,4,1,kw_843,0.,0.,0,N_mdm(utype,methodName_OPTPP_Q_NEWTON)},
		{"output",8,5,2,0,kw_844},
		{"pareto_set",8,10,4,1,kw_848,0.,0.,0,N_mdm(utype,methodName_PARETO_SET)},
		{"pof_darts",8,11,4,1,kw_734,0.,0.,0,N_mdm(utype,methodName_POF_DARTS)},
		{"polynomial_chaos",8,36,4,1,kw_778,0.,0.,0,N_mdm(utype,methodName_POLYNOMIAL_CHAOS)},
		{"psuade_moat",8,4,4,1,kw_849,0.,0.,0,N_mdm(utype,methodName_PSUADE_MOAT)},
		{"richardson_extrap",8,7,4,1,kw_850,0.,0.,0,N_mdm(utype,methodName_RICHARDSON_EXTRAP)},
		{"rkd_darts",8,11,4,1,kw_787,0.,0.,0,N_mdm(utype,methodName_RKD_DARTS)},
		{"rol",8,7,4,1,kw_851,0.,0.,0,N_mdm(utype,methodName_ROL)},
		{"sampling",8,19,4,1,kw_802,0.,0.,0,N_mdm(utype,methodName_RANDOM_SAMPLING)},
		{"snowpac",8,6,4,1,kw_853,0.,0.,0,N_mdm(utype,methodName_MIT_SNOWPAC)},
		{"soga",8,15,4,1,kw_866,0.,0.,0,N_mdm(utype,methodName_SOGA)},
		{"stoch_collocation",8,31,4,1,kw_825,0.,0.,0,N_mdm(utype,methodName_STOCH_COLLOCATION)},
		{"surrogate_based_global",8,8,4,1,kw_867,0.,0.,0,N_mdm(utype,methodName_SURROGATE_BASED_GLOBAL)},
		{"surrogate_based_local",8,16,4,1,kw_873,0.,0.,0,N_mdm(utype,methodName_SURROGATE_BASED_LOCAL)},
		{"vector_parameter_study",8,4,4,1,kw_874,0.,0.,0,N_mdm(utype,methodName_VECTOR_PARAMETER_STUDY)}
		},
	kw_876[1] = {
		{"refinement_samples",13,0,1,0,0,0.,0.,0,N_mom(ivec,refineSamples)}
		},
	kw_877[3] = {
		{"local_gradient",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_LOCAL_GRAD)},
		{"mean_gradient",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_MEAN_GRAD)},
		{"mean_value",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceNormalization_SUBSPACE_NORM_MEAN_VALUE)}
		},
	kw_878[2] = {
		{"lhs",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceSampleType_SUBMETHOD_LHS)},
		{"random",8,0,1,1,0,0.,0.,0,N_mom(utype,subspaceSampleType_SUBMETHOD_RANDOM)}
		},
	kw_879[7] = {
		{"decrease",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_DECREASE_TOLERANCE)},
		{"decrease_tolerance",10,0,3,0,0,0.,0.,0,N_mom(Real,decreaseTolerance)},
		{"exhaustive",8,0,5,0,0,0.,0.,0,N_mom(false,subspaceCVIncremental)},
		{"max_rank",9,0,4,0,0,0.,0.,0,N_mom(int,subspaceCVMaxRank)},
		{"minimum",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_MINIMUM_METRIC)},
		{"relative",8,0,1,0,0,0.,0.,0,N_mom(utype,subspaceIdCVMethod_RELATIVE_TOLERANCE)},
		{"relative_tolerance",10,0,2,0,0,0.,0.,0,N_mom(Real,relTolerance)}
		},
	kw_880[1] = {
		{"truncation_tolerance",10,0,1,0,0,0.,0.,0,N_mom(Real,truncationTolerance)}
		},
	kw_881[4] = {
		{"bing_li",8,0,1,0,0,0.,0.,0,N_mom(true,subspaceIdBingLi)},
		{"constantine",8,0,2,0,0,0.,0.,0,N_mom(true,subspaceIdConstantine)},
		{"cross_validation",8,7,4,0,kw_879,0.,0.,0,N_mom(true,subspaceIdCV)},
		{"energy",8,1,3,0,kw_880,0.,0.,0,N_mom(true,subspaceIdEnergy)}
		},
	kw_882[8] = {
		{"actual_model_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"bootstrap_samples",9,0,6,0,0,0.,0.,0,N_mom(int,numReplicates)},
		{"build_surrogate",8,1,7,0,kw_876,0.,0.,0,N_mom(true,subspaceBuildSurrogate)},
		{"dimension",9,0,5,0,0,0.,0.,0,N_mom(int,subspaceDimension)},
		{"initial_samples",9,0,2,0,0,0.,0.,0,N_mom(int,initialSamples)},
		{"normalization",8,3,8,0,kw_877},
		{"sample_type",8,2,3,0,kw_878},
		{"truncation_method",8,4,4,0,kw_881}
		},
	kw_883[1] = {
		{"collocation_ratio",10,0,1,1,0,0.,0.,0,N_mom(Real,adaptedBasisCollocRatio)}
		},
	kw_884[3] = {
		{"actual_model_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"expansion_order",9,1,2,2,kw_883,0.,0.,0,N_mom(int,adaptedBasisExpOrder)},
		{"sparse_grid_level",9,0,2,2,0,0.,0.,0,N_mom(int,adaptedBasisSparseGridLev)}
		},
	kw_885[1] = {
		{"optional_interface_responses_pointer",11,0,1,0,0,0.,0.,0,N_mom(str,optionalInterfRespPointer)}
		},
	kw_886[2] = {
		{"master",8,0,1,1,0,0.,0.,0,N_mom(type,subMethodScheduling_MASTER_SCHEDULING)},
		{"peer",8,0,1,1,0,0.,0.,0,N_mom(type,subMethodScheduling_PEER_SCHEDULING)}
		},
	kw_887[8] = {
		{"identity_response_mapping",8,0,8,0,0,0.,0.,0,N_mom(true,identityRespMap)},
		{"iterator_scheduling",8,2,2,0,kw_886},
		{"iterator_servers",0x19,0,1,0,0,0.,0.,0,N_mom(pint,subMethodServers)},
		{"primary_response_mapping",14,0,6,0,0,0.,0.,0,N_mom(RealDL,primaryRespCoeffs)},
		{"primary_variable_mapping",15,0,4,0,0,0.,0.,0,N_mom(strL,primaryVarMaps)},
		{"processors_per_iterator",0x19,0,3,0,0,0.,0.,0,N_mom(pint,subMethodProcs)},
		{"secondary_response_mapping",14,0,7,0,0,0.,0.,0,N_mom(RealDL,secondaryRespCoeffs)},
		{"secondary_variable_mapping",15,0,5,0,0,0.,0.,0,N_mom(strL,secondaryVarMaps)}
		},
	kw_888[2] = {
		{"optional_interface_pointer",11,1,1,0,kw_885,0.,0.,0,N_mom(str,interfacePointer)},
		{"sub_method_pointer",11,8,2,1,kw_887,0.,0.,0,N_mom(str,subMethodPointer)}
		},
	kw_889[2] = {
		{"exponential",8,0,1,1,0,0.,0.,0,N_mom(utype,analyticCovIdForm_EXP_L1)},
		{"squared_exponential",8,0,1,1,0,0.,0.,0,N_mom(utype,analyticCovIdForm_EXP_L2)}
		},
	kw_890[3] = {
		{"analytic_covariance",8,2,1,1,kw_889},
		{"dace_method_pointer",11,0,1,1,0,0.,0.,0,N_mom(str,subMethodPointer)},
		{"rf_data_file",11,0,1,1,0,0.,0.,0,N_mom(str,rfDataFileName)}
		},
	kw_891[2] = {
		{"karhunen_loeve",8,0,1,1,0,0.,0.,0,N_mom(utype,randomFieldIdForm_RF_KARHUNEN_LOEVE)},
		{"principal_components",8,0,1,1,0,0.,0.,0,N_mom(utype,randomFieldIdForm_RF_PCA_GP)}
		},
	kw_892[5] = {
		{"build_source",8,3,1,0,kw_890},
		{"expansion_bases",9,0,3,0,0,0.,0.,0,N_mom(int,subspaceDimension)},
		{"expansion_form",8,2,2,0,kw_891},
		{"propagation_model_pointer",11,0,5,1,0,0.,0.,0,N_mom(str,propagationModelPointer)},
		{"truncation_tolerance",10,0,4,0,0,0.,0.,0,N_mom(Real,truncationTolerance)}
		},
	kw_893[1] = {
		{"solution_level_control",11,0,1,0,0,0.,0.,0,N_mom(str,solutionLevelControl)}
		},
	kw_894[2] = {
		{"interface_pointer",11,0,1,0,0,0.,0.,0,N_mom(str,interfacePointer)},
		{"solution_level_cost",14,1,2,0,kw_893,0.,0.,0,N_mom(RealDL,solutionLevelCost)}
		},
	kw_895[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importChalUseVariableLabels)}
		},
	kw_896[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importChalUseVariableLabels)}
		},
	kw_897[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_EVAL_ID)},
		{"header",8,1,1,0,kw_896,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,importChallengeFormat_TABULAR_IFACE_ID)}
		},
	kw_898[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mom(true,importChallengeActive)},
		{"annotated",8,1,1,0,kw_895,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_897,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,importChallengeFormat_TABULAR_NONE)}
		},
	kw_899[6] = {
		{"additive",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_ADDITIVE_CORRECTION)},
		{"combined",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_COMBINED_CORRECTION)},
		{"first_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_1)},
		{"multiplicative",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_MULTIPLICATIVE_CORRECTION)},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_2)},
		{"zeroth_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_0)}
		},
	kw_900[1] = {
		{"folds",0x19,0,1,0,0,0.,0.,0,N_mom(int,refineCVFolds)}
		},
	kw_901[5] = {
		{"convergence_tolerance",10,0,3,0,0,0.,0.,0,N_mom(Real,convergenceTolerance)},
		{"cross_validation_metric",11,1,5,0,kw_900,0.,0.,0,N_mom(str,refineCVMetric)},
		{"max_function_evaluations",0x19,0,2,0,0,0.,0.,0,N_mom(int,maxFunctionEvals)},
		{"max_iterations",0x19,0,1,0,0,0.,0.,0,N_mom(int,maxIterations)},
		{"soft_convergence_limit",0x29,0,4,0,0,0.,0.,0,N_mom(int,softConvergenceLimit)}
		},
	kw_902[1] = {
		{"auto_refinement",8,5,1,0,kw_901,0.,0.,0,N_mom(true,autoRefine)}
		},
	kw_903[2] = {
		{"folds",9,0,1,0,0,0.,0.,0,N_mom(int,numFolds)},
		{"percent",10,0,1,0,0,0.,0.,0,N_mom(Real,percentFold)}
		},
	kw_904[2] = {
		{"cross_validation",8,2,1,0,kw_903,0.,0.,0,N_mom(true,crossValidateFlag)},
		{"press",8,0,2,0,0,0.,0.,0,N_mom(true,pressFlag)}
		},
	kw_905[2] = {
		{"gradient_threshold",10,0,1,1,0,0.,0.,0,N_mom(Real,discontGradThresh)},
		{"jump_threshold",10,0,1,1,0,0.,0.,0,N_mom(Real,discontJumpThresh)}
		},
	kw_906[3] = {
		{"cell_type",11,0,1,0,0,0.,0.,0,N_mom(str,decompCellType)},
		{"discontinuity_detection",8,2,3,0,kw_905,0.,0.,0,N_mom(true,decompDiscontDetect)},
		{"support_layers",9,0,2,0,0,0.,0.,0,N_mom(int,decompSupportLayers)}
		},
	kw_907[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,exportApproxFormat_TABULAR_IFACE_ID)}
		},
	kw_908[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_907,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,exportApproxFormat_TABULAR_NONE)}
		},
	kw_909[10] = {
		{"adapt_rank",8,0,9,0,0,0.,0.,0,N_mom(true,adaptRank)},
		{"kick_rank",0x29,0,7,0,0,0.,0.,0,N_mom(sizet,kickRank)},
		{"max_cross_iterations",0x29,0,10,0,0,0.,0.,0,N_mom(sizet,crossMaxIter)},
		{"max_order",0x29,0,5,0,0,0.,0.,0,N_mom(sizet,maxOrder)},
		{"max_rank",0x29,0,8,0,0,0.,0.,0,N_mom(sizet,maxRank)},
		{"max_solver_iterations",0x29,0,1,0,0,0.,0.,0,N_mdm(nnint,maxSolverIterations)},
		{"rounding_tolerance",10,0,3,0,0,0.,0.,0,N_mom(Real,roundingTolerance)},
		{"solver_tolerance",10,0,2,0,0,0.,0.,0,N_mom(Real,solverTolerance)},
		{"start_order",0x29,0,4,0,0,0.,0.,0,N_mom(sizet,startOrder)},
		{"start_rank",0x29,0,6,0,0,0.,0.,0,N_mom(sizet,startRank)}
		},
	kw_910[3] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_constant)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_linear)},
		{"reduced_quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_reduced_quadratic)}
		},
	kw_911[2] = {
		{"point_selection",8,0,1,0,0,0.,0.,0,N_mom(true,pointSelection)},
		{"trend",8,3,2,0,kw_910}
		},
	kw_912[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_913[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_912}
		},
	kw_914[4] = {
		{"constant",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_constant)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_linear)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_quadratic)},
		{"reduced_quadratic",8,0,1,1,0,0.,0.,0,N_mom(lit,trendOrder_reduced_quadratic)}
		},
	kw_915[7] = {
		{"correlation_lengths",14,0,5,0,0,0.,0.,0,N_mom(RealDL,krigingCorrelations)},
		{"export_model",8,2,6,0,kw_913,0.,0.,0,N_mom(true,exportSurrogate)},
		{"find_nugget",9,0,4,0,0,0.,0.,0,N_mom(shint,krigingFindNugget)},
		{"max_trials",0x19,0,3,0,0,0.,0.,0,N_mom(shint,krigingMaxTrials)},
		{"nugget",0x1a,0,4,0,0,0.,0.,0,N_mom(Real,krigingNugget)},
		{"optimization_method",11,0,2,0,0,0.,0.,0,N_mom(str,krigingOptMethod)},
		{"trend",8,4,1,0,kw_914}
		},
	kw_916[2] = {
		{"dakota",8,2,1,1,kw_911,0.,0.,0,N_mom(lit,surrogateType_global_gaussian)},
		{"surfpack",8,7,1,1,kw_915,0.,0.,0,N_mom(lit,surrogateType_global_kriging)}
		},
	kw_917[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importUseVariableLabels)}
		},
	kw_918[1] = {
		{"use_variable_labels",8,0,1,0,0,0.,0.,0,N_mom(true,importUseVariableLabels)}
		},
	kw_919[3] = {
		{"eval_id",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_EVAL_ID)},
		{"header",8,1,1,0,kw_918,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_HEADER)},
		{"interface_id",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,importBuildFormat_TABULAR_IFACE_ID)}
		},
	kw_920[4] = {
		{"active_only",8,0,2,0,0,0.,0.,0,N_mom(true,importBuildActive)},
		{"annotated",8,1,1,0,kw_917,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_ANNOTATED)},
		{"custom_annotated",8,3,1,0,kw_919,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_mom(utype,importBuildFormat_TABULAR_NONE)}
		},
	kw_921[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_922[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_921}
		},
	kw_923[2] = {
		{"cubic",8,0,1,1,0,0.,0.,0,N_mom(lit,marsInterpolation_cubic)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(lit,marsInterpolation_linear)}
		},
	kw_924[3] = {
		{"export_model",8,2,3,0,kw_922,0.,0.,0,N_mom(true,exportSurrogate)},
		{"interpolation",8,2,2,0,kw_923},
		{"max_bases",9,0,1,0,0,0.,0.,0,N_mom(shint,marsMaxBases)}
		},
	kw_925[2] = {
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_926[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,2,2,1,kw_925}
		},
	kw_927[4] = {
		{"basis_order",0x29,0,1,0,0,0.,0.,0,N_mom(shint,polynomialOrder)},
		{"export_model",8,2,3,0,kw_926,0.,0.,0,N_mom(true,exportSurrogate)},
		{"poly_order",0x21,0,1,0,0,0.,0.,-2,N_mom(shint,polynomialOrder)},
		{"weight_function",9,0,2,0,0,0.,0.,0,N_mom(shint,mlsWeightFunction)}
		},
	kw_928[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_929[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_928}
		},
	kw_930[5] = {
		{"export_model",8,2,4,0,kw_929,0.,0.,0,N_mom(true,exportSurrogate)},
		{"max_nodes",9,0,1,0,0,0.,0.,0,N_mom(shint,annNodes)},
		{"nodes",1,0,1,0,0,0.,0.,-1,N_mom(shint,annNodes)},
		{"random_weight",9,0,3,0,0,0.,0.,0,N_mom(shint,annRandomWeight)},
		{"range",10,0,2,0,0,0.,0.,0,N_mom(Real,annRange)}
		},
	kw_931[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_932[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_931}
		},
	kw_933[5] = {
		{"basis_order",0x29,0,1,1,0,0.,0.,0,N_mom(shint,polynomialOrder)},
		{"cubic",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_3)},
		{"export_model",8,2,2,0,kw_932,0.,0.,0,N_mom(true,exportSurrogate)},
		{"linear",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_1)},
		{"quadratic",8,0,1,1,0,0.,0.,0,N_mom(order,polynomialOrder_2)}
		},
	kw_934[4] = {
		{"algebraic_console",8,0,4,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_CONSOLE)},
		{"algebraic_file",8,0,3,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_ALGEBRAIC_FILE)},
		{"binary_archive",8,0,2,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_BINARY_ARCHIVE)},
		{"text_archive",8,0,1,0,0,0.,0.,0,N_mom(augment_utype,modelExportFormat_TEXT_ARCHIVE)}
		},
	kw_935[2] = {
		{"filename_prefix",11,0,1,0,0,0.,0.,0,N_mom(str,modelExportPrefix)},
		{"formats",8,4,2,1,kw_934}
		},
	kw_936[5] = {
		{"bases",9,0,1,0,0,0.,0.,0,N_mom(shint,rbfBases)},
		{"export_model",8,2,5,0,kw_935,0.,0.,0,N_mom(true,exportSurrogate)},
		{"max_pts",9,0,2,0,0,0.,0.,0,N_mom(shint,rbfMaxPts)},
		{"max_subsets",9,0,4,0,0,0.,0.,0,N_mom(shint,rbfMaxSubsets)},
		{"min_partition",9,0,3,0,0,0.,0.,0,N_mom(shint,rbfMinPartition)}
		},
	kw_937[3] = {
		{"all",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_all)},
		{"none",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_none)},
		{"region",8,0,1,1,0,0.,0.,0,N_mom(lit,approxPointReuse_region)}
		},
	kw_938[27] = {
		{"actual_model_pointer",11,0,4,0,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"challenge_points_file",3,4,11,0,kw_898,0.,0.,10,N_mom(str,importChallengePtsFile)},
		{"correction",8,6,9,0,kw_899},
		{"dace_method_pointer",11,1,4,0,kw_902,0.,0.,0,N_mom(str,subMethodPointer)},
		{"diagnostics",7,2,10,0,kw_904,0.,0.,11,N_mom(strL,diagMetrics)},
		{"domain_decomposition",8,3,2,0,kw_906,0.,0.,0,N_mom(true,domainDecomp)},
		{"export_approx_points_file",11,3,7,0,kw_908,0.,0.,0,N_mom(str,exportApproxPtsFile)},
		{"export_points_file",3,3,7,0,kw_908,0.,0.,-1,N_mom(str,exportApproxPtsFile)},
		{"function_train",8,10,1,1,kw_909,0.,0.,0,N_mom(lit,surrogateType_global_function_train)},
		{"gaussian_process",8,2,1,1,kw_916},
		{"import_build_points_file",11,4,6,0,kw_920,0.,0.,0,N_mom(str,importBuildPtsFile)},
		{"import_challenge_points_file",11,4,11,0,kw_898,0.,0.,0,N_mom(str,importChallengePtsFile)},
		{"import_points_file",3,4,6,0,kw_920,0.,0.,-2,N_mom(str,importBuildPtsFile)},
		{"kriging",0,2,1,1,kw_916,0.,0.,-4},
		{"mars",8,3,1,1,kw_924,0.,0.,0,N_mom(lit,surrogateType_global_mars)},
		{"metrics",15,2,10,0,kw_904,0.,0.,0,N_mom(strL,diagMetrics)},
		{"minimum_points",8,0,3,0,0,0.,0.,0,N_mom(type,pointsManagement_MINIMUM_POINTS)},
		{"moving_least_squares",8,4,1,1,kw_927,0.,0.,0,N_mom(lit,surrogateType_global_moving_least_squares)},
		{"neural_network",8,5,1,1,kw_930,0.,0.,0,N_mom(lit,surrogateType_global_neural_network)},
		{"polynomial",8,5,1,1,kw_933,0.,0.,0,N_mom(lit,surrogateType_global_polynomial)},
		{"radial_basis",8,5,1,1,kw_936,0.,0.,0,N_mom(lit,surrogateType_global_radial_basis)},
		{"recommended_points",8,0,3,0,0,0.,0.,0,N_mom(type,pointsManagement_RECOMMENDED_POINTS)},
		{"reuse_points",8,3,5,0,kw_937},
		{"reuse_samples",0,3,5,0,kw_937,0.,0.,-1},
		{"samples_file",3,4,6,0,kw_920,0.,0.,-14,N_mom(str,importBuildPtsFile)},
		{"total_points",9,0,3,0,0,0.,0.,0,N_mom(int,pointsTotal)},
		{"use_derivatives",8,0,8,0,0,0.,0.,0,N_mom(true,modelUseDerivsFlag)}
		},
	kw_939[6] = {
		{"additive",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_ADDITIVE_CORRECTION)},
		{"combined",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_COMBINED_CORRECTION)},
		{"first_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_1)},
		{"multiplicative",8,0,2,2,0,0.,0.,0,N_mom(type,approxCorrectionType_MULTIPLICATIVE_CORRECTION)},
		{"second_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_2)},
		{"zeroth_order",8,0,1,1,0,0.,0.,0,N_mom(order,approxCorrectionOrder_0)}
		},
	kw_940[3] = {
		{"correction",8,6,2,0,kw_939},
		{"model_fidelity_sequence",7,0,1,1,0,0.,0.,1,N_mom(strL,orderedModelPointers)},
		{"ordered_model_fidelities",15,0,1,1,0,0.,0.,0,N_mom(strL,orderedModelPointers)}
		},
	kw_941[2] = {
		{"actual_model_pointer",11,0,2,2,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"taylor_series",8,0,1,1}
		},
	kw_942[3] = {
		{"actual_model_pointer",11,0,2,2,0,0.,0.,0,N_mom(str,actualModelPointer)},
		{"qmea",8,0,1,1,0,0.,0.,0,N_mom(lit,surrogateType_multipoint_qmea)},
		{"tana",8,0,1,1,0,0.,0.,0,N_mom(lit,surrogateType_multipoint_tana)}
		},
	kw_943[5] = {
		{"global",8,27,2,1,kw_938},
		{"hierarchical",8,3,2,1,kw_940,0.,0.,0,N_mom(lit,surrogateType_hierarchical)},
		{"id_surrogates",13,0,1,0,0,0.,0.,0,N_mom(intsetm1,surrogateFnIndices)},
		{"local",8,2,2,1,kw_941,0.,0.,0,N_mom(lit,surrogateType_local_taylor)},
		{"multipoint",8,3,2,1,kw_942}
		},
	kw_944[12] = {
		{"active_subspace",8,8,2,1,kw_882,0.,0.,0,N_mom(lit,modelType_active_subspace)},
		{"adapted_basis",8,3,2,1,kw_884,0.,0.,0,N_mom(lit,modelType_adapted_basis)},
		{"hierarchical_tagging",8,0,5,0,0,0.,0.,0,N_mom(true,hierarchicalTags)},
		{"id_model",11,0,1,0,0,0.,0.,0,N_mom(str,idModel)},
		{"nested",8,2,2,1,kw_888,0.,0.,0,N_mom(lit,modelType_nested)},
		{"random_field",8,5,2,1,kw_892,0.,0.,0,N_mom(lit,modelType_random_field)},
		{"responses_pointer",11,0,4,0,0,0.,0.,0,N_mom(str,responsesPointer)},
		{"simulation",0,2,2,1,kw_894,0.,0.,1,N_mom(lit,modelType_simulation)},
		{"single",8,2,2,1,kw_894,0.,0.,0,N_mom(lit,modelType_simulation)},
		{"subspace",0,8,2,1,kw_882,0.,0.,-9,N_mom(lit,modelType_active_subspace)},
		{"surrogate",8,5,2,1,kw_943,0.,0.,0,N_mom(lit,modelType_surrogate)},
		{"variables_pointer",11,0,3,0,0,0.,0.,0,N_mom(str,variablesPointer)}
		},
	kw_945[2] = {
		{"exp_id",8,0,2,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_HEADER)}
		},
	kw_946[3] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_EXPER_ANNOT)},
		{"custom_annotated",8,2,1,0,kw_945,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)}
		},
	kw_947[6] = {
		{"experiment_variance_type",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,varianceType)},
		{"interpolate",8,0,5,0,0,0.,0.,0,N_rem(true,interpolateFlag)},
		{"num_config_variables",0x29,0,2,0,0,0.,0.,0,N_rem(sizet,numExpConfigVars)},
		{"num_experiments",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numExperiments)},
		{"scalar_data_file",11,3,4,0,kw_946,0.,0.,0,N_rem(str,scalarDataFileName)},
		{"variance_type",0x807,0,3,0,0,0.,0.,-5,N_rem(strL,varianceType)}
		},
	kw_948[2] = {
		{"exp_id",8,0,2,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_EVAL_ID)},
		{"header",8,0,1,0,0,0.,0.,0,N_rem(augment_utype,scalarDataFormat_TABULAR_HEADER)}
		},
	kw_949[7] = {
		{"annotated",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_EXPER_ANNOT)},
		{"custom_annotated",8,2,1,0,kw_948,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"experiment_variance_type",0x80f,0,4,0,0,0.,0.,0,N_rem(strL,varianceType)},
		{"freeform",8,0,1,0,0,0.,0.,0,N_rem(utype,scalarDataFormat_TABULAR_NONE)},
		{"num_config_variables",0x29,0,3,0,0,0.,0.,0,N_rem(sizet,numExpConfigVars)},
		{"num_experiments",0x29,0,2,0,0,0.,0.,0,N_rem(sizet,numExperiments)},
		{"variance_type",0x807,0,4,0,0,0.,0.,-4,N_rem(strL,varianceType)}
		},
	kw_950[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_951[6] = {
		{"nonlinear_equality_scale_types",0x807,0,2,0,0,0.,0.,3,N_rem(strL,nonlinearEqScaleTypes)},
		{"nonlinear_equality_scales",0x806,0,3,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqScales)},
		{"nonlinear_equality_targets",6,0,1,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqTargets)},
		{"scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,nonlinearEqScaleTypes)},
		{"scales",0x80e,0,3,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqScales)},
		{"targets",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqTargets)}
		},
	kw_952[8] = {
		{"lower_bounds",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_lower_bounds",6,0,1,0,0,0.,0.,-1,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_scale_types",0x807,0,3,0,0,0.,0.,3,N_rem(strL,nonlinearIneqScaleTypes)},
		{"nonlinear_inequality_scales",0x806,0,4,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqScales)},
		{"nonlinear_inequality_upper_bounds",6,0,2,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqUpperBnds)},
		{"scale_types",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,nonlinearIneqScaleTypes)},
		{"scales",0x80e,0,4,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqScales)},
		{"upper_bounds",14,0,2,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqUpperBnds)}
		},
	kw_953[19] = {
		{"calibration_data",8,6,6,0,kw_947,0.,0.,0,N_rem(true,calibrationDataFlag)},
		{"calibration_data_file",11,7,6,0,kw_949,0.,0.,0,N_rem(str,scalarDataFileName)},
		{"calibration_term_scale_types",0x807,0,3,0,0,0.,0.,12,N_rem(strL,primaryRespFnScaleTypes)},
		{"calibration_term_scales",0x806,0,4,0,0,0.,0.,12,N_rem(RealDL,primaryRespFnScales)},
		{"calibration_weights",6,0,5,0,0,0.,0.,14,N_rem(RealDL,primaryRespFnWeights)},
		{"field_calibration_terms",0x29,3,2,0,kw_950,0.,0.,0,N_rem(sizet,numFieldLeastSqTerms)},
		{"least_squares_data_file",3,7,6,0,kw_949,0.,0.,-5,N_rem(str,scalarDataFileName)},
		{"least_squares_term_scale_types",0x807,0,3,0,0,0.,0.,7,N_rem(strL,primaryRespFnScaleTypes)},
		{"least_squares_term_scales",0x806,0,4,0,0,0.,0.,7,N_rem(RealDL,primaryRespFnScales)},
		{"least_squares_weights",6,0,5,0,0,0.,0.,9,N_rem(RealDL,primaryRespFnWeights)},
		{"nonlinear_equality_constraints",0x29,6,9,0,kw_951,0.,0.,0,N_rem(sizet,numNonlinearEqConstraints)},
		{"nonlinear_inequality_constraints",0x29,8,8,0,kw_952,0.,0.,0,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_nonlinear_equality_constraints",0x21,6,9,0,kw_951,0.,0.,-2,N_rem(sizet,numNonlinearEqConstraints)},
		{"num_nonlinear_inequality_constraints",0x21,8,8,0,kw_952,0.,0.,-2,N_rem(sizet,numNonlinearIneqConstraints)},
		{"primary_scale_types",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,primaryRespFnScaleTypes)},
		{"primary_scales",0x80e,0,4,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnScales)},
		{"scalar_calibration_terms",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numScalarLeastSqTerms)},
		{"simulation_variance",0x80e,0,7,0,0,0.,0.,0,N_rem(RealL,simVariance)},
		{"weights",14,0,5,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnWeights)}
		},
	kw_954[4] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_bounds)},
		{"ignore_bounds",8,0,1,0,0,0.,0.,0,N_rem(true,ignoreBounds)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_relative)}
		},
	kw_955[10] = {
		{"central",8,0,6,0,0,0.,0.,0,N_rem(lit,intervalType_central)},
		{"dakota",8,4,4,0,kw_954,0.,0.,0,N_rem(lit,methodSource_dakota)},
		{"fd_gradient_step_size",6,0,7,0,0,0.,0.,1,N_rem(RealL,fdGradStepSize)},
		{"fd_step_size",14,0,7,0,0,0.,0.,0,N_rem(RealL,fdGradStepSize)},
		{"forward",8,0,6,0,0,0.,0.,0,N_rem(lit,intervalType_forward)},
		{"id_analytic_gradients",13,0,2,2,0,0.,0.,0,N_rem(intset,idAnalyticGrads)},
		{"id_numerical_gradients",13,0,1,1,0,0.,0.,0,N_rem(intset,idNumericalGrads)},
		{"interval_type",8,0,5},
		{"method_source",8,0,3},
		{"vendor",8,0,4,0,0,0.,0.,0,N_rem(lit,methodSource_vendor)}
		},
	kw_956[2] = {
		{"fd_hessian_step_size",6,0,1,0,0,0.,0.,1,N_rem(RealL,fdHessStepSize)},
		{"fd_step_size",14,0,1,0,0,0.,0.,0,N_rem(RealL,fdHessStepSize)}
		},
	kw_957[1] = {
		{"damped",8,0,1,0,0,0.,0.,0,N_rem(lit,quasiHessianType_damped_bfgs)}
		},
	kw_958[2] = {
		{"bfgs",8,1,1,1,kw_957,0.,0.,0,N_rem(lit,quasiHessianType_bfgs)},
		{"sr1",8,0,1,1,0,0.,0.,0,N_rem(lit,quasiHessianType_sr1)}
		},
	kw_959[8] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_bounds)},
		{"central",8,0,3,0,0,0.,0.,0,N_rem(true,centralHess)},
		{"forward",8,0,3,0,0,0.,0.,0,N_rem(false,centralHess)},
		{"id_analytic_hessians",13,0,5,0,0,0.,0.,0,N_rem(intset,idAnalyticHessians)},
		{"id_numerical_hessians",13,2,1,0,kw_956,0.,0.,0,N_rem(intset,idNumericalHessians)},
		{"id_quasi_hessians",13,2,4,0,kw_958,0.,0.,0,N_rem(intset,idQuasiHessians)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_relative)}
		},
	kw_960[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_961[6] = {
		{"nonlinear_equality_scale_types",0x807,0,2,0,0,0.,0.,3,N_rem(strL,nonlinearEqScaleTypes)},
		{"nonlinear_equality_scales",0x806,0,3,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqScales)},
		{"nonlinear_equality_targets",6,0,1,0,0,0.,0.,3,N_rem(RealDL,nonlinearEqTargets)},
		{"scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,nonlinearEqScaleTypes)},
		{"scales",0x80e,0,3,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqScales)},
		{"targets",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearEqTargets)}
		},
	kw_962[8] = {
		{"lower_bounds",14,0,1,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_lower_bounds",6,0,1,0,0,0.,0.,-1,N_rem(RealDL,nonlinearIneqLowerBnds)},
		{"nonlinear_inequality_scale_types",0x807,0,3,0,0,0.,0.,3,N_rem(strL,nonlinearIneqScaleTypes)},
		{"nonlinear_inequality_scales",0x806,0,4,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqScales)},
		{"nonlinear_inequality_upper_bounds",6,0,2,0,0,0.,0.,3,N_rem(RealDL,nonlinearIneqUpperBnds)},
		{"scale_types",0x80f,0,3,0,0,0.,0.,0,N_rem(strL,nonlinearIneqScaleTypes)},
		{"scales",0x80e,0,4,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqScales)},
		{"upper_bounds",14,0,2,0,0,0.,0.,0,N_rem(RealDL,nonlinearIneqUpperBnds)}
		},
	kw_963[15] = {
		{"field_objectives",0x29,3,8,0,kw_960,0.,0.,0,N_rem(sizet,numFieldObjectiveFunctions)},
		{"multi_objective_weights",6,0,4,0,0,0.,0.,13,N_rem(RealDL,primaryRespFnWeights)},
		{"nonlinear_equality_constraints",0x29,6,6,0,kw_961,0.,0.,0,N_rem(sizet,numNonlinearEqConstraints)},
		{"nonlinear_inequality_constraints",0x29,8,5,0,kw_962,0.,0.,0,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_field_objectives",0x21,3,8,0,kw_960,0.,0.,-4,N_rem(sizet,numFieldObjectiveFunctions)},
		{"num_nonlinear_equality_constraints",0x21,6,6,0,kw_961,0.,0.,-3,N_rem(sizet,numNonlinearEqConstraints)},
		{"num_nonlinear_inequality_constraints",0x21,8,5,0,kw_962,0.,0.,-3,N_rem(sizet,numNonlinearIneqConstraints)},
		{"num_scalar_objectives",0x21,0,7,0,0,0.,0.,5,N_rem(sizet,numScalarObjectiveFunctions)},
		{"objective_function_scale_types",0x807,0,2,0,0,0.,0.,2,N_rem(strL,primaryRespFnScaleTypes)},
		{"objective_function_scales",0x806,0,3,0,0,0.,0.,2,N_rem(RealDL,primaryRespFnScales)},
		{"primary_scale_types",0x80f,0,2,0,0,0.,0.,0,N_rem(strL,primaryRespFnScaleTypes)},
		{"primary_scales",0x80e,0,3,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnScales)},
		{"scalar_objectives",0x29,0,7,0,0,0.,0.,0,N_rem(sizet,numScalarObjectiveFunctions)},
		{"sense",0x80f,0,1,0,0,0.,0.,0,N_rem(strL,primaryRespFnSense)},
		{"weights",14,0,4,0,0,0.,0.,0,N_rem(RealDL,primaryRespFnWeights)}
		},
	kw_964[3] = {
		{"lengths",13,0,1,1,0,0.,0.,0,N_rem(ivec,fieldLengths)},
		{"num_coordinates_per_field",13,0,2,0,0,0.,0.,0,N_rem(ivec,numCoordsPerField)},
		{"read_field_coordinates",8,0,3,0,0,0.,0.,0,N_rem(true,readFieldCoords)}
		},
	kw_965[4] = {
		{"field_responses",0x29,3,2,0,kw_964,0.,0.,0,N_rem(sizet,numFieldResponseFunctions)},
		{"num_field_responses",0x21,3,2,0,kw_964,0.,0.,-1,N_rem(sizet,numFieldResponseFunctions)},
		{"num_scalar_responses",0x21,0,1,0,0,0.,0.,1,N_rem(sizet,numScalarResponseFunctions)},
		{"scalar_responses",0x29,0,1,0,0,0.,0.,0,N_rem(sizet,numScalarResponseFunctions)}
		},
	kw_966[4] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_bounds)},
		{"ignore_bounds",8,0,1,0,0,0.,0.,0,N_rem(true,ignoreBounds)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdGradStepType_relative)}
		},
	kw_967[8] = {
		{"central",8,0,4,0,0,0.,0.,0,N_rem(lit,intervalType_central)},
		{"dakota",8,4,2,0,kw_966,0.,0.,0,N_rem(lit,methodSource_dakota)},
		{"fd_gradient_step_size",6,0,5,0,0,0.,0.,1,N_rem(RealL,fdGradStepSize)},
		{"fd_step_size",14,0,5,0,0,0.,0.,0,N_rem(RealL,fdGradStepSize)},
		{"forward",8,0,4,0,0,0.,0.,0,N_rem(lit,intervalType_forward)},
		{"interval_type",8,0,3},
		{"method_source",8,0,1},
		{"vendor",8,0,2,0,0,0.,0.,0,N_rem(lit,methodSource_vendor)}
		},
	kw_968[7] = {
		{"absolute",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_absolute)},
		{"bounds",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_bounds)},
		{"central",8,0,3,0,0,0.,0.,0,N_rem(true,centralHess)},
		{"fd_hessian_step_size",6,0,1,0,0,0.,0.,1,N_rem(RealL,fdHessStepSize)},
		{"fd_step_size",14,0,1,0,0,0.,0.,0,N_rem(RealL,fdHessStepSize)},
		{"forward",8,0,3,0,0,0.,0.,0,N_rem(false,centralHess)},
		{"relative",8,0,2,0,0,0.,0.,0,N_rem(lit,fdHessStepType_relative)}
		},
	kw_969[1] = {
		{"damped",8,0,1,0,0,0.,0.,0,N_rem(lit,quasiHessianType_damped_bfgs)}
		},
	kw_970[2] = {
		{"bfgs",8,1,1,1,kw_969,0.,0.,0,N_rem(lit,quasiHessianType_bfgs)},
		{"sr1",8,0,1,1,0,0.,0.,0,N_rem(lit,quasiHessianType_sr1)}
		},
	kw_971[19] = {
		{"analytic_gradients",8,0,4,2,0,0.,0.,0,N_rem(lit,gradientType_analytic)},
		{"analytic_hessians",8,0,5,3,0,0.,0.,0,N_rem(lit,hessianType_analytic)},
		{"calibration_terms",0x29,19,3,1,kw_953,0.,0.,0,N_rem(sizet,numLeastSqTerms)},
		{"descriptors",15,0,2,0,0,0.,0.,0,N_rem(strL,responseLabels)},
		{"id_responses",11,0,1,0,0,0.,0.,0,N_rem(str,idResponses)},
		{"least_squares_terms",0x21,19,3,1,kw_953,0.,0.,-3,N_rem(sizet,numLeastSqTerms)},
		{"mixed_gradients",8,10,4,2,kw_955,0.,0.,0,N_rem(lit,gradientType_mixed)},
		{"mixed_hessians",8,8,5,3,kw_959,0.,0.,0,N_rem(lit,hessianType_mixed)},
		{"no_gradients",8,0,4,2,0,0.,0.,0,N_rem(lit,gradientType_none)},
		{"no_hessians",8,0,5,3,0,0.,0.,0,N_rem(lit,hessianType_none)},
		{"num_least_squares_terms",0x21,19,3,1,kw_953,0.,0.,-8,N_rem(sizet,numLeastSqTerms)},
		{"num_objective_functions",0x21,15,3,1,kw_963,0.,0.,4,N_rem(sizet,numObjectiveFunctions)},
		{"num_response_functions",0x21,4,3,1,kw_965,0.,0.,6,N_rem(sizet,numResponseFunctions)},
		{"numerical_gradients",8,8,4,2,kw_967,0.,0.,0,N_rem(lit,gradientType_numerical)},
		{"numerical_hessians",8,7,5,3,kw_968,0.,0.,0,N_rem(lit,hessianType_numerical)},
		{"objective_functions",0x29,15,3,1,kw_963,0.,0.,0,N_rem(sizet,numObjectiveFunctions)},
		{"quasi_hessians",8,2,5,3,kw_970,0.,0.,0,N_rem(lit,hessianType_quasi)},
		{"response_descriptors",7,0,2,0,0,0.,0.,-14,N_rem(strL,responseLabels)},
		{"response_functions",0x29,4,3,1,kw_965,0.,0.,0,N_rem(sizet,numResponseFunctions)}
		},
	kw_972[6] = {
		{"aleatory",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_ALEATORY_UNCERTAIN_VIEW)},
		{"all",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_ALL_VIEW)},
		{"design",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_DESIGN_VIEW)},
		{"epistemic",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_EPISTEMIC_UNCERTAIN_VIEW)},
		{"state",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_STATE_VIEW)},
		{"uncertain",8,0,1,1,0,0.,0.,0,N_vam(type,varsView_UNCERTAIN_VIEW)}
		},
	kw_973[11] = {
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
	kw_974[5] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_binomial)},
		{"initial_point",13,0,3,0,0,0.,0.,0,N_vam(IntLb,binomialUncVars)},
		{"num_trials",13,0,2,2,0,0.,0.,0,N_vam(IntLb,binomialUncNumTrials)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,binomialUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,binomialUncProbPerTrial)}
		},
	kw_975[12] = {
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
	kw_976[10] = {
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
	kw_977[8] = {
		{"csv_descriptors",7,0,4,0,0,0.,0.,4,N_vam(strL,continuousStateLabels)},
		{"csv_initial_state",6,0,1,0,0,0.,0.,4,N_vam(rvec,continuousStateVars)},
		{"csv_lower_bounds",6,0,2,0,0,0.,0.,4,N_vam(rvec,continuousStateLowerBnds)},
		{"csv_upper_bounds",6,0,3,0,0,0.,0.,4,N_vam(rvec,continuousStateUpperBnds)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,continuousStateLabels)},
		{"initial_state",14,0,1,0,0,0.,0.,0,N_vam(rvec,continuousStateVars)},
		{"lower_bounds",14,0,2,0,0,0.,0.,0,N_vam(rvec,continuousStateLowerBnds)},
		{"upper_bounds",14,0,3,0,0,0.,0.,0,N_vam(rvec,continuousStateUpperBnds)}
		},
	kw_978[8] = {
		{"ddv_descriptors",7,0,4,0,0,0.,0.,4,N_vam(strL,discreteDesignRangeLabels)},
		{"ddv_initial_point",5,0,1,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeVars)},
		{"ddv_lower_bounds",5,0,2,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeLowerBnds)},
		{"ddv_upper_bounds",5,0,3,0,0,0.,0.,4,N_vam(ivec,discreteDesignRangeUpperBnds)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteDesignRangeLabels)},
		{"initial_point",13,0,1,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeVars)},
		{"lower_bounds",13,0,2,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeLowerBnds)},
		{"upper_bounds",13,0,3,0,0,0.,0.,0,N_vam(ivec,discreteDesignRangeUpperBnds)}
		},
	kw_979[1] = {
		{"adjacency_matrix",13,0,1,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddsia)}
		},
	kw_980[7] = {
		{"categorical",15,1,3,0,kw_979,0.,0.,0,N_vam(categorical,discreteDesignSetIntCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetIntLabels)},
		{"elements",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_ddsi)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddsi)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(ivec,discreteDesignSetIntVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddsi)},
		{"set_values",5,0,2,1,0,0.,0.,-4,N_vam(newivec,Var_Info_ddsi)}
		},
	kw_981[1] = {
		{"adjacency_matrix",13,0,1,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddsra)}
		},
	kw_982[7] = {
		{"categorical",15,1,3,0,kw_981,0.,0.,0,N_vam(categorical,discreteDesignSetRealCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetRealLabels)},
		{"elements",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_ddsr)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddsr)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,discreteDesignSetRealVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddsr)},
		{"set_values",6,0,2,1,0,0.,0.,-4,N_vam(newrvec,Var_Info_ddsr)}
		},
	kw_983[7] = {
		{"adjacency_matrix",13,0,3,0,0,0.,0.,0,N_vam(newivec,Var_Info_ddssa)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteDesignSetStrLabels)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_ddss)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nddss)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteDesignSetStrVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nddss)},
		{"set_values",7,0,2,1,0,0.,0.,-4,N_vam(newsarray,Var_Info_ddss)}
		},
	kw_984[3] = {
		{"integer",0x19,7,1,0,kw_980,0.,0.,0,N_vam(pintz,numDiscreteDesSetIntVars)},
		{"real",0x19,7,3,0,kw_982,0.,0.,0,N_vam(pintz,numDiscreteDesSetRealVars)},
		{"string",0x19,7,2,0,kw_983,0.,0.,0,N_vam(pintz,numDiscreteDesSetStrVars)}
		},
	kw_985[9] = {
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
	kw_986[8] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteStateRangeLabels)},
		{"dsv_descriptors",7,0,4,0,0,0.,0.,-1,N_vam(strL,discreteStateRangeLabels)},
		{"dsv_initial_state",5,0,1,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeVars)},
		{"dsv_lower_bounds",5,0,2,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeLowerBnds)},
		{"dsv_upper_bounds",5,0,3,0,0,0.,0.,3,N_vam(ivec,discreteStateRangeUpperBnds)},
		{"initial_state",13,0,1,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeVars)},
		{"lower_bounds",13,0,2,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeLowerBnds)},
		{"upper_bounds",13,0,3,0,0,0.,0.,0,N_vam(ivec,discreteStateRangeUpperBnds)}
		},
	kw_987[7] = {
		{"categorical",15,0,3,0,0,0.,0.,0,N_vam(categorical,discreteStateSetIntCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteStateSetIntLabels)},
		{"elements",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_dssi)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndssi)},
		{"initial_state",13,0,4,0,0,0.,0.,0,N_vam(ivec,discreteStateSetIntVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndssi)},
		{"set_values",5,0,2,1,0,0.,0.,-4,N_vam(newivec,Var_Info_dssi)}
		},
	kw_988[7] = {
		{"categorical",15,0,3,0,0,0.,0.,0,N_vam(categorical,discreteStateSetRealCat)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vam(strL,discreteStateSetRealLabels)},
		{"elements",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_dssr)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndssr)},
		{"initial_state",14,0,4,0,0,0.,0.,0,N_vam(rvec,discreteStateSetRealVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndssr)},
		{"set_values",6,0,2,1,0,0.,0.,-4,N_vam(newrvec,Var_Info_dssr)}
		},
	kw_989[6] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteStateSetStrLabels)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_dsss)},
		{"elements_per_variable",0x80d,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_ndsss)},
		{"initial_state",15,0,3,0,0,0.,0.,0,N_vam(strL,discreteStateSetStrVars)},
		{"num_set_values",0x805,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_ndsss)},
		{"set_values",7,0,2,1,0,0.,0.,-4,N_vam(newsarray,Var_Info_dsss)}
		},
	kw_990[3] = {
		{"integer",0x19,7,1,0,kw_987,0.,0.,0,N_vam(pintz,numDiscreteStateSetIntVars)},
		{"real",0x19,7,3,0,kw_988,0.,0.,0,N_vam(pintz,numDiscreteStateSetRealVars)},
		{"string",0x19,6,2,0,kw_989,0.,0.,0,N_vam(pintz,numDiscreteStateSetStrVars)}
		},
	kw_991[9] = {
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
	kw_992[9] = {
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
	kw_993[8] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(deuslbl,DEUSVar_set_str)},
		{"elements",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_duss)},
		{"elements_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nduss)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,discreteUncSetStrVars)},
		{"num_set_values",5,0,1,0,0,0.,0.,-2,N_vam(newiarray,Var_Info_nduss)},
		{"set_probabilities",14,0,3,0,0,0.,0.,0,N_vam(newrvec,Var_Info_DSSp)},
		{"set_probs",6,0,3,0,0,0.,0.,-1,N_vam(newrvec,Var_Info_DSSp)},
		{"set_values",7,0,2,1,0,0.,0.,-6,N_vam(newsarray,Var_Info_duss)}
		},
	kw_994[3] = {
		{"integer",0x19,9,1,0,kw_991,0.,0.,0,N_vam(pintz,numDiscreteUncSetIntVars)},
		{"real",0x19,9,3,0,kw_992,0.,0.,0,N_vam(pintz,numDiscreteUncSetRealVars)},
		{"string",0x19,8,2,0,kw_993,0.,0.,0,N_vam(pintz,numDiscreteUncSetStrVars)}
		},
	kw_995[5] = {
		{"betas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,exponentialUncBetas)},
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(caulbl,CAUVar_exponential)},
		{"euv_betas",6,0,1,1,0,0.,0.,-2,N_vam(RealLb,exponentialUncBetas)},
		{"euv_descriptors",7,0,3,0,0,0.,0.,-2,N_vae(caulbl,CAUVar_exponential)},
		{"initial_point",14,0,2,0,0,0.,0.,0,N_vam(RealLb,exponentialUncVars)}
		},
	kw_996[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,frechetUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(rvec,frechetUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_frechet)},
		{"fuv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,frechetUncAlphas)},
		{"fuv_betas",6,0,2,2,0,0.,0.,-3,N_vam(rvec,frechetUncBetas)},
		{"fuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_frechet)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,frechetUncVars)}
		},
	kw_997[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,gammaUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(RealLb,gammaUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_gamma)},
		{"gauv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,gammaUncAlphas)},
		{"gauv_betas",6,0,2,2,0,0.,0.,-3,N_vam(RealLb,gammaUncBetas)},
		{"gauv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_gamma)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,gammaUncVars)}
		},
	kw_998[4] = {
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_geometric)},
		{"initial_point",13,0,2,0,0,0.,0.,0,N_vam(IntLb,geometricUncVars)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,geometricUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,geometricUncProbPerTrial)}
		},
	kw_999[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,gumbelUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(rvec,gumbelUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_gumbel)},
		{"guuv_alphas",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,gumbelUncAlphas)},
		{"guuv_betas",6,0,2,2,0,0.,0.,-3,N_vam(rvec,gumbelUncBetas)},
		{"guuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_gumbel)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,gumbelUncVars)}
		},
	kw_1000[11] = {
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
	kw_1001[6] = {
		{"abscissas",13,0,2,1,0,0.,0.,0,N_vam(newivec,Var_Info_hpia)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hpic)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_histogram_point_int)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(ivec,histogramPointIntUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhpip)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhpip)}
		},
	kw_1002[6] = {
		{"abscissas",14,0,2,1,0,0.,0.,0,N_vam(newrvec,Var_Info_hpra)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hprc)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(daurlbl,DAURVar_histogram_point_real)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(rvec,histogramPointRealUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhprp)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhprp)}
		},
	kw_1003[6] = {
		{"abscissas",15,0,2,1,0,0.,0.,0,N_vam(newsarray,Var_Info_hpsa)},
		{"counts",14,0,3,2,0,0.,0.,0,N_vam(newrvec,Var_Info_hpsc)},
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauslbl,DAUSVar_histogram_point_str)},
		{"initial_point",15,0,4,0,0,0.,0.,0,N_vam(strL,histogramPointStrUncVars)},
		{"num_pairs",5,0,1,0,0,0.,0.,1,N_vam(newiarray,Var_Info_nhpsp)},
		{"pairs_per_variable",13,0,1,0,0,0.,0.,0,N_vam(newiarray,Var_Info_nhpsp)}
		},
	kw_1004[3] = {
		{"integer",0x19,6,1,0,kw_1001,0.,0.,0,N_vam(pintz,numHistogramPtIntUncVars)},
		{"real",0x19,6,3,0,kw_1002,0.,0.,0,N_vam(pintz,numHistogramPtRealUncVars)},
		{"string",0x19,6,2,0,kw_1003,0.,0.,0,N_vam(pintz,numHistogramPtStrUncVars)}
		},
	kw_1005[5] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_hypergeometric)},
		{"initial_point",13,0,4,0,0,0.,0.,0,N_vam(IntLb,hyperGeomUncVars)},
		{"num_drawn",13,0,3,3,0,0.,0.,0,N_vam(IntLb,hyperGeomUncNumDrawn)},
		{"selected_population",13,0,2,2,0,0.,0.,0,N_vam(IntLb,hyperGeomUncSelectedPop)},
		{"total_population",13,0,1,1,0,0.,0.,0,N_vam(IntLb,hyperGeomUncTotalPop)}
		},
	kw_1006[2] = {
		{"lnuv_zetas",6,0,1,1,0,0.,0.,1,N_vam(RealLb,lognormalUncZetas)},
		{"zetas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncZetas)}
		},
	kw_1007[4] = {
		{"error_factors",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncErrFacts)},
		{"lnuv_error_factors",6,0,1,1,0,0.,0.,-1,N_vam(RealLb,lognormalUncErrFacts)},
		{"lnuv_std_deviations",6,0,1,1,0,0.,0.,1,N_vam(RealLb,lognormalUncStdDevs)},
		{"std_deviations",14,0,1,1,0,0.,0.,0,N_vam(RealLb,lognormalUncStdDevs)}
		},
	kw_1008[11] = {
		{"descriptors",15,0,5,0,0,0.,0.,0,N_vae(caulbl,CAUVar_lognormal)},
		{"initial_point",14,0,4,0,0,0.,0.,0,N_vam(RealLb,lognormalUncVars)},
		{"lambdas",14,2,1,1,kw_1006,0.,0.,0,N_vam(rvec,lognormalUncLambdas)},
		{"lnuv_descriptors",7,0,5,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_lognormal)},
		{"lnuv_lambdas",6,2,1,1,kw_1006,0.,0.,-2,N_vam(rvec,lognormalUncLambdas)},
		{"lnuv_lower_bounds",6,0,2,0,0,0.,0.,3,N_vam(RealLb,lognormalUncLowerBnds)},
		{"lnuv_means",6,4,1,1,kw_1007,0.,0.,3,N_vam(RealLb,lognormalUncMeans)},
		{"lnuv_upper_bounds",6,0,3,0,0,0.,0.,3,N_vam(RealUb,lognormalUncUpperBnds)},
		{"lower_bounds",14,0,2,0,0,0.,0.,0,N_vam(RealLb,lognormalUncLowerBnds)},
		{"means",14,4,1,1,kw_1007,0.,0.,0,N_vam(RealLb,lognormalUncMeans)},
		{"upper_bounds",14,0,3,0,0,0.,0.,0,N_vam(RealUb,lognormalUncUpperBnds)}
		},
	kw_1009[7] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_loguniform)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,loguniformUncVars)},
		{"lower_bounds",14,0,1,1,0,0.,0.,0,N_vam(RealLb,loguniformUncLowerBnds)},
		{"luuv_descriptors",7,0,4,0,0,0.,0.,-3,N_vae(caulbl,CAUVar_loguniform)},
		{"luuv_lower_bounds",6,0,1,1,0,0.,0.,-2,N_vam(RealLb,loguniformUncLowerBnds)},
		{"luuv_upper_bounds",6,0,2,2,0,0.,0.,1,N_vam(RealUb,loguniformUncUpperBnds)},
		{"upper_bounds",14,0,2,2,0,0.,0.,0,N_vam(RealUb,loguniformUncUpperBnds)}
		},
	kw_1010[5] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_negative_binomial)},
		{"initial_point",13,0,3,0,0,0.,0.,0,N_vam(IntLb,negBinomialUncVars)},
		{"num_trials",13,0,2,2,0,0.,0.,0,N_vam(IntLb,negBinomialUncNumTrials)},
		{"prob_per_trial",6,0,1,1,0,0.,0.,1,N_vam(rvec,negBinomialUncProbPerTrial)},
		{"probability_per_trial",14,0,1,1,0,0.,0.,0,N_vam(rvec,negBinomialUncProbPerTrial)}
		},
	kw_1011[11] = {
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
	kw_1012[3] = {
		{"descriptors",15,0,3,0,0,0.,0.,0,N_vae(dauilbl,DAUIVar_poisson)},
		{"initial_point",13,0,2,0,0,0.,0.,0,N_vam(IntLb,poissonUncVars)},
		{"lambdas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,poissonUncLambdas)}
		},
	kw_1013[9] = {
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
	kw_1014[7] = {
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_uniform)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(rvec,uniformUncVars)},
		{"lower_bounds",14,0,1,1,0,0.,0.,0,N_vam(RealLb,uniformUncLowerBnds)},
		{"upper_bounds",14,0,2,2,0,0.,0.,0,N_vam(RealUb,uniformUncUpperBnds)},
		{"uuv_descriptors",7,0,4,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_uniform)},
		{"uuv_lower_bounds",6,0,1,1,0,0.,0.,-3,N_vam(RealLb,uniformUncLowerBnds)},
		{"uuv_upper_bounds",6,0,2,2,0,0.,0.,-3,N_vam(RealUb,uniformUncUpperBnds)}
		},
	kw_1015[7] = {
		{"alphas",14,0,1,1,0,0.,0.,0,N_vam(RealLb,weibullUncAlphas)},
		{"betas",14,0,2,2,0,0.,0.,0,N_vam(RealLb,weibullUncBetas)},
		{"descriptors",15,0,4,0,0,0.,0.,0,N_vae(caulbl,CAUVar_weibull)},
		{"initial_point",14,0,3,0,0,0.,0.,0,N_vam(RealLb,weibullUncVars)},
		{"wuv_alphas",6,0,1,1,0,0.,0.,-4,N_vam(RealLb,weibullUncAlphas)},
		{"wuv_betas",6,0,2,2,0,0.,0.,-4,N_vam(RealLb,weibullUncBetas)},
		{"wuv_descriptors",7,0,4,0,0,0.,0.,-4,N_vae(caulbl,CAUVar_weibull)}
		},
	kw_1016[42] = {
		{"active",8,6,2,0,kw_972},
		{"beta_uncertain",0x19,11,13,0,kw_973,0.,0.,0,N_vam(pintz,numBetaUncVars)},
		{"binomial_uncertain",0x19,5,20,0,kw_974,0.,0.,0,N_vam(pintz,numBinomialUncVars)},
		{"continuous_design",0x19,12,4,0,kw_975,0.,0.,0,N_vam(pintz,numContinuousDesVars)},
		{"continuous_interval_uncertain",0x19,10,26,0,kw_976,0.,0.,0,N_vam(pintz,numContinuousIntervalUncVars)},
		{"continuous_state",0x19,8,29,0,kw_977,0.,0.,0,N_vam(pintz,numContinuousStateVars)},
		{"discrete_design_range",0x19,8,5,0,kw_978,0.,0.,0,N_vam(pintz,numDiscreteDesRangeVars)},
		{"discrete_design_set",8,3,6,0,kw_984},
		{"discrete_interval_uncertain",0x19,9,27,0,kw_985,0.,0.,0,N_vam(pintz,numDiscreteIntervalUncVars)},
		{"discrete_state_range",0x19,8,30,0,kw_986,0.,0.,0,N_vam(pintz,numDiscreteStateRangeVars)},
		{"discrete_state_set",8,3,31,0,kw_990},
		{"discrete_uncertain_set",8,3,28,0,kw_994},
		{"exponential_uncertain",0x19,5,12,0,kw_995,0.,0.,0,N_vam(pintz,numExponentialUncVars)},
		{"frechet_uncertain",0x19,7,16,0,kw_996,0.,0.,0,N_vam(pintz,numFrechetUncVars)},
		{"gamma_uncertain",0x19,7,14,0,kw_997,0.,0.,0,N_vam(pintz,numGammaUncVars)},
		{"geometric_uncertain",0x19,4,22,0,kw_998,0.,0.,0,N_vam(pintz,numGeometricUncVars)},
		{"gumbel_uncertain",0x19,7,15,0,kw_999,0.,0.,0,N_vam(pintz,numGumbelUncVars)},
		{"histogram_bin_uncertain",0x19,11,18,0,kw_1000,0.,0.,0,N_vam(pintz,numHistogramBinUncVars)},
		{"histogram_point_uncertain",8,3,24,0,kw_1004},
		{"hypergeometric_uncertain",0x19,5,23,0,kw_1005,0.,0.,0,N_vam(pintz,numHyperGeomUncVars)},
		{"id_variables",11,0,1,0,0,0.,0.,0,N_vam(str,idVariables)},
		{"interval_uncertain",0x11,10,26,0,kw_976,0.,0.,-17,N_vam(pintz,numContinuousIntervalUncVars)},
		{"linear_equality_constraint_matrix",14,0,37,0,0,0.,0.,0,N_vam(rvec,linearEqConstraintCoeffs)},
		{"linear_equality_scale_types",15,0,39,0,0,0.,0.,0,N_vam(strL,linearEqScaleTypes)},
		{"linear_equality_scales",14,0,40,0,0,0.,0.,0,N_vam(rvec,linearEqScales)},
		{"linear_equality_targets",14,0,38,0,0,0.,0.,0,N_vam(rvec,linearEqTargets)},
		{"linear_inequality_constraint_matrix",14,0,32,0,0,0.,0.,0,N_vam(rvec,linearIneqConstraintCoeffs)},
		{"linear_inequality_lower_bounds",14,0,33,0,0,0.,0.,0,N_vam(rvec,linearIneqLowerBnds)},
		{"linear_inequality_scale_types",15,0,35,0,0,0.,0.,0,N_vam(strL,linearIneqScaleTypes)},
		{"linear_inequality_scales",14,0,36,0,0,0.,0.,0,N_vam(rvec,linearIneqScales)},
		{"linear_inequality_upper_bounds",14,0,34,0,0,0.,0.,0,N_vam(rvec,linearIneqUpperBnds)},
		{"lognormal_uncertain",0x19,11,8,0,kw_1008,0.,0.,0,N_vam(pintz,numLognormalUncVars)},
		{"loguniform_uncertain",0x19,7,10,0,kw_1009,0.,0.,0,N_vam(pintz,numLoguniformUncVars)},
		{"mixed",8,0,3,0,0,0.,0.,0,N_vam(type,varsDomain_MIXED_DOMAIN)},
		{"negative_binomial_uncertain",0x19,5,21,0,kw_1010,0.,0.,0,N_vam(pintz,numNegBinomialUncVars)},
		{"normal_uncertain",0x19,11,7,0,kw_1011,0.,0.,0,N_vam(pintz,numNormalUncVars)},
		{"poisson_uncertain",0x19,3,19,0,kw_1012,0.,0.,0,N_vam(pintz,numPoissonUncVars)},
		{"relaxed",8,0,3,0,0,0.,0.,0,N_vam(type,varsDomain_RELAXED_DOMAIN)},
		{"triangular_uncertain",0x19,9,11,0,kw_1013,0.,0.,0,N_vam(pintz,numTriangularUncVars)},
		{"uncertain_correlation_matrix",14,0,25,0,0,0.,0.,0,N_vam(newrvec,Var_Info_ucm)},
		{"uniform_uncertain",0x19,7,9,0,kw_1014,0.,0.,0,N_vam(pintz,numUniformUncVars)},
		{"weibull_uncertain",0x19,7,17,0,kw_1015,0.,0.,0,N_vam(pintz,numWeibullUncVars)}
		},
	kw_1017[6] = {
		{"environment",0x108,15,1,1,kw_15,0.,0.,0,NIDRProblemDescDB::env_start},
		{"interface",0x308,12,5,5,kw_32,0.,0.,0,N_ifm3(start,0,stop)},
		{"method",0x308,92,2,2,kw_875,0.,0.,0,N_mdm3(start,0,stop)},
		{"model",8,12,3,3,kw_944,0.,0.,0,N_mom3(start,0,stop)},
		{"responses",0x308,19,6,6,kw_971,0.,0.,0,N_rem3(start,0,stop)},
		{"variables",0x308,42,4,4,kw_1016,0.,0.,0,N_vam3(start,0,stop)}
		};

} // namespace Dakota

#ifdef __cplusplus
extern "C" {
#endif
KeyWord Dakota_Keyword_Top = {"KeywordTop",0,6,0,0,Dakota::kw_1017};
#ifdef __cplusplus
}
#endif
