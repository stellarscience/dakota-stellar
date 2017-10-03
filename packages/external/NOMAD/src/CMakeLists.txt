## NOMAD library sources.
set(nomad_SOURCES
  Barrier.cpp
  Cache.cpp
  Cache_File_Point.cpp
  Cache_Point.cpp
  Cache_Search.cpp
  Clock.cpp
  Direction.cpp
  Directions.cpp
  Display.cpp
  Double.cpp
  Eval_Point.cpp
  Evaluator_Control.cpp
  Evaluator.cpp
  Exception.cpp
  Extended_Poll.cpp
  L_Curve.cpp
  LH_Search.cpp
  Mads.cpp
  Model_Sorted_Point.cpp
  Model_Stats.cpp
  Multi_Obj_Evaluator.cpp
  nomad.cpp
  OrthogonalMesh.cpp
  Parameter_Entries.cpp
  Parameter_Entry.cpp
  Parameters.cpp
  Pareto_Front.cpp
  Pareto_Point.cpp
  Phase_One_Evaluator.cpp
  Phase_One_Search.cpp
  Point.cpp
  Priority_Eval_Point.cpp
  Quad_Model.cpp
  Quad_Model_Evaluator.cpp
  Quad_Model_Search.cpp
  Random_Pickup.cpp
  RNG.cpp
  Signature.cpp
  Slave.cpp
  SMesh.cpp
  Speculative_Search.cpp
  Stats.cpp
  TGP_Model.cpp
  TGP_Model_Evaluator.cpp
  TGP_Model_Search.cpp
  TGP_Output_Model.cpp
  utils.cpp
  Variable_Group.cpp
  VNS_Search.cpp
  XMesh.cpp
)

## Set the NOMAD library name.
add_library(nomad ${nomad_SOURCES})

install(TARGETS nomad EXPORT ${ExportTarget} DESTINATION lib)
