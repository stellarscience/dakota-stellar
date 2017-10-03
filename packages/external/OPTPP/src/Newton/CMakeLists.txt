set(Newton_sources
  OptBaNewton.C OptBaQNewton.C OptBCEllipsoid.C OptBCFDNewton.C OptBCNewton.C
  OptBCNewtonLike.C OptBCQNewton.C OptCG.C OptConstrFDNewton.C OptConstrNewton.C
  OptConstrNewtonLike.C OptConstrQNewton.C OptDHNIPS.C OptFDNewton.C OptFDNIPS.C
  OptLBFGS.C OptNewton.C OptNewtonLike.C OptNIPS.C OptNIPSLike.C OptQNewton.C
  OptQNIPS.C
)

if(npsol_lib)
  list(APPEND Newton_sources OptNPSOL.C npsol_setup.c)
endif()
