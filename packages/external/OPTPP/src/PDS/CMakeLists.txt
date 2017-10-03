set(PDS_sources
  create_scheme.C dqrdc.c make_search.c pds.C pdschk.C pdscld.c pdscom.c
  pdsdgn.c pdsdone.c pdseql.c pdsget.c pdshrk.c pdsinit.C pdslen.c pdsopt.C
  pdsrgt.c pdsstep.C pdsupd.c pdswork.C sort.c writes.c
)

if(OPTPP_HAVE_MPI AND MPI_FOUND)
  list(APPEND PDS_sources pdsglb.c pdsgop.c pdswap.c) 
endif()
