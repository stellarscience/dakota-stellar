## CONMIN library sources.
set(conmin_SOURCES
  conmin.f		
)

## Set the CONMIN library name.
add_library(conmin ${conmin_SOURCES})

INSTALL(TARGETS conmin EXPORT ${ExportTarget} DESTINATION lib)

# WJB - NOTE: 'sp_' prefix is used to distinguish this CONMIN lib
# from the one the built as part of DAKOTA.  Work with Zack for an 
# improved design.
#if(HAVE_SP_CONMIN)
#  ## Set the SURFPACK version of CONMIN library name.
#  add_library(sp_conmin SHARED ${conmin_SOURCES})
#  add_library(sp_conmin-static STATIC ${conmin_SOURCES})
#
#  ## Ensure both the shared and static libraries have the same root name
#  set_target_properties(sp_conmin-static PROPERTIES OUTPUT_NAME "sp_conmin")
#  set_target_properties(sp_conmin-static PROPERTIES PREFIX "lib")
#endif(HAVE_SP_CONMIN)


#set_target_properties(sp_conmin PROPERTIES LINKER_LANGUAGE Fortran)

