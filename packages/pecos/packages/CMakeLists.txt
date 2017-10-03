include(ManageRelocatablePackage)
include(FftwExternalProject)

if(HAVE_FFT)
  ManageRelocatablePackage(DFFTPACK dfftpack)
  FftwExternalProject()
endif(HAVE_FFT)

ManageRelocatablePackage(LHS LHS)

ManageRelocatablePackage(SPARSE_GRID VPISparseGrid)
