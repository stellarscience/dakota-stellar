#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

include(ManageRelocatablePackage)
#include(FftwExternalProject)

if(HAVE_FFT)
  ManageRelocatablePackage(DFFTPACK dfftpack)
  #FftwExternalProject()
endif(HAVE_FFT)

ManageRelocatablePackage(LHS LHS)

ManageRelocatablePackage(SPARSE_GRID VPISparseGrid)
