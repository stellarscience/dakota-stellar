
DAKOTA_PROJ=/home/rhoope/Projects/Dakota
MUQ_SRC_DIR=/home/rhoope/temp/muq2

#module load sems-boost/1.63.0
#module load sems-hdf5/1.10.4

rm -rf CMakeCache.txt CMakeFiles

cmake \
  -D CMAKE_INSTALL_PREFIX=${DAKOTA_PROJ}/packages/external/muq2 \
  -D CMAKE_CXX_FLAGS:STRING="-L/projects/sems/install/rhel6-x86_64/sems/tpl/boost/1.63.0/gcc/4.8.4/base/lib" \
  -D MUQ_EIGEN3_DIR:FILEPATH=${DAKOTA_PROJ}/packages/external/eigen3/include/eigen3 \
  -D MUQ_BOOST_DIR:FILEPATH=/projects/sems/install/rhel6-x86_64/sems/tpl/boost/1.63.0/gcc/4.8.4/base \
  -D MUQ_HDF5_DIR:FILEPATH=/projects/sems/install/rhel6-x86_64/sems/tpl/hdf5/1.10.4/gcc/4.8.4/base \
  -D MUQ_ENABLEGROUP_MODELING_SUNDIALS_MODELS:BOOL=OFF \
  -D MUQ_USE_NLOPT:BOOL=OFF \
  -D MUQ_ENABLEGROUP_DEFAULT:BOOL=OFF \
  -D MUQ_ENABLEGROUP_SAMPLING_ALGORITHM:BOOL=ON \
  -D MUQ_ENABLEGROUP_UTILITIES_MULTIINDEX:BOOL=ON \
  -D MUQ_ENABLEGROUP_MODELING_FLANN:BOOL=OFF \
  ${MUQ_SRC_DIR}

# Make and install MUQ in the packages/external directory.
#https_proxy=http://wwwproxy.sandia.gov:80 http_proxy=http://wwwproxy.sandia.gov:80 make -j15 install
