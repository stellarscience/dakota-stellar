#!/bin/sh 

#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# Link to dakota-packages when using Pecos separately from Dakota

dakota_pkgs_dir=dakota-packages
if [ $# -gt 1 ]; then
  echo "Usage ${0} [ /path/to/dakota/packages/external ]"
  exit 1
elif [ $# -eq 1 ]; then
  dakota_pkgs_dir=${1}
fi

if [ -e ${dakota_pkgs_dir} ]; then
  echo "Not cloning dakota-packages; using existing ${dakota_pkgs_dir}."
else
  echo "Cloning dakota-packages into ${dakota_pkgs_dir}."
  git clone software.sandia.gov:/git/dakota-packages ${dakota_pkgs_dir}
fi

for pkg in dfftpack fftw LHS trilinos VPISparseGrid; do
  if [ -e ${pkg} ]; then
    echo "${pkg} exists; skipping."
  else
    if [ -d ${dakota_pkgs_dir}/${pkg} ]; then
      ln -v -s ${dakota_pkgs_dir}/${pkg} ${pkg}
    else
      echo "Link target ${dakota_pkgs_dir}/${pkg} does not exist; skipping."
    fi
  fi
done
