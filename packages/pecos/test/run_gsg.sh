#!/bin/bash

#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

if [ -z "$PECOS_INS" ]; then
    echo "Need to set PECOS_INS"
    exit 1
fi  

# runs pecos gsg driver with default params
$PECOS_INS/bin/pecos_gsg_driver >& gsg_default.log

# runs pecos gsg driver with set 1
$PECOS_INS/bin/pecos_gsg_driver -d 2 -i 0 >& gsg_set1.log


