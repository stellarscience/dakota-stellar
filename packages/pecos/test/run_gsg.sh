#!/bin/bash

if [ -z "$PECOS_INS" ]; then
    echo "Need to set PECOS_INS"
    exit 1
fi  

# runs pecos gsg driver with default params
$PECOS_INS/bin/pecos_gsg_driver >& gsg_default.log

# runs pecos gsg driver with set 1
$PECOS_INS/bin/pecos_gsg_driver -d 2 -i 0 >& gsg_set1.log


