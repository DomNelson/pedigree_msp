#!/bin/bash

# Get script directory regardless of where it was called from
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd)"

cd ~/project/msprime/lib/build

ninja && valgrind ./dev-cli simulate <(cat $SCRIPTDIR/test/example_for_ped.cfg $SCRIPTDIR/test/ll_msp_sample_config.txt $SCRIPTDIR/test/ll_pedEx.txt) -o ~/temp/trees.h5
