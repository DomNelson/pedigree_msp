set -e

cd ~/project/msprime/lib/build

ninja
valgrind ./dev-cli simulate <(cat ~/project/pedigree_msp/data/example_for_ped.cfg ~/project/pedigree_msp/data/ll_msp_sample_config.txt ~/project/pedigree_msp/data/ll_pedEx.txt) -o ~/temp/trees.h5

