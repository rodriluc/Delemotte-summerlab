#!/usr/bin/env bash

# Compute side-chain distance map, removed -trj 1CLL_pr.dcd
for file in *.pdb
do
echo ${file%.*}
python MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od Results_data/ -sc_cmap
done
# Get binary contact map
for file in *.pdb
do
echo ${file%.*}
python MD_cmaps.py -d_in Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od Results_data/ -bin -coff 0.45
done
