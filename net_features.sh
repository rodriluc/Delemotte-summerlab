#!/usr/bin/env bash
#Overall aim: to save feature vectors
python create_pdb.py #create PDB files first
# Compute side-chain distance map from all PDBs now
for file in *.pdb
do
echo $file
python MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od Results_data/ -sc_cmap
done
# Get binary contact map 
for file in *.pdb
do
echo $file
python MD_cmaps.py -d_in Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od Results_data/ -bin -coff 0.45
done
python input_netx.py #create feature vector