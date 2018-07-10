#!/usr/bin/env bash
# Aim: distance matrix for all PDBs with feature vectors
cd feature_vectors
for file in *.txt
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/feature_distances.py
done
cd ..