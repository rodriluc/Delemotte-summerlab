#!/usr/bin/env bash

#Overall aim: to save feature vectors

# Strip PDB files in correct format
cd PDB_files
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/strip_pdb.py
done
cd ..
# Create individual PDB files for hydrophobic and charged residues
cd PDB_edited
for file in *pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/create_pdb.py 
done
cd ..
# Compute side-chain distance map from all PDBs now
cd PDB_edited
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -sc_cmap
done
# Get binary contact map Origin pdb
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -d_in /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -bin -coff 0.45
done
cd ..
cd hydrophobic_files
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -sc_cmap
done
# Get binary contact map Hydrophobic pdb
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -d_in /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -bin -coff 0.45
done
cd ..
cd charged_files
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -sc_cmap
done
# Get binary contact map Charged pdb
for file in *.pdb
do
echo $file
python /afs/kth.se/home/l/u/lucier/Documents/protein_networks/MD_cmaps.py -d_in /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od /afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/ -bin -coff 0.45
done
cd ..
# Create feature vector
python input_netx.py
