#!/usr/bin/env bash
#pypath = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_edited'
#${pypath}

#Overall aim: to save feature vectors

# Strip PDB files in correct format
cd PDB_files
for file in *pdb
do
echo $file
python strip_pdb.py
done
cd ..
# Create individual PDB files for hydrophobic and charged residues
cd PDB_edited
for file in *pdb
do
echo $file
python create_pdb.py 
done
cd ..
# Compute side-chain distance map from all PDBs now
cd PDB_edited
for file in *.pdb
do
echo $file
python MD_cmaps.py -top $file -trj $file -fe ${file%.*} -od Results_data/ -sc_cmap
done
# Get binary contact map 
for file in *.pdb
do
echo $file
python ${pypath}MD_cmaps.py -d_in Results_data/distance_matrix_min_${file%.*}.txt -fe ${file%.*} -od Results_data/ -bin -coff 0.45
done
cd ..
python input_netx.py #create feature vector