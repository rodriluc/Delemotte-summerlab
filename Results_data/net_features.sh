#!/usr/bin/env bash

# Path to input_netx.py
#pypath= 'C:\Users\lucie\AppData\Local\Temp\scp31535\afs\kth.se\home\l\u\lucier\Documents\protein_networks\Results_data\input_netx.py'

# Save deature vectors
for file in *.pdb
do
echo $file
python input_netx.py -top $file -trj file -fe ${file%.*} -od Results_data/ -sc_cmap
done
