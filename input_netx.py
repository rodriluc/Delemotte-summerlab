#!/usr/bin/env python
import argparse
import networkx as nx
import numpy as np
from numpy import genfromtxt
import os
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
#PATH=$PATH :/home/l/u/lucier/.local/bin/ #for bash jupyter notebook need to link path
import scipy
import scipy.sparse
import mdtraj as md
import itertools
import pandas as pd
import re

python_path = os.path.dirname(__file__);

next_folder = '';
parent_folder = '';
for i in range(len(python_path)-1):
	next_folder+=python_path[i];
	if python_path[i]=='/':
		parent_folder += next_folder;
		next_folder = '';
		
sys.path.append(python_path);
sys.path.append(parent_folder);

import MD_cmaps

path_data = '/home/lrodriguez/Documents/Delemotte-summerlab/Results_data/'
path_pdb = '/home/lrodriguez/Documents/Delemotte-summerlab/PDB_edited/'
feat_path = '/home/lrodriguez/Documents/Delemotte-summerlab/feature_vectors/'
path_hydro = '/home/lrodriguez/Documents/Delemotte-summerlab/hydrophobic_files/'
path_charg = '/home/lrodriguez/Documents/Delemotte-summerlab/charged_files/'
path_cmap = '/home/lrodriguez/Documents/Delemotte-summerlab/network_files/'

os.path.splitext(path_pdb)[0]
from os.path import basename

#class input_netx():

file_end_name = ''
save_folder = ''
inpnetx = []

def install(package):
	subprocess.call([sys.executable, "-m", "pip", "install", package])
	
def __init__():
	return

def name_base():
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			#base_list.append(base)
			return base

def input_nx(): 
	base = name_base()
	origin_list = []
	print os.listdir(path_data)
	for file in os.listdir(path_data):
		print file
		#for i in range(len(base)):
		origin2_list = []
		if file.startswith('cmap_processed_') and not file.startswith('cmap_processed_hydrophobic') and not file.startswith('cmap_processed_charged'): 
			#print file
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) #usecols=(0, -1) blank line at end
			B = np.matrix(np.array(A))
			G = nx.from_numpy_matrix(B)
			origin2_list.append(G)

			nx.draw(G)

			plt.show()
			plt.savefig(path_cmap+'cmap_'+base+'.svg') #'cmap.svg'
			print 'Figure has been saved as *.svg' #called twice so prints twice
			#print origin2_list
			return origin2_list
	
def alpha_content():
	base = name_base()
	alpha_list = []
	for file in os.listdir(path_pdb):
		#print file
		#if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'): 
			
		iter_file = open((path_pdb+file), 'r')
		lines = iter_file.readlines()
		cryst = ('CRYST1')
		for line in lines:
			col = line.split()
			if cryst in line:
				alpha_list.append(col[1])
	return alpha_list
			
def beta_content():
	base = name_base()
	beta_list = []
	for file in os.listdir(path_pdb):
		#print file
		#if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'): 
			
		iter_file = open((path_pdb+file), 'r')
		lines = iter_file.readlines()
		cryst = ('CRYST1')
		for line in lines:
			col = line.split()
			if cryst in line:
				beta_list.append(col[2])
	return beta_list

def load_hydrophobic():
	base = name_base()
	hydrophobic_list = []
	for file in os.listdir(path_data):
		#for i in range(len(base)):
		if file.startswith('cmap_processed_hydrophobic'): #'cmap_processed_hydrophobic'
			#print file
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) 
			B = np.matrix(np.array(A))
			H = nx.from_numpy_matrix(B)
			nx.draw(H)
			plt.show()
			plt.savefig(path_cmap+'cmap_hydrophobic_'+base+'.svg')
			hydrophobic_list.append(H)
	return hydrophobic_list
				
def load_charged():
	base = name_base()
	charged_list = []
	for file in os.listdir(path_data):
		#for i in range(len(base)):
		if file.startswith('cmap_processed_charged'): #'cmap_processed_charged'
			#print file
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) 
			B = np.matrix(np.array(A))
			C = nx.from_numpy_matrix(B)
			nx.draw(C)
			plt.show()
			plt.savefig(path_cmap+'cmap_charged_'+base+'.svg')
			charged_list.append(C)
	return charged_list
			
def laplacian_walk():
	G = input_nx()
	lap_list = []

	for file in os.listdir(path_pdb):
		for i in range(len(G)):
			#print i
			lap = nx.normalized_laplacian_matrix(G[i], weight='weight')
			lap = np.eye(lap.shape[0])-lap #random walk: identity matrix - lap
			eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
			lap_list.append(eigenvalues[1]-eigenvalues[0]) #difference 1-second largest
	return lap_list
	
def features_vector(): #self
	base = name_base()
	G = input_nx()
	L = laplacian_walk()
	H = load_hydrophobic()
	C = load_charged()
	A = alpha_content()
	B = beta_content()

	for file in os.listdir(path_pdb):
		print file
		#for i, base in enumerate(base):
		#for i in range(len(base)):
			#print base[i]
		with open ((feat_path+('features_'+base+'.txt')), 'w') as w: #base[i], works well with file just additional .pdb ending
			for i in range(len(G)):
				print len(G)
				print(str(float(nx.number_of_edges(G[i]))/float(nx.number_of_nodes(G[i])))+'\n') #avg degree
				w.write(str(nx.average_shortest_path_length(G[i]))+'\n') #avg shortest path length
				w.write(str(nx.diameter(G[i]))+'\n') #diameter (max. shortest path)
				w.write(str(nx.radius(G[i]))+'\n') #radius (min. shortest path)
				w.write(str(nx.average_clustering(G[i]))+'\n') #global clustering coefficient
				#***Number of quasi-rigid domains
				w.write(str(nx.degree_assortativity_coefficient(G[i]))+'\n') #assortativity coefficient
				w.write(str(nx.global_reaching_centrality(G[i]))+'\n') #nx.degree_centrality(G)
				#***Residue intrinsic dimensionality (may be used to compute 6.?)
			for i in range(len(L)):
				#print len(L)
				w.write(str(L[i])+'\n') #Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
			for i in range(len(A)):
				w.write(str(A[i])+'\n') #alpha content 
			for i in range(len(B)):
				w.write(str(B[i])+'\n') #beta content
			
			for i in range(len(H)):
				w.write(str(float(nx.number_of_edges(H[i]))/float(nx.number_of_nodes(H[i])))+'\n') #***Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
				w.write(str(nx.average_clustering(H[i]))+'\n') #***Average local clustering coefficient of hydrophobic residues
				w.write(str(nx.global_reaching_centrality(H[i]))+'\n') #global, Average local reaching centrality of hydrophobic residues
			
			for i in range(len(C)):
				w.write(str(float(nx.number_of_edges(C[i]))/float(nx.number_of_nodes(C[i])))+'\n') #***Average degree of charged residues (R,D,E,H,K)
				w.write(str(nx.average_clustering(C[i]))+'\n') #***Average local clustering coefficient of charged residues
				w.write(str(nx.global_reaching_centrality(C[i]))+'\n') #global, Average local reaching centrality of charged residues
	
	print 'Features vector *.txt file has been saved! (original PDB, hydrophobic, charged)'

def main(parser):
	parser.add_argument('-pdb','--pdb_file',help='Input 1 PDB file (.pdb)',type=str,default='')
	parser.add_argument('-fe','--file_end_name',type=str,help='Output file end name (optional)', default='')
	parser.add_argument('-od','--out_directory',type=str,help='The directory where data should be saved (optional)',default='')
	parser.add_argument('fv','--feature_vector',type=str,help='Feature vector function')
	
	args = parser.parse_args()
	
	self.save_folder = args.out_directory;
	self.file_end_name = args.file_end_name;
	
	#for file in os.listdir(path_data):
		#if file.endswith('.pdb'):
			#self.traj = md.load_pdb(path_data+file) #args.topology_file
	if args.feature_vector:
			self.features_vector(self)

if __name__ == '__main__':
	print(input_nx()) 
	#name_base()
	#alpha_content()
	#beta_content()
	#load_hydrophobic()
	#load_charged()
	#laplacian_walk()
	#features_vector()
	
	#inpnetx_obj = input_netx()
	#inpnetx_obj.main(parser)

