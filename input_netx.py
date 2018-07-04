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



def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			#f = os.path.basename(file_nm)
			basename = file.split('.')[:-1]
			base =''.join(basename)
		#return base
		base_list.append(base)
	return base_list

def input_nx(path_data,file): 
	base_list = name_base(path_pdb,file)
	origin_list = []
	#print os.listdir(path_data)
	#for file in os.listdir(path_data):
		#print file
		#for i in range(len(base)):

	#if file.startswith('cmap_processed_') and not file.startswith('cmap_processed_hydrophobic') and not file.startswith('cmap_processed_charged'): 
	#print file
	#if os.path.isfile(os.path.join(path_data,i)) and i.startswith('cmap_processed_') and not i.startswith('cmap_processed_hydrophobic') and not i.startswith('cmap_processed_charged') in i:
	for i in range(len(base_list)):
		A = np.loadtxt(path_data+'cmap_processed_'+base_list[i]+'.txt', dtype=float, unpack=True) 
		B = np.matrix(np.array(A))
		G = nx.from_numpy_matrix(B)
		origin_list.append(G)

		nx.draw(G)

		plt.show()
		plt.savefig(path_cmap+'cmap_'+base_list[i]+'.svg') #'cmap.svg'
		print 'Figure has been saved as *.svg' #called twice so prints twice
		#print origin_list
		#print len(G)
		return G

def alpha_content(path_pdb,file):
	base = name_base(path_pdb,file)
	alpha_list = []
	for file in os.listdir(path_pdb):
	#iter_file = open((path_pdb+file), 'r')
		with open(path_pdb+file) as iter_file:
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines:
				col = line.split()
				if cryst in line:
					return col[1]
			#alpha_list.append(col[1])
	#return alpha_list
			
def beta_content(path_pdb,file):
	base = name_base(path_pdb,file)
	beta_list = []
	for file in os.listdir(path_pdb):	
	#iter_file = open((path_pdb+file), 'r')
		with open(path_pdb+file) as iter_file:
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines:
				col = line.split()
				if cryst in line:
					return col[2]
			#beta_list.append(col[2])
	#return beta_list

def load_hydrophobic(path_data,file):
	base_list = name_base(path_pdb,file)
	hydrophobic_list = []
	#for file in os.listdir(path_data):
		#for i in range(len(base)):
	#if file.startswith('cmap_processed_hydrophobic'): #'cmap_processed_hydrophobic'
		#print file\
	for i in range(len(base_list)):
		A = np.loadtxt(path_data+'cmap_processed_hydrophobic_'+base_list[i]+'.txt', dtype=float, unpack=True) 
		B = np.matrix(np.array(A))
		H = nx.from_numpy_matrix(B)
		nx.draw(H)
		plt.show()
		plt.savefig(path_cmap+'cmap_hydrophobic_'+base_list[i]+'.svg')
		#hydrophobic_list.append(H)
		return H
				
def load_charged(path_data,file):
	base_list = name_base(path_pdb,file)
	charged_list = []
	#for file in os.listdir(path_data):
		#for i in range(len(base)):
	#if file.startswith('cmap_processed_charged'): #'cmap_processed_charged'
		#print file
	for i in range(len(base_list)):
		A = np.loadtxt(path_data+'cmap_processed_charged_'+base_list[i]+'.txt', dtype=float, unpack=True) 
		B = np.matrix(np.array(A))
		C = nx.from_numpy_matrix(B)
		nx.draw(C)
		plt.show()
		plt.savefig(path_cmap+'cmap_charged_'+base_list[i]+'.svg')
		#charged_list.append(C)
		return C
			
def laplacian_walk(G):
	G = input_nx(path_data,file)
	lap_list = []

	#for file in os.listdir(path_pdb):
	#for i in range(len(G)):
		#print i
	lap = nx.normalized_laplacian_matrix(G, weight='weight')
	lap = np.eye(lap.shape[0])-lap #random walk: identity matrix - lap
	eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
	lap_sum = (eigenvalues[1]-eigenvalues[0]) #difference 1-second largest
	return lap_sum
	
def features_vector(file): #self

	#for file in os.listdir(path_pdb):
	base_list = name_base(path_pdb,file)
	
	#print file
	#for i, base in enumerate(base):
	#for i in range(len(base)):
			#print base[i]
	for i in range(len(base_list)):
		G = input_nx(path_data,file)
		L = laplacian_walk(G)
		H = load_hydrophobic(path_data,file)
		C = load_charged(path_data,file)
		A = alpha_content(path_pdb,file)
		B = beta_content(path_pdb,file)
		#print len(base_list)
		with open ((feat_path+('features_'+base_list[i]+'.txt')), 'w') as w: 
			#for i in range(len(G)):
			w.write(str(float(nx.number_of_edges(G))/float(nx.number_of_nodes(G)))+'\n') #avg degree
			w.write(str(nx.average_shortest_path_length(G))+'\n') #avg shortest path length
			w.write(str(nx.diameter(G))+'\n') #diameter (max. shortest path)
			w.write(str(nx.radius(G))+'\n') #radius (min. shortest path)
			w.write(str(nx.average_clustering(G))+'\n') #global clustering coefficient
			#***Number of quasi-rigid domains
			w.write(str(nx.degree_assortativity_coefficient(G))+'\n') #assortativity coefficient
			w.write(str(nx.global_reaching_centrality(G))+'\n') #nx.degree_centrality(G)
			#***Residue intrinsic dimensionality (may be used to compute 6.?)
			#for i in range(len(L)):
			#print len(L)
			w.write(str(L)+'\n') #Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
			#for i in range(len(A)):
			w.write(str(A)+'\n') #alpha content 
			#for i in range(len(B)):
			w.write(str(B)+'\n') #beta content
	
			#for i in range(len(H)):
			w.write(str(float(nx.number_of_edges(H))/float(nx.number_of_nodes(H)))+'\n') #***Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
			w.write(str(nx.average_clustering(H))+'\n') #***Average local clustering coefficient of hydrophobic residues
			w.write(str(nx.global_reaching_centrality(H))+'\n') #global, Average local reaching centrality of hydrophobic residues
	
			#for i in range(len(C)):
			w.write(str(float(nx.number_of_edges(C))/float(nx.number_of_nodes(C)))+'\n') #***Average degree of charged residues (R,D,E,H,K)
			w.write(str(nx.average_clustering(C))+'\n') #***Average local clustering coefficient of charged residues
			w.write(str(nx.global_reaching_centrality(C))+'\n') #global, Average local reaching centrality of charged residues
	
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
	#print(input_nx(file)) 
	#name_base(path_pdb,file)
	#alpha_content(file)
	#beta_content(file)
	#load_hydrophobic(file)
	#load_charged(file)
	#laplacian_walk(G)
	features_vector(file)
	
	#inpnetx_obj = input_netx()
	#inpnetx_obj.main(parser)

