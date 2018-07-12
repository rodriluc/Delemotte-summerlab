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

path_data = '/home/lrodriguez/Delemotte-summerlab/Results_data/'
path_pdb = '/home/lrodriguez/Delemotte-summerlab/PDB_edited/'
feat_path = '/home/lrodriguez/Delemotte-summerlab/feature_vectors/'
path_hydro = '/home/lrodriguez/Delemotte-summerlab/hydrophobic_files/'
path_charg = '/home/lrodriguez/Delemotte-summerlab/charged_files/'
path_cmap = '/home/lrodriguez/Delemotte-summerlab/network_files/'

os.path.splitext(path_pdb)[0]
from os.path import basename


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
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def input_nx(path_data,base): 

	A = np.loadtxt(path_data+'cmap_processed_'+base+'.txt', dtype=float, unpack=True) 
	B = np.matrix(np.array(A))
	G = nx.from_numpy_matrix(B)

	nx.draw(G)

	plt.show()
	plt.savefig(path_cmap+'cmap_'+base+'.svg') #'cmap.svg'
	print 'Figure has been saved as *.svg' #called twice so prints twice
	return G

def alpha_content(path_pdb,base):
	
	for file in os.listdir(path_pdb):
		with open(path_pdb+base+'.pdb') as iter_file:
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines:
				col = line.split()
				if cryst in line:
					return col[1]
			
def beta_content(path_pdb,base):
	
	for file in os.listdir(path_pdb):	
		with open(path_pdb+base+'.pdb') as iter_file:
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines:
				col = line.split()
				if cryst in line:
					return col[2]

def load_hydrophobic(path_data,base):
	
	A = np.loadtxt(path_data+'cmap_processed_hydrophobic_'+base+'.txt', dtype=float, unpack=True) 
	B = np.matrix(np.array(A))
	H = nx.from_numpy_matrix(B)
	nx.draw(H)
	plt.show()
	plt.savefig(path_cmap+'cmap_hydrophobic_'+base+'.svg')
	return H
				
def load_charged(path_data,base):
	
	A = np.loadtxt(path_data+'cmap_processed_charged_'+base+'.txt', dtype=float, unpack=True) 
	B = np.matrix(np.array(A))
	C = nx.from_numpy_matrix(B)
	nx.draw(C)
	plt.show()
	plt.savefig(path_cmap+'cmap_charged_'+base+'.svg')
	return C
			
def laplacian_walk(G):

	lap = nx.normalized_laplacian_matrix(G, weight='weight')
	lap = np.eye(lap.shape[0])-lap #random walk: identity matrix - lap
	eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
	lap_sum = (eigenvalues[1]-eigenvalues[0]) #difference 1-second largest
	return lap_sum
	
def features_vector(file): #self

	base_list = name_base(path_pdb,file)

	for i in range(len(base_list)):
		G = input_nx(path_data,base_list[i])
		L = laplacian_walk(G)
		H = load_hydrophobic(path_data,base_list[i])
		C = load_charged(path_data,base_list[i])
		A = alpha_content(path_pdb,base_list[i])
		B = beta_content(path_pdb,base_list[i])
		with open ((feat_path+('features_'+base_list[i]+'.txt')), 'w') as w: 
			#w.write(str(float(nx.number_of_edges(G))/float(nx.number_of_nodes(G)))+'\n') #removed avg degree
			w.write(str(nx.average_shortest_path_length(G)/(((log(float(nx.number_of_nodes(G))))-0.5772)/(log((float(nx.number_of_edges(G)))/float(nx.number_of_nodes(G)))+0.5))+'\n') #avg shortest path length
			w.write(str(log((float(nx.number_of_nodes(G))))/log((float(nx.number_of_edges(G))/float(nx.number_of_nodes(G)))))+'\n') #diameter (max. shortest path)
			w.write(str(nx.radius(G)/nx.average_shortest_path_length(G))+'\n') #radius (min. shortest path)
			w.write(str(nx.average_clustering(G)/((nx.average_shortest_path_length(G)/float(nx.number_of_nodes(G)))))+'\n') #clustering coefficient of random graph
			#***Number of quasi-rigid domains
			w.write(str(nx.degree_assortativity_coefficient(G))+'\n') #assortativity coefficient
			w.write(str(nx.global_reaching_centrality(G))+'\n') #nx.degree_centrality(G)
			#***Residue intrinsic dimensionality (may be used to compute 6.?)
			w.write(str(L)+'\n') #Normalized laplacian walk from function, N = D^{-1/2} L D^{-1/2}
			w.write(str(A/(A+B))+'\n') #alpha content 
			w.write(str(B/(A+B))+'\n') #beta content
	
			w.write(str(float(nx.number_of_edges(H))/float(nx.number_of_nodes(H)))+'\n') #***Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
			w.write(str(nx.average_clustering(H)/((nx.average_shortest_path_length(H)/float(nx.number_of_nodes(H)))))+'\n') #clustering coefficient of random graph - hydrophobic residues
			w.write(str(nx.global_reaching_centrality(H))+'\n') #global, Average local reaching centrality of hydrophobic residues
	
			w.write(str(float(nx.number_of_edges(C))/float(nx.number_of_nodes(C)))+'\n') #***Average degree of charged residues (R,D,E,H,K)
			w.write(str(nx.average_clustering(C)/((nx.average_shortest_path_length(C)/float(nx.number_of_nodes(C)))))+'\n') #clustering coefficient of random graph - charged residues
			w.write(str(nx.global_reaching_centrality(C))+'\n') #global, Average local reaching centrality of charged residues
	
	print 'Features vector *.txt file has been saved! (original PDB, hydrophobic, charged)'


if __name__ == '__main__':
	#print(input_nx(file)) 
	#name_base(path_pdb,file)
	#alpha_content(file)
	#beta_content(file)
	#load_hydrophobic(file)
	#load_charged(file)
	#laplacian_walk(G)
	features_vector(file)

