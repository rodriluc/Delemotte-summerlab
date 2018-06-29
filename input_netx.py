#!/usr/bin/env python

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

path_data = 'Results_data/'
path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_edited/'
feat_path = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/feature_vectors/'
path_hydro = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/hydrophobic_files/'
path_charg = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/charged_files/'

os.path.splitext(path_pdb)[0]
from os.path import basename

def install(package):
	subprocess.call([sys.executable, "-m", "pip", "install", package])
	
def __init__():
	return

def input_nx(): 
	
	for file in os.listdir(path_data):
		#print file
		if file.startswith('cmap_processed_'): #'/cmap_processed_'

			#A = genfromtxt(file, delimiter=' ') #usecols=range(0,ncols-1)
			# can also use genfromtxt(file) instead of loadtxt
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) #usecols=(0, -1) blank line at end
			B = np.matrix(np.array(A))
			G = nx.from_numpy_matrix(B)
			#print B
			#G = nx.Graph(G)
			nx.draw(G)

			plt.show()
			plt.savefig('cmap.svg' )
			print 'Figure has been saved as *.svg' #called twice so prints twice
			return G
			
def attributes_graph(n = None):
	G = input_nx()
	#G = nx.Graph(G)
	num_nodes = nx.number_of_nodes(G)
	num_edges = nx.number_of_edges(G)
	avg = float(num_edges)/float(num_nodes)
	print '***Statistic attributes of networks***', '\n'
	print 'Number of nodes: ',num_nodes
	print 'Number of edges: ',num_edges
	print 'Average degree: ',avg
	print 'Average shortest path length: ',nx.average_shortest_path_length(G)
	print 'Diameter (max. shortest path): ',nx.diameter(G)
	print 'Radius (min. shortest path): ',nx.radius(G)
	print 'Average local clustering: ',nx.clustering(G) #local
	print 'Average global clustering: ',nx.average_clustering(G) #global
	#Eigenvalue of laplacian, L = D - A
	print 'Laplacian matrix of G', '\n'
	print nx.laplacian_matrix(G, weight='weight')
	#Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
	print 'Normalized laplacian matrix of G', '\n'
	print nx.normalized_laplacian_matrix(G, weight='weight')
	print 'Degree assortativity coefficient: ',nx.degree_assortativity_coefficient(G)
	
def name_base():
	
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'):
			#base = file[:4]
			#return base
			basename = file.split('.')[:-1]
			base =''.join(basename)
			return base
		
def alpha_content():

	base = name_base()
	for file in os.listdir(path_pdb):
		#print file
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'): 
			#print file
			iter_file = open((path_pdb+file), 'r')
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines:
				#print line
				if cryst in line:
					#print line
					for i in range(len(line)):
						#return line[12:18]
						return line[9:15]
			
def beta_content():
	base = name_base()
	for file in os.listdir(path_pdb):
		#print file
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'): 
			#print file
			iter_file = open((path_pdb+file), 'r')
			lines = iter_file.readlines()
			cryst = ('CRYST1')
			for line in lines: 
				#print line
				if cryst in line:
					#print line
					for i in range(len(line)):
						#return line[21:27]
						return line[18:24]

def load_hydrophobic():
	base = name_base()

	for file in os.listdir(path_data):
		if file.startswith('cmap_processed_hydrophobic'): #'cmap_processed_hydrophobic'
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) 
			B = np.matrix(np.array(A))
			H = nx.from_numpy_matrix(B)
			nx.draw(H)
			plt.show()
			plt.savefig('cmap_hydrophobic.svg' )
			#print 'Figure has been saved as *.svg' 
			return H
				
def load_charged():
	base = name_base()
	
	for file in os.listdir(path_data):
		if file.startswith('cmap_processed_charged'): #'cmap_processed_charged'
			A = np.loadtxt(path_data+file, dtype=float, unpack=True) 
			B = np.matrix(np.array(A))
			C = nx.from_numpy_matrix(B)
			nx.draw(C)
			plt.show()
			plt.savefig('cmap_charged.svg' )
			#print 'Figure has been saved as *.svg' 
			return C
	
def features_vector():
	base = name_base()
	G = input_nx()
	H = load_hydrophobic()
	C = load_charged()
	A = alpha_content()
	B = beta_content()

	with open ((feat_path+('features_'+base+'.txt')), 'w') as w:
		for file in os.listdir(path_pdb):
			w.write(str(float(nx.number_of_edges(G))/float(nx.number_of_nodes(G)))+'\n') #avg degree
			w.write(str(nx.average_shortest_path_length(G))+'\n') #avg shortest path length
			w.write(str(nx.diameter(G))+'\n') #diameter (max. shortest path)
			w.write(str(nx.radius(G))+'\n') #radius (min. shortest path)
			#w.write(str(nx.clustering(G))) #local clustering coefficient (individual nodes)
			w.write(str(nx.average_clustering(G))+'\n') #global clustering coefficient
			#***Number of quasi-rigid domains
			#w.write(str(nx.laplacian_matrix(G, weight='weight'))+'\n') #Eigenvalue of laplacian, L = D - A
			w.write(str(nx.normalized_laplacian_matrix(G, weight='weight'))+'\n') #Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
			w.write(str(nx.degree_assortativity_coefficient(G))+'\n') #assortativity coefficient
			w.write(str(nx.degree_centrality(G))+'\n') #global reaching centrality
			#***Residue intrinsic dimensionality (may be used to compute 6.?)
			#***Secondary structure content (helix content + beta strand content)
			w.write(str(A)+'\n') #alpha content
			w.write(str(B)+'\n') #beta content
		for file in os.listdir(path_hydro):
			#***Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
			w.write(str(float(nx.number_of_edges(H))/float(nx.number_of_nodes(H)))+'\n')
			#***Average local clustering coefficient of hydrophobic residues
			w.write(str(nx.clustering(H))+'\n')
			#***Average local reaching centrality of hydrophobic residues
			w.write(str(nx.local_reaching_centrality(H,3))+'\n') #needs second parameter, the node
		for file in os.listdir(path_charg):
			#***Average degree of charged residues (R,D,E,H,K)
			w.write(str(float(nx.number_of_edges(C))/float(nx.number_of_nodes(C)))+'\n')
			#***Average local clustering coefficient of charged residues
			w.write(str(nx.clustering(C))+'\n')
			#***Average local reaching centrality of charged residues
			w.write(str(nx.local_reaching_centrality(C,3))+'\n')
			
			'''if file.endswith('.pdb'):
				if file.startswith('hydrophobic'):
					#***Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
					w.write(str(float(nx.number_of_edges(H))/float(nx.number_of_nodes(H)))+'\n')
					#***Average local clustering coefficient of hydrophobic residues
					w.write(str(nx.clustering(H))+'\n')
					#***Average local reaching centrality of hydrophobic residues
					w.write(str(nx.local_reaching_centrality(H,3))+'\n') #needs second parameter, the node
				elif file.startswith('charged'):
					#***Average degree of charged residues (R,D,E,H,K)
					w.write(str(float(nx.number_of_edges(C))/float(nx.number_of_nodes(C)))+'\n')
					#***Average local clustering coefficient of charged residues
					w.write(str(nx.clustering(C))+'\n')
					#***Average local reaching centrality of charged residues
					w.write(str(nx.local_reaching_centrality(C,3))+'\n')
				else:
					w.write(str(float(nx.number_of_edges(G))/float(nx.number_of_nodes(G)))+'\n') #avg degree
					w.write(str(nx.average_shortest_path_length(G))+'\n') #avg shortest path length
					w.write(str(nx.diameter(G))+'\n') #diameter (max. shortest path)
					w.write(str(nx.radius(G))+'\n') #radius (min. shortest path)
					#w.write(str(nx.clustering(G))) #local clustering coefficient (individual nodes)
					w.write(str(nx.average_clustering(G))+'\n') #global clustering coefficient
					#***Number of quasi-rigid domains
					#w.write(str(nx.laplacian_matrix(G, weight='weight'))+'\n') #Eigenvalue of laplacian, L = D - A
					w.write(str(nx.normalized_laplacian_matrix(G, weight='weight'))+'\n') #Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
					w.write(str(nx.degree_assortativity_coefficient(G))+'\n') #assortativity coefficient
					w.write(str(nx.degree_centrality(G))+'\n') #global reaching centrality
					#***Residue intrinsic dimensionality (may be used to compute 6.?)
					#***Secondary structure content (helix content + beta strand content)
					w.write(str(A)+'\n') #alpha content
					w.write(str(B)+'\n') #beta content'''
	
	print 'Features vector *.txt file has been saved! (original PDB, hydrophobic, charged)'

if __name__ == '__main__':
	#install('networkx') #satisfied
	#print(input_nx()) 
	#print(attributes_graph())
	#name_base()
	#alpha_content()
	#beta_content()
	#load_hydrophobic()
	#load_charged()
	features_vector()
