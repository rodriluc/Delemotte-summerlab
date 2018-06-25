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
path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/'


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


def features_vector():
	G = input_nx()
	g = nx.DiGraph(G)
	num_nodes = nx.number_of_nodes(G)
	num_edges = nx.number_of_edges(G)
	
	#Average degree of hydrophobic residues (F,M,W,I,V,L,P,A)
	for file in os.listdir(path_pdb):
		#print file
		if file.endswith('.pdb'): 
			iter_file = open(file)
			#with open(file) as iter_file:
			lines = iter_file.readlines()
			head = lines[0]
			#print iter_file
			base = file[:4] 
			top = md.load_pdb(file).topology
			with open ('hydrophobic_'+base+'.pdb', 'w') as w:
				w.write(head)
				#print top
				hydrophobic_list = top.select('resname PHE MET TRP ILE VAL LEU PRO ALA')
				new_list = [x+1 for x in hydrophobic_list]
				#print new_list
				T = [lines[i] for i in new_list] #new_list
				for line in T:
					if line.rstrip():
						w.write(line) 
				w.write('END')
				w.write('\n')
				print 'Hydrophobic file saved!'
					

	'''for file in os.listdir(path_data):
		if file.endswith('.pdb'):
			base = file[:4]
			with open ('features_'+base+'.txt', 'w') as w:
				top = md.load_pdb(file).topology
				#resid = ['PHE','MET','TRP','ILE','VAL','LEU','PRO','ALA']
				#res_list = []
				
				#for residues in top.topology.residues:
				hydrophobic_list = top.select('resname PHE MET TRP ILE VAL LEU PRO ALA')
				print hydrophobic_list
				#for line in file:
					#if re.match(hydrophobic_list, line):
						
				#cm = cm(hydrophobic_list)
				#print cm
				#for residues in top.select('resname PHE'):
					#print residues
					#res_list.append(residue)
					#hydrophobic_list = [name for name in res_list if (name[0:3] in resid)]
					#if residues == top.select('resname PHE').all():
						#print residues
					#if residue != hydrophobic_list:
						#hydrophobic_list.append(residue)
						#pho = residue #everything but hydrophobic
						#G = nx.from_numpy_matrix(B)
				pho = nx.Graph(Phe)
				#pho.add_nodes_from(hydrophobic_list) #remove_nodes_from
				nx.draw(pho)
				plt.savefig('hydrophobic.svg')
				print(top.select('resname PHE'))
				print(nx.number_of_edges(pho))
				print(nx.number_of_nodes(pho))
				w.write('hydrophobic avg degree below'+'\n') #test where result is
				w.write(str(float(nx.number_of_edges(pho))/float(nx.number_of_nodes(pho)))+'\n')
					
					#else:
				w.write(str(float(num_edges)/float(num_nodes))+'\n') #avg degree
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

				
				top = [residue for residue in top.chain().residues if residue == hydrophobic_resid]
				print (top)
				topology = md.load_pdb(file).topology
				table, bonds = topology.to_dataframe()
				print(table.head())
				a = [d for n,d in G.nodes_iter(data=True)]
				
				g = nx.Graph()
				g.add_nodes_from(self.atoms)
				#g.add_edges_from(self.bonds)
				#print(top.select('resname Phe'))
				if 	top.select('resname == Phe'):
					top == 1 #assign 1 if hydrophobic
				else:
					top == 0 #assign 0 if hydrophilic
				g = nx.Graph(top)
				w.write(str(float(nx.number_of_edges(g))/float(nx.number_of_nodes(g)))+'\n')'''
				#***Average degree of charged residues (R,D,E,H,K)
				#***Average local clustering coefficient of hydrophobic residues
				#***Average local clustering coefficient of charged residues
				#***Average local reaching centrality of hydrophobic residues
				#***Average local reaching centrality of charged residues
				#***Secondary structure content (helix content + beta strand content)
	
	print 'Features vector *.txt file has been saved!'''

if __name__ == '__main__':
	#install('networkx') #satisfied
	#print(input_nx()) 
	#print(attributes_graph())
	print(features_vector())
