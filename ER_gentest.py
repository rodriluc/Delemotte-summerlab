#!/usr/bin/env
from __future__ import division
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
import scipy
import scipy.sparse
import mdtraj as md
import itertools
import time
from igraph import*



import MD_cmaps

path_data = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/Results_data/'
path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_edited/'
feat_path = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/feature_vectors/'
path_cmap = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/network_files/'

os.path.splitext(path_pdb)[0]
from os.path import basename
	
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
	plt.savefig(path_cmap+'cmap_'+base+'.svg')
	#print 'Network has been saved as *.svg' 
	return G
			
def laplacian_walk(G):

	lap = nx.normalized_laplacian_matrix(G, weight='weight')
	lap = np.eye(lap.shape[0])-lap #random walk: identity matrix - lap
	eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
	lap_sum = (eigenvalues[1]-eigenvalues[0]) #difference 1-second largest
	return lap_sum

def erdos_renyi_grc(G):
	
	grc_test = []
	N = nx.number_of_nodes(G)
	E = nx.number_of_edges(G)
	prob_edges = ((N*(N-1))/2)/E
	for i in range(500):
		er_graph = nx.erdos_renyi_graph(N,prob_edges,seed=10) #number of nodes, probability for edge creation
		grc = nx.global_reaching_centrality(er_graph)
		grc_test.append(grc)
	print grc_test
	#return np.average(grc_test)
	
def erdos_renyi_lap(G): #ER build graphs based on individual proteins, so number of nodes and probability of those edges being created

	lap_test = []
	lap_test1 = []
	lap_test2 = []
	N = nx.number_of_nodes(G)
	E = nx.number_of_edges(G)
	prob_edges = (((N*(N-1))/2)/E)/100 #divide by 100 float not percentage
	#print prob_edges
	'''start_time = time.time()
	for i in range(500):
		er_graph = nx.gnp_random_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test1.append(lap_sum)
	print lap_test1
	print ("--- %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	for i in range(500):
		er_graph = nx.fast_gnp_random_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test2.append(lap_sum)
	print lap_test2
	print ("--- %s seconds ---" % (time.time() - start_time))'''
	start_time = time.time()
	for i in range(500):
		er_graph = nx.erdos_renyi_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test.append(lap_sum)
	print lap_test
	print ("--- %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	for i in range(500):
		er_graph = nx.erdos_renyi_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues, eigenvectors = np.linalg.eigh(lap)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test1.append(eigenvalues)
	print lap_test1
	print ("--- %s seconds ---" % (time.time() - start_time))
	# using IGRAPH
	start_time = time.time()
	for i in range(500):
		er_graph = Graph.Erdos_Renyi(n=N,p=prob_edges, directed=True, loops=False) #number of nodes, probability for edge creation
		lap = er_graph.laplacian(normalized=True)
		e = np.linalg.eigvals(lap)
		lap_test2.append(e)
	print lap_test2
	print ("--- %s seconds ---" % (time.time() - start_time))
	start_time = time.time()

	for i in range(500):
		er_graph = nx.gnp_random_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		'''lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test1.append(lap_sum)'''
	print er_graph
	print ("--- %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	for i in range(500):
		er_graph = nx.fast_gnp_random_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		'''lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test2.append(lap_sum)'''
	print er_graph
	print ("--- %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	for i in range(500):
		er_graph = nx.erdos_renyi_graph(N,prob_edges,seed=None) #number of nodes, probability for edge creation
		'''lap = nx.normalized_laplacian_matrix(er_graph, weight='weight')
		lap = np.eye(lap.shape[0])-lap 
		eigenvalues,eigenvectors = scipy.sparse.linalg.eigsh(lap,k=2)
		lap_sum = (eigenvalues[1]-eigenvalues[0]) 
		lap_test.append(lap_sum)'''
	print er_graph
	print ("--- %s seconds ---" % (time.time() - start_time))

	#return np.average(lap_test)
	
def features_vector(file): 

	base_list = name_base(path_pdb,file)

	for i in range(len(base_list)):
		G = input_nx(path_data,base_list[i])
		L = laplacian_walk(G)
		er_lap = erdos_renyi_lap(G)
		er_grc = erdos_renyi_grc(G)

		print base_list[i]
		try:
			'''print 'GRC'
			print(str(nx.global_reaching_centrality(G))+'\n') 
			print(str(er_grc)+'\n') 
			print(str(nx.global_reaching_centrality(G)/er_grc)+'\n') '''
			print 'Ion Spectra'
			print(str(L)+'\n') 
			print(str(er_lap)+'\n')
			print(str(L/er_lap)+'\n')

		except nx.exception.NetworkXError: 
			pass #graph is not connected


if __name__ == '__main__':
	#laplacian_walk(G)
	#erdos_renyi_grc()
	#erdos_renyi_lap()
	features_vector(file)

