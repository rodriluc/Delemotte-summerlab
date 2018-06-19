import networkx as nx
import numpy as np
from numpy import genfromtxt
import os
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from MD_cmaps import MD_cmaps
import subprocess
import sys
#PATH=$PATH :/home/l/u/lucier/.local/bin/ #for bash jupyter notebook need to link path
import scipy
import scipy.sparse
'''if nodelist is None:
	nodelist = G.nodes()
A = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, weight=weight, format='csr')
n,m = A.shape
diags = A.sum(axis=1).flatten()
D = scipy.sparse.spdiags(diags, [0], m, n, format='csr')
L = D - A
with scipy.errstate(divide='ignore'):
	diags_sqrt = 1.0/scipy.sqrt(diags)
diags_sqrt[scipy.isinf(diags_sqrt)] = 0
DH = scipy.sparse.spdiags(diags_sqrt, [0], m, n, format='csr')
return DH.dot(L.dot(DH))'''

def install(package):
	subprocess.call([sys.executable, "-m", "pip", "install", package])

def input_nx(): 
	#cmap = MD_cmaps().getContactMap()
	for file in os.listdir('/afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/'):
		if file.startswith('cmap_processed_'): #distance_matrix or cmap_processed

			#A = genfromtxt(file, delimiter=' ') #usecols=range(0,ncols-1)
			'''A = np.fromfile(file, dtype=float)
			B = A.shape
			print(B)
			G = nx.from_numpy_matrix(A)'''
			
			# can also use genfromtxt(file) instead of loadtxt
			A = np.loadtxt(file, dtype=float, unpack=True) #usecols=(0, -1) blank line at end
			B = np.matrix(np.array(A))
			#print(B)
			G = nx.from_numpy_matrix(B)
			#G = nx.Graph(G)
			nx.draw(G)
			
			#label nodes with amino acid residues
			'''labels = {}    
			for node in G.nodes():
				if node in dict_labels:
					#set the node name as the key and the label as its value 
					labels[node] = node
			nx.draw_networkx_labels(G,pos,labels,font_size=16,font_color='r')'''
			
			plt.show()
			plt.savefig('cmap.svg' )
			print 'Figure has been saved as *.svg' #called twice so prints twice
			return G
			
def attributes_graph(n = None):
	G = input_nx()
	g = nx.DiGraph(G)
	num_nodes = nx.number_of_nodes(G)
	num_edges = nx.number_of_edges(G)
	print '***Statistic attributes of networks***', '\n'
	print 'Number of nodes: ',num_nodes
	print 'Number of edges: ',nx.number_of_edges(G)
	print 'Average degree: ', float(num_edges)/float(num_nodes)
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

	print float(num_edges)/float(num_nodes) #avg degree
	print nx.average_shortest_path_length(G) #avg shortest path length
	print nx.diameter(G) #diameter (max. shortest path)
	print nx.radius(G) #radius (min. shortest path)
	print nx.clustering(G) #local clustering coefficient (individual nodes)
	print nx.average_clustering(G) #global clustering coefficient
	
	print nx.laplacian_matrix(G, weight='weight') #Eigenvalue of laplacian, L = D - A
	print nx.normalized_laplacian_matrix(G, weight='weight') #Normalized laplacian matrix, N = D^{-1/2} L D^{-1/2}
	print nx.degree_assortativity_coefficient(G) #assortativity coefficient
	

if __name__ == '__main__':
	#install('networkx') #satisfied
	print(input_nx()) 
	print(attributes_graph())
	print(features_vector())
