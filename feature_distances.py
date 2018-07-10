# Aim: to calculate feature distances
import os
import fileinput
from glob import glob
import linecache
import math
from itertools import combinations
import filemapper as fm
import numpy as np
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
import itertools
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt

path_fv = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/feature_vectors/'
path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_files/'

def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list
			
def create_array(file):
	base_list = name_base(path_pdb,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_fv+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			#temp.extend(fv)
			x = np.array(fv)
			#y = np.append([x],[x], axis=0)
			temp.append(x)
			
	return temp
			#print x
			#dst = distance.euclidean(z,z)
			#print dst
			#m, n = np.meshgrid(x, x)
			#out = abs(m-n) 
			#print m,n

def calc_dist1():
	base_list = name_base(path_pdb,file)
	vec = create_array(file)
	#print vec

	vec = np.asarray(vec)
	#print vec
	
	#mat = [[[[abs(vec[i][l]-vec[j][l] for l in len(vec[i]))] for k in len(vec)] for j in len(vec)] for i in len(vec)]
	#print mat
	
	y = pdist(vec, 'euclidean')
	#print y
	z = squareform(y)
	print z
	
def k_means():
	z = calc_dist1()
	
	kmeans = KMeans(n_clusters = 2, random_state = 0).fit(z)
	plt.show()
	#print kmeans
	
def spectral_cluster():
	
					
if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_array(file)
	#calc_dist1()
	k_means()
	#spectral_cluster()

