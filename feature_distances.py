# Aim: to calculate feature distances
import os
import fileinput
from glob import glob
import linecache
import math
from itertools import combinations, cycle, islice
import filemapper as fm
import numpy as np
np.set_printoptions(threshold=np.inf)

from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
import itertools
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn import metrics

from matplotlib import pyplot as plt
import pandas as pd
from functools import partial
import seaborn
import collections
import spectral
#from spectral import utils, affinity, clustering 


path_fv = '/home/lrodriguez/Delemotte-summerlab_ERnorm_100ER_4.5A/feature_vectors/'
path_pdb = '/home/lrodriguez/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_files/'

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
			#print x
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

	return z
	
def k_means():
	#z = calc_dist1()
	vec = create_array(file)
	vec = np.asarray(vec)

	#PCA-reduced data
	reduced_data = PCA(n_components=2).fit_transform(vec) #n_components=3

	model = KMeans(4)
	model.fit(reduced_data) #vec for k-means full
	clust_labels = model.labels_
	centroids = model.cluster_centers_
	
	
	kmeans = pd.DataFrame(clust_labels)
	#print(vec.shape)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(reduced_data[:,0],reduced_data[:,1],s=50,c=clust_labels) #when pca comp. = 3 then 1,2
	ax.set_title('K-means clustering')
	#plt.colorbar(scatf)
	plt.show()

	np.savetxt('cluster_indices.txt',clust_labels)

	
def spectral_cluster():

	Z = calc_dist1()
	#vec = np.asarray(vec)	
	#vec = np.matrix(vec) #precomputed, affinity matrix should pass through Gaussian kernel

	#gk = np.exp(-Z**2/(2.*delta**2))
	#print gk

	'''sc = SpectralClustering(n_clusters=4, affinity='precomputed')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(sc[:,0],sc[:,1],s=50,c=clust_labels)
	ax.set_title('Spectral clustering')
	plt.show()'''

	#Spectral Clustering example
	sc = SpectralClustering(eigen_solver='arpack',affinity='precomputed',assign_labels='discretize').fit(Z) #affinity='rbf'
	print 'Spectral Clustering'
	print collections.Counter(sc.labels_)
	print metrics.silhouette_score(Z, sc.labels_)

	reduced_data = PCA(n_components=2).fit_transform(Z)
	#plot_2d_data(reduced_data,sc.labels_)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(reduced_data[:,0],reduced_data[:,1],s=50,c=sc.labels_)
	ax.set_title('Spectral clustering')
	plt.show()
	
					
if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_array(file)
	#calc_dist1()
	k_means()
	#spectral_cluster()








