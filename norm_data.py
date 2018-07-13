# Aim: create barplots and possibly normalize data based on random graph

import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing
from sklearn.preprocessing import Normalizer
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
import itertools
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn import metrics
import collections

path_fv = '/home/lrodriguez/Delemotte-summerlab/feature_vectors/'
path_pdb = '/home/lrodriguez/Delemotte-summerlab/PDB_edited/'

def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list
	#print base_list		
def create_vlist():
	base_list = name_base(path_pdb,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_fv+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			x = np.array(fv)

			temp.append(x)
			
	return temp
	

def create_barplots():
	base_list = name_base(path_pdb,file)
	fv_list = create_vlist()
	
	#min_max_scaler = preprocessing.MinMaxScaler()
	#fv_minmax = min_max_scaler.fit_transform(fv_list) #normalize data set 0-1, scales features to range

	fv_list = np.matrix(fv_list)
	fv_norm = preprocessing.normalize(fv_list,norm='l1') #normalize data, sum of array equals to 1


	Tot = 15
	Cols = 5
	Rows = Tot//Cols
	Rows += Tot%Cols
	Position = range(1,Tot+1)	
	
	for item in fv_norm[:15]:
		width = 0.5
		N = len(item)
		x = range(N)

	fig = plt.figure(1)
	for i in range(Tot):
		ax=fig.add_subplot(Rows,Cols,Position[i])
		ax.bar(x, fv_norm[i], width, align='center')
		plt.ylim([0,0.9])

	plt.show()


def cluster_norm():
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)
	
	vec = preprocessing.normalize(vec,norm='l1')
	#PCA-reduced data
	reduced_data = PCA(n_components=2).fit_transform(vec) #n_components=3

	model = KMeans(4) #3?
	model.fit(reduced_data) #vec for k-means full
	clust_labels = model.labels_
	centroids = model.cluster_centers_
	print 'K-means Clustering'
	print collections.Counter(clust_labels)
	print metrics.silhouette_score(vec, clust_labels)

	kmeans = pd.DataFrame(clust_labels)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(reduced_data[:,0],reduced_data[:,1],s=50,c=clust_labels) #when pca comp. = 3 then 1,2
	ax.set_title('K-means clustering Normalized')
	plt.show()

	np.savetxt('norm_cluster_indices.txt',clust_labels)

def calc_dist1():
	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	#Normalize arrays before calculating distance
	vec = np.matrix(vec) #might not need to force matrix
	vec = preprocessing.normalize(vec,norm='l1')

	y = pdist(vec, 'euclidean')
	z = squareform(y)
	return z
	
def spectral_cluster():
	Z = calc_dist1()

	sc = SpectralClustering(eigen_solver='arpack',affinity='precomputed',assign_labels='discretize').fit(Z) #affinity='rbf'
	print 'Spectral Clustering'
	print collections.Counter(sc.labels_)
	print metrics.silhouette_score(Z, sc.labels_)

	reduced_data = PCA(n_components=3).fit_transform(Z) #2

	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(reduced_data[:,1],reduced_data[:,2],s=50,c=sc.labels_)
	ax.set_title('Spectral clustering Normalized')
	plt.show()


if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_vlist()
	#create_barplots()
	cluster_norm()
	#spectral_cluster()

