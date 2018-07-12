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

path_fv = '/home/lrodriguez/Delemotte-summerlab/feature_vectors/'
path_pdb = '/home/lrodriguez/Delemotte-summerlab/PDB_files/'

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
		plt.ylim([0,0.66])

	plt.show()

def create_array(file):
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

def cluster_norm():
	vec = create_array(file)
	vec = np.asarray(vec)
	vec = np.matrix(vec)
	
	vec = preprocessing.normalize(vec,norm='l1')
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
	ax.set_title('K-means clustering: Normalized Data')
	plt.show()


if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_vlist()
	create_barplots()
	#cluster_norm()

