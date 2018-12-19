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
from sklearn.cluster import KMeans, SpectralClustering, SpectralCoclustering, AffinityPropagation
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn import metrics
import collections
from scipy.cluster.hierarchy import fclusterdata
from sklearn import datasets
import seaborn as sns
import hdbscan
from mpl_toolkits.mplot3d import Axes3D

path_fv = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/feature_vectors/'
path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_edited/'
path_gpcr = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_gpcr/'
path_ic = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_ionchannel/'
path_enz = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_enzyme/'
path_kin = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_kinase/'

def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list
	#print base_list

def gpcr_base(path_gpcr,file):
	base_list = []
	for file in os.listdir(path_gpcr):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def ic_base(path_ic,file):
	base_list = []
	for file in os.listdir(path_ic):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list	

def enz_base(path_enz,file):
	base_list = []
	for file in os.listdir(path_enz):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list

def kin_base(path_kin,file):
	base_list = []
	for file in os.listdir(path_kin):
		if file.endswith('.txt'):
			basename = file.split('_')[-1].split('.')[0]
			base =''.join(basename)
			base_list.append(base)
	return base_list
		
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
	
def create_gpcr_list():
	base_list = gpcr_base(path_gpcr,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_gpcr+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			x = np.array(fv)

			temp.append(x)
			
	return temp

def create_ic_list():
	base_list = ic_base(path_ic,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_ic+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			x = np.array(fv)

			temp.append(x)
			
	return temp

def create_enz_list():
	base_list = enz_base(path_enz,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_enz+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			x = np.array(fv)

			temp.append(x)
			
	return temp

def create_kin_list():
	base_list = kin_base(path_kin,file)
	temp = []
	for i in range(len(base_list)):
		
		with open (path_kin+'features_'+base_list[i]+'.txt') as fv:
			lines = fv.readlines()
			fv = [x.strip() for x in lines]
			fv = [float(i) for i in fv]

			x = np.array(fv)

			temp.append(x)
			
	return temp

def create_barplots():
	base_list = name_base(path_pdb,file)
	fv_list = create_vlist()

	fv_list = np.matrix(fv_list)
	fv_norm = preprocessing.normalize(fv_list,norm='l1')


	Tot = 63
	Cols = 6
	Rows = Tot//Cols
	Rows += Tot%Cols
	Position = range(1,Tot+1)	
	
	for item in fv_norm:
		width = 0.5
		N = len(item)
		x = range(N)

	fig = plt.figure(1)
	for i in range(Tot):
		ax=fig.add_subplot(Rows,Cols,Position[i])
		fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.9)
		ax.bar(x, fv_norm[i], width, align='center')
		ax.set_title(base_list[i])
		plt.ylim([0,0.1])

	plt.show()

def avg_barplot():
	fv_list = create_vlist()

	fv_list = np.matrix(fv_list)
	x = preprocessing.normalize(fv_list,norm='l1')	
	
	for item in x:
		width = 0.5
		N = len(item)
		A = range(N)
	x = np.array(x)
	y = x.mean(axis=0)
	
	error = x.std(axis=0)
	#print error
	barwidth = 0.25
	Y = ['avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew']
	X = np.arange(len(Y))
	fig,ax = plt.subplots()
	ax.bar(A, y, yerr=error, align='center')
	plt.xticks(X,Y)
	plt.ylim([0,0.08])

	plt.show()

def avg_byfeat_barplot():
	gpcr = create_gpcr_list()
	ic = create_ic_list()
	enz = create_enz_list()
	kin = create_kin_list()

	gpcr = np.matrix(gpcr)
	x = preprocessing.normalize(gpcr,norm='l1')	
	for item in x:
		width = 0.5
		N = len(item)
		A = range(N)
	x = np.array(x)
	y = x.mean(axis=0)
	error = x.std(axis=0)

	ic = np.matrix(ic)
	x1 = preprocessing.normalize(ic,norm='l1')	
	for item in x1:
		width = 0.5
		N1 = len(item)
		A1 = range(N1)
	x1 = np.array(x1)
	y1 = x1.mean(axis=0)
	error1 = x1.std(axis=0)

	enz = np.matrix(enz)
	x2 = preprocessing.normalize(enz,norm='l1')	
	for item in x2:
		width = 0.5
		N2 = len(item)
		A2 = range(N2)
	x2 = np.array(x2)
	y2 = x2.mean(axis=0)
	error2 = x2.std(axis=0)

	kin = np.matrix(kin)
	x3 = preprocessing.normalize(kin,norm='l1')	
	for item in x3:
		width = 0.5
		N3 = len(item)
		A3 = range(N3)
	x3 = np.array(x3)
	y3 = x3.mean(axis=0)
	error3 = x3.std(axis=0)
	
	barwidth = 0.20
	Y = ['avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew']
	X = np.arange(len(Y))
	X1 = [z +barwidth for z in X]
	X2 = [z +barwidth for z in X1]
	X3 = [z +barwidth for z in X2]
	fig,ax = plt.subplots()
	ax.bar(X, y, yerr=error, color = 'b', width=barwidth, edgecolor='white', label='GPCR') #gpcr
	ax.bar(X1, y1, yerr=error1, color = 'g', width=barwidth, edgecolor='white',label='Ion channel') #ion channels
	ax.bar(X2, y2, yerr=error2, color = 'r', width=barwidth, edgecolor='white',label='Enzyme') #enzyme
	ax.bar(X3, y3, yerr=error3, color = 'y', width=barwidth, edgecolor='white',label='Kinase') #kinase
	plt.xticks([w + barwidth for w in range(N)],Y)
	plt.ylim([0,0.08])
	plt.legend()

	plt.show()

def cluster_norm():

	base_list = name_base(path_pdb,file)

	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())

	vec = create_vlist()
	vec = np.asarray(vec)

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) #purple
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) #blue
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) #green
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) #yellow

	fam_array = np.asarray(fam_list) #cluster label

	vec = np.asarray(vec)
	vec = np.matrix(vec)
	#vec = preprocessing.normalize(vec,norm='l1')
	#print vec

	print(vec.shape)

	#PCA-reduced data
	pca = PCA(n_components=2).fit(vec)
	reduced_data = pca.transform(vec) #n_components=3
	print('Exp var: ' +str(np.sum(pca.explained_variance_ratio_)))
	x = reduced_data[:,0]
	x = (x-np.min(x))
	x = (x/np.max(x))
	#print reduced_data.shape()
	#print x.shape()
	reduced_data[:,0] = x

	x = reduced_data[:,1]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,1] = x

	'''x = reduced_data[:,2]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,2] = x'''

	model = KMeans(4).fit(reduced_data) #vec for k-means full
	y_pred = model.predict(reduced_data)
	clust_labels = model.labels_
	centroids = model.cluster_centers_

	'''print 'K-means Clustering'
	print collections.Counter(clust_labels)
	print metrics.silhouette_samples(vec, clust_labels)
	print metrics.silhouette_score(vec, clust_labels)'''
	#print metrics.mutual_info_score(vec, clust_labels)
	#print metrics.homogeneity_completeness_v_measure(vec, clust_labels)

	fig, ax=plt.subplots(1,2,figsize=(8,4)) 
	ax[0].scatter(reduced_data[:,0], reduced_data[:,1],c=fam_array) 
	ax[1].scatter(reduced_data[:,0], reduced_data[:,1],c=y_pred) #when pca comp. = 3 then 1,2

	ax[0].set_title('K-means clustering by protein structure')
	ax[1].set_title('K-means clustering (predicted)')

	ax[0].scatter(centroids[:,0], centroids[:,1],marker='x',s=80, linewidths=3,color='r')
	ax[1].scatter(centroids[:,0], centroids[:,1],marker='x',s=80, linewidths=3,color='r')

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

	sc = SpectralClustering(n_clusters=4,eigen_solver='arpack',affinity='rbf',assign_labels='discretize').fit(Z) 
	#sc = AffinityPropagation(damping=0.5, affinity='precomputed').fit(Z) #affinity='rbf'
	print 'Spectral Clustering'
	print collections.Counter(sc.labels_)
	print metrics.silhouette_samples(Z, sc.labels_)
	print metrics.silhouette_score(Z, sc.labels_)

	reduced_data = PCA(n_components=2).fit_transform(Z) #2

	fig = plt.figure()
	ax = fig.add_subplot(111)
	scatf = ax.scatter(reduced_data[:,0],reduced_data[:,1],s=50,c=sc.labels_)
	ax.set_title('Spectral clustering Normalized')

	plt.show()

	np.savetxt('norm_spec.cluster_indices.txt',sc.labels_)

def ClusterIndicesNumpy(clustNum, labels_array): #numpy 
	return np.where(labels_array == clustNum)[0]

def avg_feat_cluster(): #average of feature in each cluster
	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)
	
	vec = preprocessing.normalize(vec,norm='l1')

	#reduced_data = PCA(n_components=2).fit_transform(vec) 
	km = KMeans(4).fit(vec) #reduced_data

	clust_0 = ClusterIndicesNumpy(0, km.labels_) # all samples in cluster number __
	clust_1 = ClusterIndicesNumpy(1, km.labels_)
	clust_2 = ClusterIndicesNumpy(2, km.labels_)
	clust_3 = ClusterIndicesNumpy(3, km.labels_)
	vec0 = vec[clust_0] #extract all of my cluster number __ data points	
	vec1 = vec[clust_1] 
	vec2 = vec[clust_2]
	vec3 = vec[clust_3]

# Create average feat barplot for each cluster	(need to change vec1 for each cluster 1-4)

	for item in vec0:
		width = 0.5
		N4 = len(item)
		A4 = range(N4)
	x4 = np.array(vec0)
	y4 = vec0.mean(axis=0)
	error4 = x4.std(axis=0)
	for item in vec1:
		width = 0.5
		N5 = len(item)
		A5 = range(N5)
	x5 = np.array(vec1)
	y5 = vec1.mean(axis=0)
	error5 = x5.std(axis=0)
	for item in vec2:
		width = 0.5
		N6 = len(item)
		A6 = range(N6)
	x6 = np.array(vec2)
	y6 = vec2.mean(axis=0)
	error6 = x6.std(axis=0)
	for item in vec3:
		width = 0.5
		N7 = len(item)
		A7 = range(N7)
	x7 = np.array(vec3)
	y7 = vec3.mean(axis=0)
	error7 = x7.std(axis=0)

	fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True)
	ax4,ax5,ax6,ax7 = axes.flatten()
	for ax in [ax4,ax5,ax6,ax7]:
		for label in ax.get_xticklabels():
			label.set_rotation(90)
		Y = ['avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew']
		plt.xticks(np.arange(len(Y)), ('avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew'))
		plt.ylim([0,1])
		#plt.setp(ax.get_xticklabels(), rotation = 90)
	ax4.bar(A4, y4, yerr=error4, align='center')
	ax5.bar(A5, y5, yerr=error5, align='center')
	ax6.bar(A6, y6, yerr=error6, align='center')
	ax7.bar(A7, y7, yerr=error7, align='center')

	plt.show()

def avg_famiy_cluster(): #fraction family in each cluster
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())

	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)
	#vec = preprocessing.normalize(vec,norm='l1')

	reduced_data = PCA(n_components=2).fit_transform(vec)
	km = KMeans(4).fit(reduced_data)  #vec
	
	labels = km.labels_	
	centroids = km.cluster_centers_

	'''centers_sort = np.sort(centroids)
	clust_0 = centers_sort[0]
	clust_1 = centers_sort[1]
	clust_2 = centers_sort[2]
	clust_3 = centers_sort[3]
	vec1 = vec[np.where(labels==clust_0)]
	print vec1'''
	
	clust_0 = ClusterIndicesNumpy(0, km.labels_) # all samples in cluster number __
	clust_1 = ClusterIndicesNumpy(1, km.labels_)
	clust_2 = ClusterIndicesNumpy(2, km.labels_)
	clust_3 = ClusterIndicesNumpy(3, km.labels_)
	vec0 = vec[clust_0] #extract all of my cluster number __ data points	
	vec1 = vec[clust_1] 
	vec2 = vec[clust_2]
	vec3 = vec[clust_3]

	fig,ax = plt.subplots()

	temp0 = []
	count1 = 0
	for i in vec0.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	#a = float(count1)/float(len(vec0))
	temp0.append(count1)
	count2 = 0
	for i in vec0.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	#b = float(count2)/float(len(vec0))
	temp0.append(count2)
	count3 = 0
	for i in vec0.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	#c = float(count3)/float(len(vec0))
	temp0.append(count3)
	count4 = 0
	for i in vec0.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	#d = float(count4)/float(len(vec0))
	temp0.append(count4)
	N0 = len(np.asarray(temp0))
	A0 = range(N0)

	temp1 = []
	count1 = 0
	for i in vec1.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	#a = float(count1)/float(len(vec1))
	temp1.append(count1)
	count2 = 0
	for i in vec1.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	#b = float(count2)/float(len(vec1))
	temp1.append(count2)
	count3 = 0
	for i in vec1.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	#c = float(count3)/float(len(vec1))
	temp1.append(count3)
	count4 = 0
	for i in vec1.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	#d = float(count4)/float(len(vec1))
	temp1.append(count4)
	N1 = len(np.asarray(temp1))
	A1 = range(N1)

	temp2 = []
	count1 = 0
	for i in vec2.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	#a = float(count1)/float(len(vec2))
	temp2.append(count1)
	count2 = 0
	for i in vec2.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	#b = float(count2)/float(len(vec2))
	temp2.append(count2)
	count3 = 0
	for i in vec2.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	#c = float(count3)/float(len(vec2))
	temp2.append(count3)
	count4 = 0
	for i in vec2.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	#d = float(count4)/float(len(vec2))
	temp2.append(count4)
	N2 = len(np.asarray(temp2))
	A2 = range(N2)

	temp3 = []
	count1 = 0
	for i in vec3.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	#a = float(count1)/float(len(vec3))
	temp3.append(count1)
	count2 = 0
	for i in vec3.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	#b = float(count2)/float(len(vec3))
	temp3.append(count2)
	count3 = 0
	for i in vec3.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	#c = float(count3)/float(len(vec3))
	temp3.append(count3)
	count4 = 0
	for i in vec3.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	#d = float(count4)/float(len(vec3))
	temp3.append(count4)
	N3 = len(np.asarray(temp3))
	A3 = range(N3)
	

	fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True)
	ax0,ax1,ax2,ax3 = axes.flatten()
	plt.xticks(np.arange(4),('GPCR','Ion channel','Kinase','Enzyme'))
	plt.ylim([0,80])
	ax0.bar(A0, np.asarray(temp0), align='center')
	ax0.set_title('Avg. protein structure for cluster 0')
	ax1.bar(A1, np.asarray(temp1), align='center') 
	ax1.set_title('Avg. protein structure for cluster 1')
	ax2.bar(A2, np.asarray(temp2), align='center') 
	ax2.set_title('Avg. protein structure for cluster 2')
	ax3.bar(A3, np.asarray(temp3), align='center') 
	ax3.set_title('Avg. protein structure for cluster 3')

	plt.show()

def avg_famfeat_cluster(): #fraction family in each cluster
	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())

	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)
	#vec = preprocessing.normalize(vec,norm='l1')

	reduced_data = PCA(n_components=2).fit_transform(vec)
	x = reduced_data[:,0]
	x = (x-np.min(x))
	x = (x/np.max(x))
	#print reduced_data.shape()
	#print x.shape()
	reduced_data[:,0] = x

	x = reduced_data[:,1]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,1] = x

	'''x = reduced_data[:,2]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,2] = x'''

	km = KMeans(4).fit(reduced_data) #vec
	#labels = np.argsort(km.labels_)
	#print labels	
	centroids = km.cluster_centers_

	'''centers_sort = np.sort(centroids)
	clust_0 = centers_sort[0]
	clust_1 = centers_sort[1]
	clust_2 = centers_sort[2]
	clust_3 = centers_sort[3]

	vec0 = []
	if len(clust_0) >200:
		vec0.append(vec[clust_0])
	elif len(clust_1) >200:
		vec0.append(vec[clust_1])
	elif len(clust_2) >200:
		vec0.append(vec[clust_2])
	else:
		vec0.append(vec[clust_3])
	vec1=[]
	if len(clust_0) >40:
		vec1.append(vec[clust_1])
	elif len(clust_1) >40:
		vec1.append(vec[clust_1])
	elif len(clust_2) >40:
		vec1.append(vec[clust_2])
	else:
		vec1.append(vec[clust_3])
	vec2=[]
	if len(clust_0)>20:
		vec2.append(vec[clust_0])
	elif len(clust_1) >20:
		vec2.append(vec[clust_1])
	elif len(clust_2) >20:
		vec2.append(vec[clust_2])
	else:
		vec2.append(vec[clust_3])
	vec3=[]
	if len(clust_0)>2:
		vec3.append(vec[clust_0])
	elif len(clust_1) >2:
		vec3.append(vec[clust_1])
	elif len(clust_2) >2:
		vec3.append(vec[clust_2])
	else:
		vec3.append(vec[clust_3])'''
	
	clust_0 = ClusterIndicesNumpy(0, km.labels_) 
	clust_1 = ClusterIndicesNumpy(1, km.labels_)
	clust_2 = ClusterIndicesNumpy(2, km.labels_)
	clust_3 = ClusterIndicesNumpy(3, km.labels_)
	
	vec0 = vec[clust_0] 
	vec1 = vec[clust_1] 
	vec2 = vec[clust_2]
	vec3 = vec[clust_3]
	print len(vec0)
	print len(vec1)
	print len(vec2)
	print len(vec3)

	fig,ax = plt.subplots()

	temp0 = []
	count1 = 0
	for i in vec0.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	a = float(count1)/float(len(vec0))
	temp0.append(a)
	count2 = 0
	for i in vec0.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	b = float(count2)/float(len(vec0))
	temp0.append(b)
	count3 = 0
	for i in vec0.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	c = float(count3)/float(len(vec0))
	temp0.append(c)
	count4 = 0
	for i in vec0.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	d = float(count4)/float(len(vec0))
	temp0.append(d)
	N0 = len(np.asarray(temp0))
	A0 = range(N0)

	temp1 = []
	count1 = 0
	for i in vec1.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	a = float(count1)/float(len(vec1))
	temp1.append(a)
	count2 = 0
	for i in vec1.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	b = float(count2)/float(len(vec1))
	temp1.append(b)
	count3 = 0
	for i in vec1.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	c = float(count3)/float(len(vec1))
	temp1.append(c)
	count4 = 0
	for i in vec1.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	d = float(count4)/float(len(vec1))
	temp1.append(d)
	N1 = len(np.asarray(temp1))
	A1 = range(N1)

	temp2 = []
	count1 = 0
	for i in vec2.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	a = float(count1)/float(len(vec2))
	temp2.append(a)
	count2 = 0
	for i in vec2.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	b = float(count2)/float(len(vec2))
	temp2.append(b)
	count3 = 0
	for i in vec2.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	c = float(count3)/float(len(vec2))
	temp2.append(c)
	count4 = 0
	for i in vec2.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	d = float(count4)/float(len(vec2))
	temp2.append(d)
	N2 = len(np.asarray(temp2))
	A2 = range(N2)

	temp3 = []
	count1 = 0
	for i in vec3.tolist():
		for j in gpcr.tolist():
			if i == j:
				count1+=1 
	a = float(count1)/float(len(vec3))
	temp3.append(a)
	count2 = 0
	for i in vec3.tolist():
		for j in ic.tolist():
			if i == j:
				count2+=1 
	b = float(count2)/float(len(vec3))
	temp3.append(b)
	count3 = 0
	for i in vec3.tolist():
		for j in kin.tolist():
			if i == j:
				count3+=1 
	c = float(count3)/float(len(vec3))
	temp3.append(c)
	count4 = 0
	for i in vec3.tolist():
		for j in enz.tolist():
			if i == j:
				count4+=1 
	d = float(count4)/float(len(vec3))
	temp3.append(d)
	N3 = len(np.asarray(temp3))
	A3 = range(N3)

	fig1, axes1 = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True)
	ax0,ax1,ax2,ax3 = axes1.flatten()
	plt.xticks(np.arange(4),('GPCR','Ion channel','Kinase','Enzyme'))
	plt.ylim([0,1])

	major_ticks = np.arange(0, 1, 0.1)
	ax0.set_yticks(major_ticks)
	ax0.grid(which='major',alpha=0.2)
	ax1.set_yticks(major_ticks)
	ax1.grid(which='major',alpha=0.2)
	ax2.set_yticks(major_ticks)
	ax2.grid(which='major',alpha=0.2)
	ax3.set_yticks(major_ticks)
	ax3.grid(which='major',alpha=0.2)

	ax0.bar(A0, np.asarray(temp0), align='center')
	ax0.set_title('Avg. protein structure for cluster 0')
	ax1.bar(A1, np.asarray(temp1), align='center') 
	ax1.set_title('Avg. protein structure for cluster 1')
	ax2.bar(A2, np.asarray(temp2), align='center') 
	ax2.set_title('Avg. protein structure for cluster 2')
	ax3.bar(A3, np.asarray(temp3), align='center') 
	ax3.set_title('Avg. protein structure for cluster 3')
	
	
	vec0 = preprocessing.normalize(vec0,norm='l1')
	vec1 = preprocessing.normalize(vec1,norm='l1')
	vec2 = preprocessing.normalize(vec2,norm='l1')
	vec3 = preprocessing.normalize(vec3,norm='l1')

	for item in vec0:
		width = 0.5
		N4 = len(item)
		A4 = range(N4)
	x4 = np.array(vec0)
	y4 = vec0.mean(axis=0)
	error4 = x4.std(axis=0)
	for item in vec1:
		width = 0.5
		N5 = len(item)
		A5 = range(N5)
	x5 = np.array(vec1)
	y5 = vec1.mean(axis=0)
	error5 = x5.std(axis=0)
	for item in vec2:
		width = 0.5
		N6 = len(item)
		A6 = range(N6)
	x6 = np.array(vec2)
	y6 = vec2.mean(axis=0)
	error6 = x6.std(axis=0)
	for item in vec3:
		width = 0.5
		N7 = len(item)
		A7 = range(N7)
	x7 = np.array(vec3)
	y7 = vec3.mean(axis=0)
	error7 = x7.std(axis=0)

	fig2, axes2 = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True)
	ax4,ax5,ax6,ax7 = axes2.flatten()
	for ax in [ax4,ax5,ax6,ax7]:
		for label in ax.get_xticklabels():
			label.set_rotation(90)
		Y = ['avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew']
		plt.xticks(np.arange(len(Y)), ('avg.sp','d','r','CC','DAC','GRC','Lap','A','B','avg.deg.H', 'CC_H', 'GRC_H', 'avg.deg.C', 'CC_C', 'GRC_C','var','skew'))
		plt.ylim([0,0.07])
		
	ax4.bar(A4, y4, yerr=error4, align='center')
	ax4.set_title('Avg. feature for cluster 0')
	ax5.bar(A5, y5, yerr=error5, align='center')
	ax5.set_title('Avg. feature for cluster 1')
	ax6.bar(A6, y6, yerr=error6, align='center')
	ax6.set_title('Avg. feature for cluster 2')
	ax7.bar(A7, y7, yerr=error7, align='center')
	ax7.set_title('Avg. feature for cluster 3')

	fig1.show()
	fig2.show()

	plt.show()

def cluster_3D():

	base_list = name_base(path_pdb,file)

	gpcr = np.asarray(create_gpcr_list())
	ic = np.asarray(create_ic_list())
	enz = np.asarray(create_enz_list())
	kin = np.asarray(create_kin_list())

	vec = create_vlist()
	vec = np.asarray(vec)

	fam_list = []
	for item in vec.tolist():
		for a in gpcr.tolist():
			if item == a:
				fam_list.append(0) #purple
		for b in ic.tolist():
			if item == b:
				fam_list.append(1) #blue
		for c in kin.tolist():
			if item == c:
				fam_list.append(2) #green
		for d in enz.tolist():
			if item == d:
				fam_list.append(3) #yellow

	fam_array = np.asarray(fam_list) #cluster label

	vec = np.asarray(vec)
	vec = np.matrix(vec)
	#vec = preprocessing.normalize(vec,norm='l1')

	#PCA-reduced data
	pca = PCA(n_components=3).fit(vec)
	reduced_data = pca.transform(vec) #n_components=3
	print('Exp var: ' +str(np.sum(pca.explained_variance_ratio_)))
	x = reduced_data[:,0]
	x = (x-np.min(x))
	x = (x/np.max(x))

	reduced_data[:,0] = x

	x = reduced_data[:,1]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,1] = x

	x = reduced_data[:,2]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,2] = x

	model = KMeans(4).fit(reduced_data) 
	y_pred = model.predict(reduced_data)
	clust_labels = model.labels_
	centroids = model.cluster_centers_

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(reduced_data[:,0], reduced_data[:,1], reduced_data[:,2], c=fam_array) 
	#ax.scatter(reduced_data[:,1], reduced_data[:,2],c=fam_array) 
	#ax.scatter(reduced_data[:,0], reduced_data[:,1],c=y_pred) 


	#ax.scatter(centroids[:,0], centroids[:,1],marker='x',s=80, linewidths=3,color='r')
	#ax.scatter(centroids[:,0], centroids[:,1],marker='x',s=80, linewidths=3,color='r')

	plt.show()

	


if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_vlist()
	#create_barplots()
	#avg_barplot()
	#avg_byfeat_barplot()
	cluster_norm()
	#spectral_cluster()
	#avg_feat_cluster()
	#avg_famiy_cluster()
	#avg_famfeat_cluster()
	#cluster_3D()



	
