from __future__ import division
#Kernel-density estimate using Gaussian kernels

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
import sklearn.datasets as ds
from sklearn.decomposition import PCA
import seaborn as sns
import random
from sklearn.model_selection import train_test_split
from sklearn.neighbors.kde import KernelDensity
from sklearn.grid_search import GridSearchCV
import math 
from mayavi import mlab
from sklearn.preprocessing import StandardScaler
from sklearn.datasets.base import Bunch

sns.set(color_codes=True)

path_fv = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/feature_vectors/'
path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_edited/'
path_gpcr = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_gpcr/'
path_ic = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_ionchannel/'
path_enz = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_enzyme/'
path_kin = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/fv_kinase/'


#----------------------------------------------------------------------
#Necessary data
def name_base(path_pdb,file):
	base_list = []
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			base_list.append(base)
	return base_list

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
#----------------------------------------------------------------------
#PCA weight vector
base_list = name_base(path_pdb,file)
vec = create_vlist()
vec = np.asarray(vec)
vec = np.matrix(vec)

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

fam_array = np.asarray(fam_list)
#print fam_array

n_comp = 2
pca = PCA(n_components=n_comp)

data = pd.DataFrame(data=vec, index=range(vec.shape[0]), columns=range(vec.shape[1]))
z_scaler = StandardScaler() #Normalize
z_data = z_scaler.fit_transform(data)
pca_data = pca.fit_transform(z_data)

reduced_data = pca.fit_transform(vec) 
print pca.explained_variance_ratio_ #Amount of variance: PC1 explains 99% and PC2 explains 0.0064%
print (abs(pca.components_))

def pca_weight(reduced_data,coeff,labels=None):
	x = reduced_data[:,0]
	y = reduced_data[:,1]
	n = coeff.shape[0]
	#print n
	scalex = 1.0/(x.max() - x.min())
	scaley = 1.0/(y.max() - y.min())
	colors = ['r','g','y','b']
	plt.scatter(x * scalex,y * scaley, color=[colors[lc] for lc in fam_array])
	for i in range(n):
		plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
		if labels is None:
			plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Feat"+str(i+1), color = 'g', ha = 'center', va = 'center')
		else:
			plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')
		
plt.xlim(-0.5,1.5)
plt.ylim(-0.5,1)
plt.xlabel("PC{}".format(1))
plt.ylabel("PC{}".format(2))
plt.grid()

pca_weight(reduced_data[:,0:2],np.transpose(pca.components_[0:2, :]))
plt.show()

#----------------------------------------------------------------------
#Inverse PCA and biases
def heatmap_pca():
	n_comp =2
	pca = PCA(n_components=n_comp)
	data = pd.DataFrame(data=vec, index=range(vec.shape[0]), columns=range(vec.shape[1]))

	z_scaler = StandardScaler() #Normalize
	z_data = z_scaler.fit_transform(data)

	pca_data = pca.fit_transform(z_data)
	pca_inv_data = pca.inverse_transform(np.eye(n_comp))

	fig = plt.figure(figsize=(10, 6.5))
	sns.heatmap(pca.inverse_transform(np.eye(n_comp)), cmap="hot", cbar=True, vmin=0, vmax=1)
	'''ax = sns.heatmap(pca.inverse_transform(np.eye(n_comp)), cbar=False, vmin=0, vmax=1)
	cbar = ax.figure.colorbar(ax.collections[0])
	cbar.set_ticks([0,1])'''
	plt.ylabel('principal component', fontsize=20)
	plt.xlabel('original feature index', fontsize=20)
	#plt.xlim(1,17)
	#plt.tick_params(axis='both', which='major', labelsize=18)
	#plt.tick_params(axis='both', which='minor', labelsize=12)

	fig = plt.figure(figsize=(10, 6.5))
	plt.plot(pca_inv_data.mean(axis=0), '--o', label = 'mean')
	plt.plot(np.square(pca_inv_data.std(axis=0)), '--o', label = 'variance')
	plt.legend(loc='lower right')
	plt.ylabel('feature contribution', fontsize=20)
	plt.xlabel('feature index', fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18)
	plt.tick_params(axis='both', which='minor', labelsize=12)
	plt.xlim([0, 16])
	plt.legend(loc='lower left', fontsize=18)
	plt.show()
#----------------------------------------------------------------------
#Two-dimensional data KDE
def gauss_kde():

	pca = PCA(n_components=2)
	reduced_data = pca.fit_transform(vec) 

	x = reduced_data[:,0]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,0] = x

	x = reduced_data[:,1]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,1] = x

	'''x = reduced_data[:,2]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,2] = x'''

	x = reduced_data[:,0]
	y = reduced_data[:,1]
	#z = reduced_data[:,2]
	xmin, xmax = x.min(), x.max()
	ymin, ymax = y.min(),y.max()
	#zmin, ymax = y.min(),y.max()
	#Kernel density estimate on data using h = Scott's rule (bandwidth)
	bw = 0.21
	X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
	positions = np.vstack([X.ravel(), Y.ravel()])
	values = np.vstack([x, y])
	kernel = stats.gaussian_kde(values, bw_method = bw)
	Z = np.reshape(kernel(positions).T, X.shape)

	#Plot the results
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
		  extent=[xmin, xmax, ymin, ymax])
	ax.plot(x, y, 'k.', markersize=2)
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([ymin, ymax])

	'''fig = plt.figure()
	ax = fig.gca()
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	# Contourf plot
	cfset = ax.contourf(X, Y, Z, cmap='Blues')
	## Or kernel density estimate plot instead of the contourf plot
	#ax.imshow(np.rot90(Z), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
	# Contour plot
	cset = ax.contour(X, Y, Z, colors='k')
	# Label plot
	ax.clabel(cset, inline=1, fontsize=10)
	#ax.set_xlabel('Y1')
	#ax.set_ylabel('Y0')'''


	plt.show()

#----------------------------------------------------------------------
#Test for sigma
def sigma_test():
	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)

	pca = PCA(n_components=3).fit(vec)
	reduced_data = pca.transform(vec)
	
	random.shuffle(reduced_data)
	X_train, X_val = train_test_split(reduced_data, test_size=0.50)
	
	# Max log likelihood
	n = len(X_train)
	#fin = []
	temp = []
	#for y in range(10):
	for sig in np.arange(0.01,2,0.1): #increments of 0.25
		
		for item_val in X_val:
			#print item_val
			for item_train in X_train: 
				#print item_train
				d_xval2 = (np.exp(-((np.power((item_val-item_train),2))/(2*np.power(sig,2)))))
		d_xval = (1/(n*np.sqrt(2*np.pi*sig)))*(np.sum(d_xval2))
		log_d = np.log(d_xval)
		#print log_d
		sum_log = np.sum(log_d)
		temp.append([sig, sum_log])
	#print temp
		#print 'log likelihood using',str(sig), '=', str(sum_log) #sum_log
	opt_sig = max(temp, key=lambda x:x[1])
	print opt_sig #return

def avg_sigma():
	from itertools import repeat, starmap
	sigma = sigma_test()
	fin = []
	for i in repeat(10):
		fin.append(sigma)
	print fin
	#median = np.median(fin)
	#print median

#----------------------------------------------------------------------
#Histogram visual - NOT finished
def hist_plot():
	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)

	pca = PCA(n_components=2).fit(vec)
	reduced_data = pca.transform(vec) 

	x = reduced_data[:,0]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,0] = x

	x = reduced_data[:,1]
	x = (x-np.min(x))
	x = (x/np.max(x))
	reduced_data[:,1] = x

	x = reduced_data[:,0]
	y = reduced_data[:,1]

	fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(21,5))
	nbins=24

	#2D Histogram 
	axes[0].set_title('2D Histogram')
	axes[0].hist2d(x,y,bins=nbins, cmap=plt.cm.jet) # bins = square root choice

	#Gaussian KDE
	bw = 1.21
	xy = np.vstack([x,y])
	kde = stats.gaussian_kde(xy, bw_method = bw)
	xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
	zi = kde(np.vstack([xi.flatten(), yi.flatten()]))
	axes[1].set_title('Gaussian KDE')
	axes[1].pcolormesh(xi,yi,zi.reshape(xi.shape), cmap=plt.cm.jet)

	#2D Density with shading
	axes[2].set_title('2D Density with shading')
	axes[2].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.jet)
	 
	#2D Density with shading and contour
	axes[3].set_title('2D Density with shading + contour')
	axes[3].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.jet)
	axes[3].contour(xi, yi, zi.reshape(xi.shape) )
	
	plt.show()


def sns_plot():
	base_list = name_base(path_pdb,file)
	vec = create_vlist()
	vec = np.asarray(vec)
	vec = np.matrix(vec)

	pca = PCA(n_components=3).fit(vec)
	reduced_data = pca.transform(vec) 

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

	x = reduced_data[:,0]
	y = reduced_data[:,1]
	z = reduced_data[:,2]
	xmin, xmax = x.min(), x.max()
	ymin, ymax = y.min(),y.max()
	zmin, zmax = z.min(),z.max()

	#Density scatter plot
	bw = 1.09
	xyz = np.vstack([x,y,z])
	kde = stats.gaussian_kde(xyz, bw_method = bw)
	density = kde(xyz)

	figure = mlab.figure('DensityPlot')
	pts = mlab.points3d(x, y, z, density, scale_mode='none', scale_factor=0.07)

	#Evaluate the gaussian kde on a grid
	'''bw = 1.92
	X, Y, Z = np.mgrid[xmin:xmax:100j, ymin:ymax:100j, zmin:zmax:100j]
	values = np.vstack([x, y, z])
	kde = stats.gaussian_kde(values, bw_method = bw)
	coords = np.vstack([item.ravel() for item in [X,Y,Z]])
	density = kde(coords).reshape(X.shape)

	figure = mlab.figure('DensityPlot')

	grid = mlab.pipeline.scalar_field(X, Y, Z, density)
	min = density.min()
	max=density.max()
	mlab.pipeline.volume(grid, vmin=min, vmax=min + .5*(max-min))'''

	mlab.axes()
	mlab.show()

if __name__ == '__main__':
	#pca_weight()
	#heatmap_pca()
	gauss_kde()
	#sigma_test()
	#avg_sigma()
	#hist_plot()
	#sns_plot()






























