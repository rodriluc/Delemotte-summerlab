# Aim: create barplots and possibly normalize data based on random graph

import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing

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
	
	min_max_scaler = preprocessing.MinMaxScaler()
	fv_minmax = min_max_scaler.fit_transform(fv_list) #normalize data set 0-1
	#print fv_minmax

	Tot = 15
	Cols = 5
	Rows = Tot//Cols
	Rows += Tot%Cols
	Position = range(1,Tot+1)	
	
	for item in fv_minmax[:15]:
		width = 0.5
		N = len(item)
		x = range(N)

	fig = plt.figure(1)
	for i in range(Tot):

		ax=fig.add_subplot(Rows,Cols,Position[i])
		ax.bar(x, fv_minmax[i], width, align='center')
		plt.ylim([0,1])

	plt.show()


if __name__ == '__main__':
	#name_base(path_pdb,file)
	#create_vlist()
	print(create_barplots())


