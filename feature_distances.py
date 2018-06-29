# Aim: to calculate feature distances
import os
import fileinput
from glob import glob
import linecache
import math

path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/'

def find_originfeat():
	
	for file in os.listdir(path_pdb):
		if file.endswith('.txt') and file.startswith('features'): 
			#print file
			iter_file = open(file)
			lines = iter_file.readlines()
			#lines = map(lambda it: it.strip(), lines)
			lines = [item.strip() for item in lines]
			avg_degree = lines[0]
			#print avg_degree
			avg_path = lines[1]
			diam = lines[2]
			rad = lines[3]
			clust_coeff = lines[4]
			lap_matrix = lines[5:55] #might need to convert to ACTUAL matrix
			ass_coeff = lines[56]
			glob_rc = lines[57]
			a_cont = lines[58]
			b_cont = lines[59]
			
			return avg_degree, avg_path, diam, rad, clust_coeff, lap_matrix, ass_coeff, glob_rc, a_cont, b_cont
			
def find_hydrofeat():
			
	for file in os.listdir(path_pdb):
		if file.endswith('.txt') and file.startswith('features'): 
			iter_file = open(file)
			lines = iter_file.readlines()
			avg_hdegree = lines[60]
			hlclust_coeff = lines [61]
			hloc_rc = lines[62]
			
			return avg_hdegree, hlclust_coeff, hloc_rc
		
def find_chargfeat():
	
	for file in os.listdir(path_pdb):
		if file.endswith('.txt') and file.startswith('features'): 
			iter_file = open(file)
			lines = iter_file.readlines()
			avg_cdegree = lines[63]
			clclust_coeff = lines [64]
			cloc_rc = lines[65]
			
			return avg_cdegree, clclust_coeff, cloc_rc
		
def common_feat(): #compiles one feature for all PDBs into another list
	#origin = find_originfeat()[0][1][2][3][4][5][6][7][8][9][10]
	#hydro = find_hydrofeat()[0][1][2]
	#charg = find_chargfeat()[0][1][2]
	new_feat = []
	for file in os.listdir(path_pdb):
		
		'''fnames = glob('features*') #feature vectors, fnames is list
		for line in fileinput.input(fnames): #line is string
			for i in line:
				print i[0]'''
		if file.startswith('features'): 
			with open(file) as f:
				f=f.readlines()[0]
				new_feat.append(f)
				#print item
			#line = linecache.getline(fnames, 2)
			#new_feat.append(line[0])
			#print new_feat


def calc_distance1():
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'): 
			for element in file:
				d = math.sqrt()
	
					
if __name__ == '__main__':
	#print(find_originfeat())
	#print(find_hydrofeat())
	#print(find_chargfeat())
	common_feat()
	#calc_distance1()