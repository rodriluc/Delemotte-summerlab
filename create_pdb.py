#script creates and saves PDB file containing hydrophobic/ charged residues

path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/'
import mdtraj as md
import os

def name_base():
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'):
			base = file[:7]
			return base

def create_hydrophobic():
	base = name_base()
	
	for file in os.listdir(path_pdb):
		if file == base+'.pdb': 
			iter_file = open(file)
			lines = iter_file.readlines()
			head = lines[0]
			top = md.load_pdb(file).topology
			with open ('hydrophobic_'+base+'.pdb', 'w') as w:
				w.write(head)
				hydrophobic_list = top.select('resname PHE MET TRP ILE VAL LEU PRO ALA')
				new_list = [x+1 for x in hydrophobic_list]
				T = [lines[i] for i in new_list] #new_list
				for line in T:
					if line.rstrip():
						w.write(line) 
				w.write('END')
				w.write('\n')
	print 'Hydrophobic file saved!'
				
def create_charged():
	base = name_base()
	
	for file in os.listdir(path_pdb):
		#print file
		if file == base+'.pdb': 
			iter_file = open(file)
			lines = iter_file.readlines()
			head = lines[0]
			top = md.load_pdb(file).topology
			with open ('charged_'+base+'.pdb', 'w') as w:
				w.write(head)
				charged_list = top.select('resname ARG LYS ASP GLU HIS') #do i inculde histidine
				new_list = [x+1 for x in charged_list]
				T = [lines[i] for i in new_list] 
				for line in T:
					if line.rstrip():
						w.write(line) 
				w.write('END')
				w.write('\n')
	print 'Charged file saved!'


if __name__ == '__main__':
	name_base()
	create_hydrophobic()
	create_charged()

