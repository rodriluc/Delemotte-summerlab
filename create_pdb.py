#script creates and saves PDB file containing hydrophobic/ charged residues

path_pdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_edited/'
path_hydro = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/hydrophobic_files/'
path_charg = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/charged_files/'
import mdtraj as md
import os

def name_base():
	for file in os.listdir(path_pdb):
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			return base

def create_hydrophobic():
	base = name_base()
	
	for file in os.listdir(path_pdb):
		#print file
		if file == base+'.pdb': 
			#print file
			iter_file = open((path_pdb+file), 'r')
			lines = iter_file.readlines()
			head = lines[0]
			top = md.load_pdb(path_pdb+file).topology
			with open ((path_hydro+('hydrophobic_'+base+'.pdb')), 'w') as w:
				w.write(head)
				hydrophobic_list = top.select('resname PHE MET TRP ILE VAL LEU PRO ALA')
				#print hydrophobic_list
				new_list = [x+1 for x in hydrophobic_list]
				T = [lines[i] for i in hydrophobic_list] #new_list
				for line in T: #T
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
			iter_file = open((path_pdb+file), 'r')
			lines = iter_file.readlines()
			head = lines[0]
			top = md.load_pdb(path_pdb+file).topology
			with open ((path_charg+('charged_'+base+'.pdb')), 'w') as w:
				w.write(head)
				charged_list = top.select('resname ARG LYS ASP GLU HIS') #do i inculde histidine
				new_list = [x+1 for x in charged_list]
				T = [lines[i] for i in charged_list] 
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

