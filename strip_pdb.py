#Aim: to strip PDBs, format CRYST1 + ATOM + END
import os

path_originpdb = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_files/'
path_edit = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/PDB_edited/'

def name_base():
	for file in os.listdir(path_originpdb):
		if file.endswith('.pdb') and not file.startswith('hydrophobic') and not file.startswith('charged'):
			basename = file.split('.')[:-1]
			base =''.join(basename)
			return base

def strip_pdb():
	base = name_base()

	for file in os.listdir(path_originpdb):
		#print file
		#if file == base+'.pdb': 
		if file.endswith('.pdb'):
			#print file
			with open((path_originpdb+file), 'r') as lines:
				#print file
				#lines = iter_file.readlines()
				cryst = ('CRYST1')
				with open ((path_edit+(file)), 'w') as w:
					for line in lines:
						#print line
						if cryst in line:
							#print line
							w.write(line)
						elif line.startswith('ATOM'):
							w.write(line)
					w.write('END')
				print 'Stripped PDB ready!'



	
if __name__ == '__main__':
	strip_pdb()
