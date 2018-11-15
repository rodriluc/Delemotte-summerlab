#script creates and saves PDB file containing hydrophobic/ charged residues

path_pdb = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/PDB_edited/'
path_hydro = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/hydrophobic_files/'
path_charg = '/data2/LucieR/Delemotte-summerlab_ERnorm_100ER_4.5A/charged_files/'
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

		iter_file = open((path_pdb+file), 'r')
		lines = iter_file.readlines()
		head = lines[0]
		top = md.load_pdb(path_pdb+file).topology
		with open ((path_hydro+('hydrophobic_'+file)), 'w') as w:
			w.write(head)
			hydrophobic_list = top.select('resname PHE or resname APHE or resname BPHE or resname CPHE or resname MET or resname AMET or resname BMET or resname TRP or resname ATRP or resname BTRP or resname ILE or resname AILE or resname BILE or resname CILE or resname VAL or resname AVAL or resname BVAL or resname CVAL or resname LEU or resname ALEU or resname BLEU or resname CLEU or resname PRO or resname APRO or resname BPRO or resname CPRO or resname ALA or resname AALA or resname BALA')
			#print hydrophobic_list
			new_list = [x for x in hydrophobic_list]
			T = [lines[i+1] for i in hydrophobic_list] #new_list
			for line in T: #T
				if line.rstrip():
					w.write(line) 
			w.write('END')
			w.write('\n')
			print 'Hydrophobic file saved!'
				
def create_charged():
	base = name_base()
	
	for file in os.listdir(path_pdb):

		iter_file = open((path_pdb+file), 'r')
		lines = iter_file.readlines()
		head = lines[0]
		top = md.load_pdb(path_pdb+file).topology
		with open ((path_charg+('charged_'+file)), 'w') as w:
			w.write(head)
			charged_list = top.select('resname ARG or resname AARG or resname BARG or resname CARG or resname LYS or resname ALYS or resname BLYS or resname CLYS or resname ASP or resname AASP or resname BASP or resname CASP or resname GLU or resname AGLU or resname BGLU or resname CGLU or resname HIS or resname AHIS or resname BHIS or resname CHIS') 
			new_list = [x for x in charged_list]
			T = [lines[i+1] for i in charged_list] 
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


