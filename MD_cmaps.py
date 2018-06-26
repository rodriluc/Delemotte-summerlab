import argparse
import sys
import os
import time

python_path = os.path.dirname(__file__);

next_folder = '';
parent_folder = '';
for i in range(len(python_path)-1):
	next_folder+=python_path[i];
	if python_path[i]=='/':
		parent_folder += next_folder;
		next_folder = '';
		
sys.path.append(python_path);
sys.path.append(parent_folder);

import mdtraj as md
import numpy as np
from scipy import stats
from scipy.spatial.distance import squareform, pdist, cdist
from joblib import Parallel, delayed

path_data = '/afs/kth.se/home/l/u/lucier/Documents/protein_networks/Results_data/'

def unwrap_cmap_loop(arg,**kwarg):
	return MD_cmaps.distance_matrix_loop(*arg,**kwarg);

def unwrap_cmap_semi_bin_loop(arg,**kwarg):
	return MD_cmaps.distance_matrix_semi_bin_loop(*arg,**kwarg);

class MD_cmaps():
	
	file_end_name = '';
	save_folder = '';
	cmap = [];
	
	def __init__(self):
		return;

	def getContactMap(self, distanceMatrix, cutoff, doNormalization=False, makeBinary=False):
		# Construct a contact map from a distance matrix using a cutoff. Contact map can be made binary as well.
		nRows = len(distanceMatrix[::,0]);
		contactMap = distanceMatrix; 
		
		for i in range(0, nRows):
			tmpVec = distanceMatrix[i,::];		
			contactMap[i,::] = tmpVec * (tmpVec <= cutoff);
		
		# Normalize contact map
		if doNormalization:
			contactMap = contactMap/cutoff;
		
		# Make the contacts binary
		if makeBinary:
			contactMap = contactMap > 0;
            
		np.savetxt(self.save_folder + 'cmap_processed_' + self.file_end_name + '.txt',contactMap);
		print('Saved cmap!') #savetxt
		print(contactMap.shape)
		return;

	def getAllCalphaInverseDistances(self,traj, startID=-1,endID=-1):
		print('Compute all C_alpha inverse distances.')
		nFrames = int(traj.n_frames);
		nResidues = int(traj.n_residues);
		allInds = [];

		if startID == -1: 
			# Construct the distance matrices (nFrames-residue-residue) with the distance 
			# between residues defined as the minimum distance between C_alphas of the two residues.
			# Do atom selections, save list with all heavy atoms.
			for i in range(0,nResidues):
				# Use resid so that multichains can be analyzed also.
				query = "protein and resid " + str(i) + "and name CA";
				tmpInd = traj.topology.select(query);
				allInds.append(tmpInd);
		
		else:
			# Construct the distance matrices (nFrames-residue-residue) with the distance 
			# between residues defined as the minimum distance between C_alphas of the two residues.

			# Do atom selections, save list with all heavy atoms.
			for i in range(startID,endID):
				query = "protein and resid " + str(i) + "and name CA";
				tmpInd = traj.topology.select(query);
				allInds.append(tmpInd);
		
		nResidues = int(len(allInds));	
		
		distanceMatrices = np.zeros((nFrames,nResidues,nResidues));
		
		# Compute distance matrix
		for i in range(0,nResidues):
			for j in range(i+1,nResidues):

				# Get all atom pairs			
				atom_pairs = np.zeros((1,2));
				if len(allInds[i] != 0) and len(allInds[j] != 0):
					atom_pairs[0,0] = allInds[i];
					atom_pairs[0,1] = allInds[j];

					distances = md.compute_distances(traj, atom_pairs, periodic=False);

					if len(distances) == 0:
						print('The chosen residue does not exist!');

					# The distance between residues is min distance between all heavy atoms. 
					# Take residual to get rid of cut-off.
					minDistance = np.min(distances,axis=1);
					distanceMatrices[::,i,j] = 1/minDistance;
					distanceMatrices[::,j,i] = 1/minDistance;
		return distanceMatrices;
	
	def getAtomPairs(self, inds1, inds2):
		# Construct array with all pairs
		atom_pairs = np.zeros((len(inds1)*len(inds2),2));
		counter = 0;
		for k in range(0,len(inds1)):
			for l in range(0,len(inds2)):
				atom_pairs[counter,0] = inds1[k];
				atom_pairs[counter,1] = inds2[l];
				counter += 1;
		atom_pairs = atom_pairs[0:counter-1,::];
		
		return atom_pairs;

	def getInverseCalphaDistanceMatrix(self, traj, allInds):
		# Compute one C-alpha distance matrix.
		# - traj is a one-frame trajectory.
		
		# Construct the distance matrix [nResidues x nResidues] with the distance 
		# between residues defined as the minimum distance between C_alphas of the two residues.
		
		nResidues = int(len(allInds));	
		
		distanceMatrix = np.zeros((nResidues,nResidues));
		
		# Compute distance matrix
		for i in range(0,nResidues):
			for j in range(i+1,nResidues):
				# Get all atom pairs			
				atom_pairs = np.zeros((1,2));
				if len(allInds[i] != 0) and len(allInds[j] != 0):
					atom_pairs[0,0] = allInds[i];
					atom_pairs[0,1] = allInds[j];

					distances = md.compute_distances(traj, atom_pairs, periodic=False);

					if len(distances) == 0:
						print('The chosen residue does not exist!');

					# The distance between residues is min distance between all heavy atoms. 
					# Take residual to get rid of cut-off.
					minDistance = np.min(distances,axis=1);
					distanceMatrix[i,j] = 1/minDistance;
					distanceMatrix[j,i] = 1/minDistance;
		return distanceMatrix;

	def getAllSideChainMinDistances(self,traj, startID=-1, endID=-1):

		# Construct the distance matrices (nFrames-residue-residue) with the distance 
		# between residues defined as the minimum distance between all heavy atoms of the two residues.
		nFrames = int(traj.n_frames);
		nResidues = int(traj.n_residues);
		allInds = [];
		if startID == -1:
			# Do atom selections, save list with all heavy atoms.
			for i in range(0,nResidues):
				# OBS! query is now done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and resid " + str(i) + " and !(type H)";
				tmpInd = traj.topology.select(query);
				allInds.append(tmpInd);
		else:
			# Construct the distance matrices (nFrames-residue-residue) with the distance 
			# between residues defined as the minimum distance between heavy atoms of the two residues.

			# Do atom selections, save list with all heavy atoms.
			for i in range(startID,endID):
				# OBS! query is here done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and residue " + str(i) + "and name CA";
				tmpInd = traj.topology.select(query);
				allInds.append(tmpInd);		
		
		nResidues = int(len(allInds));	

		distanceMatrices = np.zeros((nFrames,nResidues,nResidues));

		# Compute distance matrix
		for i in range(0,nResidues):
			for j in range(i,nResidues):

				atomInd1 = allInds[i];
				atomInd2 = allInds[j];

				atom_pairs = self.getAtomPairs(atomInd1, atomInd2);
				
				distances = md.compute_distances(traj, atom_pairs, periodic=False);
				
				if len(distances) == 0:
					print('The chosen residue does not exist!');

				# The distance between residues is min distance between all heavy atoms. Take mean over all frames.
				distanceMatrices[::,i,j] = np.min(distances,axis=1);
				distanceMatrices[::,j,i] = np.min(distances,axis=1);
		return distanceMatrices;


	def computeFrameToFrameSideChainContacts(self, traj, startID, endID, query='protein'):

		atom_indices = traj.topology.select(query);
		traj.atom_slice(atom_indices, inplace=True);
		print('Compute frame to frame sidechain contact map difference');
		print(traj);
		
		# Compute frame-frame residue contact map norm.
		nFrames = int(traj.n_frames);		
		
		frame2frameContacts = np.zeros((nFrames,nFrames));
		
		distanceMatrices = self.getAllSideChainMinDistances(traj,startID,endID);
		
		print('Compute frame to frame distances');
		for i in range(0,nFrames):
			print(str(i)+'/'+str(nFrames)); 
			tmpMap1 = self.getContactMap(distanceMatrices[i,::,::],0.8);
			for j in range(i+1,nFrames):
				tmpMap2 = self.getContactMap(distanceMatrices[j,::,::],0.8);
				frame2frameContacts[i,j] = np.linalg.norm((tmpMap1-tmpMap2),2);
		
		frame2frameContacts = frame2frameContacts + frame2frameContacts.T;
		
		print(frame2frameContacts);
		distances = squareform(frame2frameContacts);
		# Save distance matrix to file
		np.savetxt(self.save_folder + 'frame_to_frame_side_chain_contacts_' + self.file_end_name + '.txt',distances);
		return;

	
	def computeFrameToFrameCalpaContactsMemory(self, traj, query='protein'):
		
		atom_indices = traj.topology.select(query);
		traj.atom_slice(atom_indices, inplace=True);
		print('Compute frame to frame C_alpha map difference');
		print('Atom query: ' + query);
		print(traj);
		
		# Compute frame-frame residue contact map norm.
		nFrames = int(traj.n_frames);
		
		frame2frameContacts = np.zeros((nFrames,nFrames));
		
		allInds = [];
		# Do atom selections, save list with all heavy atoms.
		for i in range(0,int(traj.n_residues)):
			# Use resid so that multichains can be analyzed also.
			query = "protein and name CA and resid " + str(i);
			tmpInd = traj.topology.select(query);
			allInds.append(tmpInd);
		
		print('Compute frame to frame distances');
		for i in range(0,nFrames):
			print(str(i+1)+'/'+str(nFrames)); 
			tmpMap1 = self.getInverseCalphaDistanceMatrix(traj[i],allInds)
			for j in range(i+1,nFrames):
				print(str(j+1)+'/'+str(nFrames));
				tmpMap2 = self.getInverseCalphaDistanceMatrix(traj[j],allInds)
				frame2frameContacts[i,j] = np.linalg.norm((tmpMap1-tmpMap2),2);
		
		frame2frameContacts = frame2frameContacts + frame2frameContacts.T;
		
		print(frame2frameContacts);
		distances = squareform(frame2frameContacts);
		# Save distance matrix to file
		np.savetxt(self.save_folder + 'frame_to_frame_CA_contacts_' + self.file_end_name + '.txt',distances);
		return;
	
	
	def computeFrameToFrameCalphaContacts(self, traj, query='protein'):
		
		atom_indices = traj.topology.select(query);
		traj.atom_slice(atom_indices, inplace=True);
		print('Compute frame to frame C_alpha map difference');
		print('Atom query: ' + query);
		print(traj);
		
		# Compute frame-frame residue contact map norm.
		nFrames = int(traj.n_frames);		
		
		frame2frameContacts = np.zeros((nFrames,nFrames));
		
		distanceMatrices = self.getAllCalphaInverseDistances(traj);
		
		print('Compute frame to frame distances');
		for i in range(0,nFrames):
			print(str(i)+'/'+str(nFrames)); 
			tmpMap1 = distanceMatrices[i,::,::];
			for j in range(i+1,nFrames):
				tmpMap2 = distanceMatrices[j,::,::];
				frame2frameContacts[i,j] = np.linalg.norm((tmpMap1-tmpMap2),2);
		
		frame2frameContacts = frame2frameContacts + frame2frameContacts.T;
		
		print(frame2frameContacts);
		distances = squareform(frame2frameContacts);
		# Save distance matrix to file
		np.savetxt(self.save_folder + 'frame_to_frame_CA_contacts_' + self.file_end_name + '.txt',distances);
		return;
	

	def distance_matrix_loop(self,i):
		print(str(i+1)+'/'+str(self.nResidues));
		for j in range(i+1,self.nResidues): # i before

			atomInd1 = self.allInds[i];
			atomInd2 = self.allInds[j];
			
			atom_pairs = self.getAtomPairs(atomInd1, atomInd2);
			
			distances = md.compute_distances(self.traj, atom_pairs, periodic=False);
			#may have to compute distances myself if pdb_load error
            
			#print distances;
			#print '---------------------------------------------------'

			if len(distances) == 0:
				print('The chosen residue does not exist!');

			# The distance between residues is min distance between all heavy atoms. Take mean over all frames.
			
			self.distanceMatrix[i,j] = np.mean(np.min(distances,axis=1),axis=0);
			self.distanceMatrix[j,i] = np.mean(np.min(distances,axis=1),axis=0);
			
		return;

	def distance_matrix_semi_bin_loop(self,i):
		# Compute a semi-binary contact map. Residue pair within the cutoff (5 angstrom) is a contact. Outside, the "degree" of contact decreases with a gaussian (for smoothness).
		print(str(i+1)+'/'+str(self.nResidues));
		
		std_dev = 0.1667; # 1.667 angstrom standard deviation => gives 1e-5 weight at 0.8 nm.
		cutoff = 0.45 # within 4.5 angstrom, the weight is 1.
		
		# Compute normalizing factor
		cutoff_value = np.exp(-cutoff**2/(2*std_dev**2))
		
		for j in range(i+1,self.nResidues):
			
			atomInd1 = self.allInds[i];
			atomInd2 = self.allInds[j];
			
			atom_pairs = self.getAtomPairs(atomInd1, atomInd2);
			
			distances = md.compute_distances(self.traj, atom_pairs, periodic=False);
			minDistances =	np.min(distances,axis=1);		
			gaussians = np.exp(-minDistances**2/(2*std_dev**2))/cutoff_value;
			gaussians[minDistances < cutoff] = 1.0
			
			if len(distances) == 0:
				print('The chosen residue does not exist!');
			
			# The distance between residues is min distance between all heavy atoms. Take mean over all frames.
			self.distanceMatrix[i,j] = np.mean(gaussians,axis=0);
			self.distanceMatrix[j,i] = np.mean(gaussians,axis=0);
			
		return;
	
	def computeAverageSideChainMinDistanceMap(self, startID, endID, query='protein'):
		# Construct the average distance matrix (residue-residue) with the distance between residues defined as the minimum distance between all heavy atoms of the two residues.
		atom_indices = self.traj.topology.select(query);
		self.traj.atom_slice(atom_indices, inplace=True);
		print('Compute average sidechain contact map');
		print('Atom query: ' + query);
		print(self.traj);
		
		nFrames = int(self.traj.n_frames);
		self.allInds = [];
		if startID == -1:
			# Do atom selections, save list with all heavy atoms.
			for i in range(0,self.nResidues):
				# OBS! query is now done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and !(type H) and resid " + str(i);
				tmpInd = self.traj.topology.select(query);
				self.allInds.append(tmpInd);
		else:
			# Construct the distance matrices (nFrames-residue-residue) with the distance 
			# between residues defined as the minimum distance between heavy atoms of the two residues.

			# Do atom selections, save list with all heavy atoms.
			for i in range(startID,endID):
				# OBS! query is here done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and !(type H) and residue " + str(i);
				tmpInd = self.traj.topology.select(query);
				self.allInds.append(tmpInd);		
		
		self.nResidues = int(len(self.allInds));	

		self.distanceMatrix = np.zeros((self.nResidues,self.nResidues));
		#Parallel(n_jobs=-1, backend="threading")(delayed(unwrap_cmap_loop)(i) for i in zip([self]*self.nResidues,range(self.nResidues)));
		
		# Compute distance matrix
		
		for i in range(0,self.nResidues):
			print(str(i+1)+'/'+str(self.nResidues));
			self.distance_matrix_loop(i)
			'''for j in range(i,self.nResidues):

				atomInd1 = allInds[i];
				atomInd2 = allInds[j];
				
				atom_pairs = self.getAtomPairs(atomInd1, atomInd2);
				#print atom_pairs;
				distances = md.compute_distances(traj, atom_pairs, periodic=False);
				
				#print distances;
				#print '---------------------------------------------------'

				if len(distances) == 0:
					print('The chosen residue does not exist!');

				# The distance between residues is min distance between all heavy atoms. Take mean over all frames.
				distanceMatrix[i,j] = np.mean(np.min(distances,axis=1),axis=0);
				distanceMatrix[j,i] = np.mean(np.min(distances,axis=1),axis=0);
				#print distanceMatrix[i,j];'''
		#distanceMatrix = np.mean(self.getAllSideChainMinDistances(traj,startID,endID),axis=0);
		
		print(self.distanceMatrix);
		self.cmap = self.distanceMatrix;

		# Save distance matrix to file
		np.savetxt(self.save_folder + 'distance_matrix_min_' + self.file_end_name + '.txt',squareform(self.distanceMatrix));
		#print(self.save_folder + 'distance_matrix_min_' + self.file_end_name + '.txt')
		print('Data saved to file!');
		return;

	def computeAverageSideChainSemiBinCmap(self, startID, endID, query='protein'):
		# Construct the average distance matrix (residue-residue) with the distance between residues defined as the minimum distance between all heavy atoms of the two residues.
		atom_indices = self.traj.topology.select(query);
		self.traj.atom_slice(atom_indices, inplace=True);
		print('Compute average semi binary sidechain contact map');
		print('Atom query: ' + query);
		print(self.traj);
		
		nFrames = int(self.traj.n_frames);
		self.allInds = [];
		if startID == -1:
			# Do atom selections, save list with all heavy atoms.
			for i in range(0,self.nResidues):
				# OBS! query is now done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and !(type H) and resid " + str(i);
				tmpInd = self.traj.topology.select(query);
				self.allInds.append(tmpInd);
		else:
			# Construct the distance matrices (nFrames-residue-residue) with the distance 
			# between residues defined as the minimum distance between heavy atoms of the two residues.

			# Do atom selections, save list with all heavy atoms.
			for i in range(startID,endID):
				# OBS! query is here done on "residue". Can be a problem for multi-chain proteins. 
				# Then, switch to resid or feed chains separately.
				query = "protein and !(type H) and resid " + str(i);
				tmpInd = self.traj.topology.select(query);
				self.allInds.append(tmpInd);		
		
		self.nResidues = int(len(self.allInds));	
		
		self.distanceMatrix = np.zeros((self.nResidues,self.nResidues));
		
		Parallel(n_jobs=28, backend="threading")(delayed(unwrap_cmap_semi_bin_loop)(i) for i in zip([self]*self.nResidues,range(self.nResidues)));
		
		print(self.distanceMatrix);
		self.cmap = self.distanceMatrix;
		
		# Save distance matrix to file
		np.savetxt(self.save_folder + 'distance_matrix_semi_bin_' + self.file_end_name + '.txt',squareform(self.distanceMatrix));

		print('Data saved to file!');
		return;
	
	
	def computeCalpaCmapDistanceToFrame1(self, traj, query='protein',do_one_resid=False,resid=0):
		
		atom_indices = traj.topology.select(query);
		traj.atom_slice(atom_indices, inplace=True);
		print('Compute frame to frame C_alpha map difference');
		print('Atom query: ' + query);
		print(traj);
		
		# Compute frame-frame residue contact map norm.
		nFrames = int(traj.n_frames);		
		
		distance_to_frame1 = np.zeros(nFrames);
		
		distanceMatrices = self.getAllCalphaDistances(traj);


		if do_one_resid:
			distanceMatrices = distanceMatrices[::,::,resid];
		
		if do_one_resid:
			reference_map = distanceMatrices[0,::];
		else:
			reference_map = distanceMatrices[0,::,::];
		print('Compute cmap distance to frame 1');
		for i in range(1,nFrames):
			if np.mod(i,10)==0:
				print(str(i)+'/'+str(nFrames)); 
			if do_one_resid: 
				tmpMap1 = distanceMatrices[i,::];
			else:
				tmpMap1 = distanceMatrices[i,::,::];
			
			distance_to_frame1[i] = np.linalg.norm((tmpMap1-reference_map),2);
		
		print(distance_to_frame1);
		
		# Save distance matrix to file
		if do_one_resid:
			np.savetxt(self.save_folder + 'cmap_CA_distance_to_frame1_resid' +str(resid) + self.file_end_name + '.txt',distance_to_frame1);
		else:
			np.savetxt(self.save_folder + 'cmap_CA_distance_to_frame1_' + self.file_end_name + '.txt',distance_to_frame1);
		
		return;
	
	def main(self,parser):

		parser.add_argument('-sc_ff','--frame_frame_side_chain_cmap',help='Flag for computing frame-frame side-chain distance map (optional).',action='store_true');
		parser.add_argument('-sc_cmap','--side_chain_cmap',help='Flag for computing average side-chain distance map (optional).',action='store_true');
		parser.add_argument('-sc_cmap_semi_bin','--side_chain_cmap_semi_binary',help='Flag for computing average side-chain distance map with binary/Gaussian kernels with standard deviation 5 Angstrom. (optional).',action='store_true');
		parser.add_argument('-ca_ff','--frame_frame_CA_cmap',help='Flag for computing frame-frame C_alpha distance map (optional)',action='store_true');
		parser.add_argument('-ca_ff_memory','--frame_frame_CA_cmap_memory',help='Flag for computing frame-frame C_alpha distance map using less memory (slower than -ca_ff) (optional)',action='store_true');
		parser.add_argument('-q','--query',help='Query for analyzing trajectory, e.g. -protein-, or -protein and noh-',default='protein');
		parser.add_argument('-si','--startID',help='Start residue ID if only checking between certain residues (only for single-chain atm)',default=-1);
		parser.add_argument('-ei','--endID',help='End residue ID if only checking between certain residues (only for single-chain atm)',default=-1);
		parser.add_argument('-bin','--binary_cmap',help='Make a binary contact map (optional)',action='store_true');
		parser.add_argument('-d_in','--distance_matrix_in',help='Input distance matrix - can be used if making binary cmap without computing a new cmap.',default='');
		parser.add_argument('-coff','--cutoff',help='Cutoff used to binarize the cmap.',default=0.5);
		parser.add_argument('-cmap_diff','--cmap_difference_to_start',help='Contact map distance to starting frame.', action='store_true');
		
		parser.add_argument('-cmap_diff_1_resid','--cmap_difference_to_start_one_resid',help='Contact map distance to starting frame for specific resid to all others.', action='store_true');
		parser.add_argument('-top','--topology_file',help='Input 1 topology file (.gro, .pdb, etc)',type=str,default='') #removed nargs=''
		parser.add_argument('-trj','--trajectory_files',help='Input trajectory files (.xtc, .dcd, etc)',nargs='+',default='')

		parser.add_argument('-build','--build_subunits',help='Superpose the sub-units (optional).',action='store_true')
		
		parser.add_argument('-multitraj','--multiple_trajectories',help='Flag for reading multiple trajectories. Need as many arguments in -top as in -trj',action='store_true')

		parser.add_argument('-fe','--file_end_name',type=str,help='Output file end name (optional)', default='')
		parser.add_argument('-od','--out_directory',type=str,help='The directory where data should be saved (optional)',default='')

		parser.add_argument('-nhtrj','--nat_holo_traj',help='Input native holo trajectory file (optional)',default='')
		parser.add_argument('-nhtop','--nat_holo_top',help='Input native holo topology file (optional)',default='')

		parser.add_argument('-natrj','--nat_apo_traj',help='Input native apo trajectory file (optional)',default='')
		parser.add_argument('-natop','--nat_apo_top',help='Input native apo topology file (optional)',default='')

		parser.add_argument('-dt','--dt',help='Keep every dt frame.',default=1)
		parser.add_argument('-downsample','--downsample_and_save',help='Downsample and save downsampled trajectories. The trajectories will be treated as continuum but saved as separate parts.',action='store_true')
		

      #load PDB file as traj
		args = parser.parse_args()
		
		startID = int(args.startID);		
		endID = int(args.endID);
		
		self.save_folder = args.out_directory;
		self.file_end_name = args.file_end_name;
		
		#for file in os.listdir(path_data):
			#if file.endswith('.pdb'):
				#self.traj = md.load_pdb(path_data+file) #args.topology_file
		if args.topology_file != '':
			self.traj = md.load_pdb(args.topology_file);
			print(self.traj)
			self.nResidues = int(self.traj.n_residues);
			
	 		if args.frame_frame_side_chain_cmap:
				self.computeFrameToFrameSideChainContacts(self.traj,args.query);
			
			if args.frame_frame_CA_cmap:
				self.computeFrameToFrameCalphaContacts(self.traj, args.query);
	
			if args.frame_frame_CA_cmap_memory:
				self.computeFrameToFrameCalpaContactsMemory(self.traj, args.query);
			
			if args.side_chain_cmap:
				self.computeAverageSideChainMinDistanceMap(startID, endID, args.query);
			
			if args.side_chain_cmap_semi_binary:
				self.computeAverageSideChainSemiBinCmap(startID, endID, args.query);
			
			if args.cmap_difference_to_start:
				self.computeCalpaCmapDistanceToFrame1(self.traj, args.query, args.cmap_difference_to_start_one_resid, startID);
			

		if args.binary_cmap:
			if args.distance_matrix_in != '':
				self.cmap = squareform(np.loadtxt(args.distance_matrix_in));
				print(self.cmap.shape)
			self.getContactMap(self.cmap, float(args.cutoff), makeBinary=True);
	

if __name__ == '__main__':
    
	parser = argparse.ArgumentParser(epilog='Residue-residue distance maps. Annie Westerlund 2017.');
	cmaps_obj = MD_cmaps();
	
	cmaps_obj.main(parser);
