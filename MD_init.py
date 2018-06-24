#!/usr/bin/env 

import mdtraj as md
import argparse
import numpy as np
import os

class MD_initializer:

	save_folder = ''
	file_name_end = ''
	sub_units = []
	
	def __init__(self):
		return
	
	def initialize_trajectory(self, parser):
		parser = self.setParserArguments(parser)
		args = parser.parse_args()
		
		# Get command line input parameters
		self.save_folder = args.out_directory
		self.file_name_end = args.file_end_name
		
		# Put / at end of out directory if not present. Check so that folder exists, otherwise construct it.
		if self.save_folder !='':
			if self.save_folder[-1] != '/':
				self.save_folder += '/'
				args.out_directory = self.save_folder
		
			if not os.path.exists(self.save_folder):
				os.makedirs(self.save_folder)

		print('Saving output files in directory: ' + self.save_folder)

		self.figure_counter = 1
		
		if not(args.build_subunits or args.downsample_and_save or args.multiple_trajectories):
			# Get the main trajectory
			traj = self.getTrajectory(args.topology_file[0], args.trajectory_files, float(args.dt))	
			self.setNativeTrajectories(args,traj)
		elif args.multiple_trajectories:
			traj = self.getMultipleTrajectories(args.topology_file, args.trajectory_files, float(args.dt))
		elif args.build_subunits:
			# Build sub-units if specified
			traj, self.sub_units = self.buildTrajectoryFromSubunits(args.topology_file)		

		elif args.downsample_and_save:
			self.downsampleTrajectories(args.topology_file, args.trajectory_files, float(args.dt))
			return [],[]
		
		print('File end name: ' + self.file_name_end)
		
		return traj, args

	def setParserArguments(self,parser):
		parser.add_argument('-top','--topology_file',help='Input 1 topology file (.gro, .pdb, etc)',type=str,default='',nargs='+')
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
		return parser

	def setNativeTrajectories(self,args,traj):
		# Set native trajectories
		if args.nat_holo_top != '':
			if args.nat_holo_traj != '':
				nat_traj = md.load(args.nat_holo_traj, top = args.nat_holo_top, stride=1)
				self.native_state_holo = nat_traj
				print('Holo native trajectory set')
			else:
				self.native_state_holo = traj
				print('Warning! Native holo state is the same as input trajectory. Use -nhtop and -nhtrj to set topology and trajectory files')
		else:
			self.native_state_holo = traj
			print('Warning! Native holo state is the same as input trajectory. Use -nhtop and -nhtrj to set topology and trajectory files.')



		if args.nat_apo_top != '':
			if args.nat_apo_traj != '':
				nat_traj = md.load(args.nat_apo_traj, top = args.nat_apo_top, stride=1)
				self.native_state_apo = nat_traj
				print('Apo native trajectory set')
			else:
				self.native_state_apo = traj
				print('Warning! Native apo state is the same as input trajectory. use -natop and -natrj to set topology and trajectory files.')
		else:
			self.native_state_apo = traj
			print('Warning! Native apo state is the same as input trajectory. Use -natop and -natrj to set topology and trajectory files.')
		return

	def getTrajectory(self, topology_file, trajectory_files, dt):
		trajectoryString = "Trajectory files: "
		topologyString = "Topology files: " + topology_file
	
		# Print file names in string
		for i in range(0,len(trajectory_files)):
			trajectoryString += trajectory_files[i] + " "
		
		print(topologyString)
		print(trajectoryString)
		
		
		# Stack trajectories
		traj = md.load(trajectory_files[0], top = topology_file, stride=1)
		
		
		for i in range(1,len(trajectory_files)):
			print("Stacking extra trajectory: " + str(i))
			
			trajectory = md.load(trajectory_files[i], top = topology_file, stride=1)
			# Keep every dt:th frame (dt=1 by default)
			traj = traj.join(trajectory)
			print("Number of frames: " + str(traj.n_frames))
		
		# Keep every dt:th frame (dt=1 by default)
		traj = traj[::int(dt)]
		return traj

	def downsampleTrajectories(self, topology_file, trajectory_files, dt):
		# Read trajectories and slice continuously. Can downsample trajectories that are too large to read all at once. 
		# The downsampled trajectories are saved in self.saveFolder.
		trajectoryString = "Trajectory files: "
		topologyString = "Topology files: " + topology_file
		
		# Print file names in string
		for i in range(0,len(trajectory_files)):
			trajectoryString += trajectory_files[i] + " "
		
		print(trajectoryString)
		print(topologyString)
		
		print("Slicing and saving trajectory: " + str(1))
		# Slice and save trajectories
		trajectory = md.load(trajectory_files[0], top = topology_file, stride=1)
		
		remainder = float(tmpTrajectory.n_frames)%dt
		
		# Keep every dt:th frame (dt=1 by default)
		trajectory = trajectory[::int(dt)]		
		trajectory.save_dcd(self.save_folder + self.file_name_end + str(0)+ '.dcd')

		if remainder == 0:
			start_frame = 0
		else:
			start_frame = dt - remainder
		
		print('Remainder: ' +str(remainder))
		print('start_frame: '+str(start_frame))
		
		for i in range(1,len(trajectory_files)):
			print("Slicing and saving trajectory: " + str(i+1))
			
			trajectory = md.load(trajectory_files[i], top = topology_file, stride=1)
			print(trajectory.n_frames)
			remainder = (float(trajectory.n_frames)-start_frame)%dt		
			trajectory = trajectory[start_frame::int(dt)]
			
			if remainder == 0:
				start_frame = 0
			else:
				start_frame = dt - remainder
			# Keep every dt:th frame (dt=1 by default). 
			# Start from start frame to treat trajectory as continuum.

			trajectory.save_dcd(self.save_folder + self.file_name_end + str(i) + '.dcd')
			
			print('Remainder: ' + str(remainder))
			print('start_frame: '+ str(start_frame))
		return		

	def getMultipleTrajectories(self, topology_files, trajectory_files, dt):
		
		trajs = []
		for i in range(len(topology_files)):
			trajs.append(self.getTrajectory(topology_files[i], [trajectory_files[i]], dt))
		
		return trajs

	def buildTrajectoryFromSubunits(self, topology_files):
		# Read a number of pdb files. Add them to the same molecule (different chains).
		print('Join sub-units')
		traj = md.load_pdb(topology_files[0])
		sub_units = [traj]
		for i in range(1,len(topology_files)):
			trajectory = md.load_pdb(topology_files[i], stride=1)
			traj = traj.stack(trajectory)
			sub_units.append(trajectory)		

		print('Write to file')
		traj.save_pdb(self.save_folder + self.file_name_end + ".pdb")		
		print(self.save_folder + self.file_name_end + ".pdb")
		return traj, sub_units
	
	def getSubUnits(self):
		return self.sub_units

	def getFileEndName(self):
		return self.file_name_end

	def getSaveFolder(self):
		return self.save_folder
	
	def getNativeStateHolo(self):
		return self.native_state_holo

	def saveTime(self):
		# Set parameters
		if traj.n_frames > 1:
			self.simulation_time = traj.time * traj.timestep / 1000.0 # Time [nanoseconds]
		else:
			self.simulation_time = np.zeros(1)
		
		# Save simulation time
		np.savetxt(self.save_folder + "time_" + self.file_name_end + ".txt" , self.simulation_time)
		return

	def filterTrajectories(self, traj, startID, endID):
		# Trims the trajectories given start and end resids.
		selection = 'residue '+ str(startID) + ' to ' + str(endID)
		
		indices = traj.topology.select(selection)
		traj = traj.atom_slice(indices)
		
		indices = self.native_state_holo.topology.select(selection)
		self.native_state_holo = self.native_state_holo.atom_slice(indices)
		
		indices = self.native_state_apo.topology.select(selection)
		self.native_state_apo = self.native_state_apo.atom_slice(indices)
		
		return traj

	def getStartAndEndResidues(self, traj):
		# Get the highest start residue and smallest end residue IDs for the trajectories.
		startID_vector = np.zeros(3)
		endID_vector = np.zeros(3)
		
		# Get start and end for input trajectory
		res = traj.topology.residue(0)
		tmpID = filter(str.isdigit, str(res))
		startID_vector[0] = int(tmpID) + 1
		
		res = traj.topology.residue(-1)
		tmpID = filter(str.isdigit, str(res))
		endID_vector[0] = int(tmpID) - 1
		
		# Get start and end for native holo trajectory				
		res = self.native_state_holo.topology.residue(0)
		tmpID = filter(str.isdigit, str(res))
		startID_vector[1] = int(tmpID) + 1 
		
		res = self.native_state_holo.topology.residue(-1)
		tmpID = filter(str.isdigit, str(res))
		endID_vector[1] = int(tmpID) - 1
		
		# Get start and end for native apo trajectory
		res = self.native_state_apo.topology.residue(0)
		tmpID = filter(str.isdigit, str(res))
		startID_vector[2] = int(tmpID) + 1
		
		res = self.native_state_apo.topology.residue(-1)
		tmpID = filter(str.isdigit, str(res))
		endID_vector[2] = int(tmpID) - 1
		
		# Compute the start and end ID
		startID = int(np.max(startID_vector))
		endID = int(np.min(endID_vector))
		
		return startID, endID
