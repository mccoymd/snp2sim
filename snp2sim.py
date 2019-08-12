#!/usr/bin/env python

import os
import sys
import argparse
import copy
import random
import multiprocessing
import yaml
import shutil
import subprocess
import csv
import time
import logging
import glob
import re

#Class that parses the command line and Yaml Config arguments
class argParse():
	#reads command line args, then config args
	#Config args supercede in case of conflict
	#Args are added to argParse parameters
	def __init__(self):
		self.requiredArgs = ['protein', 'mode']
		self.commandargs = _parseCommandLine()
		self.__dict__.update(self.commandargs.__dict__)
		if (self.config):
			self.args = yaml.safe_load(open(self.config))
			self.args = {k: v for k, v in self.args.items() if v is not None}
			self.__dict__.update(self.args)

		if self.verbose:
			logging.basicConfig(level=logging.DEBUG,
					format='%(asctime)s [%(name)-s] %(levelname)-s: %(message)s',
					datefmt='%Y-%m-%d %H:%M:%S')
		else:
			logging.basicConfig(level=logging.INFO,
					format='%(asctime)s [%(name)-s] %(levelname)-s: %(message)s',
					datefmt='%Y-%m-%d %H:%M:%S')
		self.logger = logging.getLogger('BASE')
	#Verifies all required general and module-specific args are provided
	#Throws assertionError if required parameter is not provided
	def checkRequiredArgs(self):
		for arg in self.requiredArgs:
			if not (hasattr(self, arg) and getattr(self, arg)):
				self.logger.error(arg + " not specified! " + arg + " is required.")
				sys.exit(1)
		if self.mode == "varMDsim":
			if not self.simID:
				self.simID = str(random.randint(1,1000))
				self.logger.debug("varMDsim ID: %s", self.simID)
			if not self.newStruct:
				self.logger.error("Template PDB not specified for simulation.")
				sys.exit(1)
			if not self.simLength:
				self.logger.error("Simulation length not specified.")
				sys.exit(1)
		if self.mode == "varScaffold":
			if not (hasattr(self, "scaffID") and getattr(self, "scaffID")):
				self.logger.error("scaffID not specified.")
				sys.exit(1)
			if not (hasattr(self, "alignmentResidues") and getattr(self, "alignmentResidues")):
				self.logger.error("Alignment Residues not specified.")
				sys.exit(1)
			if not (hasattr(self, "clusterResidues") and getattr(self, "clusterResidues")):
				self.logger.error("Cluster Residues not specified.")
				sys.exit(1)

		if self.mode == "varAnalysis":
			if not (hasattr(self, "analysisVariants") and self.analysisVariants):
				parameters.logger.error("Variants for analysis not specified.")
				sys.exit(1)

	#Sets the defaults for all other parameters, and initializes any other variables needed for the workflow.
	def setDefault(self):

		#programDir is the path for the main script and analysis scripts
		#as well as topology files and other supllemental files
		self.programDIR = os.path.abspath(__file__)
		self.programDIR = os.path.dirname(self.programDIR)

		#runDir is the path to save all results, and reference previous results or modules
		if not self.runDIR:
			self.runDIR = "/opt/snp2sim_results"
		if not os.path.isdir(self.runDIR):
				os.makedirs(self.runDIR)
		self.logger.debug("Run Directory: %s", self.runDIR)
		self.logger.debug("Installation Directory: %s", self.programDIR)
		if "." in self.protein:
			self.protein = self.protein.replace(".", "_")
		self.logger.info("Protein system: %s", self.protein)
		#set paths for tools
		if not self.VMDpath:
			self.VMDpath = "vmd"
		if not self.NAMDpath:
			self.NAMDpath = "namd2"
		if not self.PYTHONSHpath:
		   #make sure pythonsh (from AutoDockTools) has been added to path
			self.PYTHONSHpath = "/opt/mgltools_x86_64Linux2_1.5.6/bin/pythonsh"
		if not self.ADTpath:
			#must have path for "prepare_xxx.py" scripts from autodock tools
			self.ADTpath = "/opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
		if not self.VINApath:
			self.VINApath = "/opt/vina/autodock_vina_1_1_2_linux_x86/bin/vina"

		#Creates the variant parameter, a string that contains a unique ID for the variants used in the run
		#Covers the cases of the variants being provided as a list of AA and residue numbers
		#Defaults to WT if not variant provided
		if self.variant != "wt":
			if isinstance(self.varAA,list) and isinstance(self.varResID,list):
				if len(self.varAA) == len(self.varResID):
					self.variant = [str(self.varResID[x]) + self.varAA[x] for x in range(len(self.varAA))]
					self.logger.info("Variant: %s", ", ".join(self.variant))
					self.variant = "_".join(self.variant)

				else:
					self.logger.warning("Number of ResIDs does not match number of variant AAs.")
					self.logger.warning("Using WT structure for simulation")
					self.variant = "wt"
			elif self.varResID and self.varAA:
				self.variant = str(self.varResID) + self.varAA
				self.logger.info("Variant: %s", self.variant)
			else:
				self.logger.warning("varResID and varAA not specified")
				self.logger.warning("Using WT structure for simulation")
				self.variant = "wt"
		else:
			self.logger.info("Using WT structure for simulation")


		#Points to the supplemental files needed to create the starter files for the MD sim.
		self.simTopology = ("%s/simParameters/top_all36_prot.rtf" % self.programDIR,
								  "%s/simParameters/toppar_water_ions_namd.str" % self.programDIR)
		self.simParameters = ("%s/simParameters/par_all36_prot.prm" % self.programDIR,
									"%s/simParameters/toppar_water_ions_namd.str" %self.programDIR)

		#by default uses all cores
		if not self.simProc:
			self.simProc = multiprocessing.cpu_count()    
		self.logger.debug("Using %d cores", self.simProc)

		self.inputDIR = self.runDIR
		self.varsDIR = os.path.join(self.runDIR, "variantSimulations", self.protein) + "/"
		self.analysisDIR = os.path.join(self.runDIR, "variantSimulations", self.protein, "analysis/")
		self.runDIR = os.path.join(self.runDIR, "variantSimulations", self.protein, self.variant) + "/"
		self.structDIR = os.path.join(self.runDIR, "structures/")
		self.binDIR = os.path.join(self.runDIR, "bin/")
		self.configDIR = os.path.join(self.runDIR, "config/")
		self.resultsDIR = os.path.join(self.runDIR, "results/")
		self.drugBindingDIR = os.path.join(self.resultsDIR, "drugBinding/")
		self.trajDIR = os.path.join(self.resultsDIR, "trajectory/")
		self.scaffoldDIR = os.path.join(self.resultsDIR, "scaffold/")

		shutil.copy(self.config, self.configDIR)

		#creates paths to all the structures used in the prep for MD sim
		self.simPDB = self.structDIR + "%s.%s.pdb" \
							% (self.protein,self.variant)
		self.simPSF = self.structDIR + "%s.%s.psf" \
							% (self.protein,self.variant)
		self.templatePDB = self.structDIR + "%s.template.pdb" \
								 % (self.protein)
		self.wtPDB = self.structDIR + "%s.wt.UNSOLVATED.pdb" \
						   % (self.protein)
		self.wtPSF = self.structDIR + "%s.wt.UNSOLVATED.psf" \
						   % (self.protein)
							  
		self.varPrefix = self.structDIR + "%s.%s.UNSOLVATED" \
							   % (self.protein,self.variant)
		self.varPDB = self.structDIR + "%s.%s.UNSOLVATED.pdb" \
							% (self.protein,self.variant)
		self.varPSF = self.structDIR + "%s.%s.UNSOLVATED.psf" \
							% (self.protein,self.variant)

		#paths to the TCL scripts used to generate the starter structures
		self.wtStructTCL = self.binDIR + "%s.wt.genStructFiles.tcl" \
								 % (self.protein)
		self.varStructTCL = self.binDIR + "%s.%s.genStructFiles.tcl" \
								  % (self.protein,self.variant)
		self.singleRunTCL = self.binDIR + "%s.%s.genSingleRunPDB.tcl" \
								  % (self.protein,self.variant)
		self.solvTCL = self.binDIR + "%s.%s.genSolvStruct.tcl" \
							 % (self.protein,self.variant)
		self.solvBoundary = self.configDIR + "%s.%s.solvBoundary.txt" \
								  % (self.protein,self.variant)
		self.structPrefix = self.structDIR + "%s.%s" \
								  % (self.protein,self.variant)

		#path to NAMD parameters
		self.NAMDconfig = self.configDIR + "%s.%s.%s.NAMD" \
								% (self.protein, self.variant, self.simID)
		self.NAMDout = self.trajDIR + "%s.%s.%s" \
							 % (self.protein, self.variant, self.simID)

		#path to the drug libraries in the program directory
		self.drugLibPath = "%s/drugLibraries/%s/" % \
				   (self.programDIR, self.drugLibrary)

				
		#hardcoding bindingID as drugLibrary
		#TODO - refactor to remove "bindingID"
		self.bindingID = self.drugLibrary

		if not self.vinaExh or not isinstance(self.vinaExh, int):
			self.vinaExh = 50
		
		#Autodock config
		if self.bindingID:
			self.drugBindConfig = self.configDIR + "%s.autodock" \
										% (self.bindingID)
		

		#allows user to input trajectories in PDB format from external source
		if self.loadPDBtraj:
			variantDIR = self.trajDIR
			if not os.path.exists(variantDIR):
				os.makedirs(variantDIR)

			for pdbFile in self.loadPDBtraj:
				pdbNewLoc = self.trajDIR + os.path.basename(pdbFile)
				os.system("cp %s %s" % (pdbFile,pdbNewLoc))


def _parseCommandLine():

	parser = argparse.ArgumentParser(
		#version="0.1",
		description="snp2sim - Molecular Simulation of Somatic Variation",
		epilog="written by Matthew McCoy, mdm299@georgetown.edu"
	)


	#general parameters
	parser.add_argument("--mode",
						help="varMDsim, varScaffold, drugSearch, or varAnalysis",
						action="store",
						type=str,
						)
	parser.add_argument("--protein",
						help="name of protein system",
						action="store",
						type=str,
						)
	parser.add_argument("--variant",
						help="variant as \"wt\" or \"x###x\"",
						action="store",
						type=str,
						)
	parser.add_argument("--runDIR",
						help="path to directory to store results",
						action="store",
						type=str,
						)
	parser.add_argument("--clean",
						help="remove all files from previous run with same protein name",
						action="store_true",
						)
	parser.add_argument("--cgcRun",
						help="if Run is on CGC platform, move pdb and log to snp2sim root for processing",
						action="store_true",
						)

	#VarMDsim options

	parser.add_argument("--newStruct",
						help="path to cleaned PDB file (protein structure w/ cannonical aa)",
						action="store",
						type=str,
						)
	parser.add_argument("--bindingTemplate",
						help="path to PDB file used to create search space to align scaffold",
						action="store",
						type=str,
						)
	parser.add_argument("--simLength",
						help="varTraj simulation length in ns",
						action="store",
						type=float,
						)
	parser.add_argument("--varResID",
						help="variant residue ID from PDB template",
						action="store",
						type=str,
						)
	parser.add_argument("--varAA",
						help="amino acid to change",
						action="store",
						type=str,
						)
	parser.add_argument("--simID",
						help="ID number for the MD sim run",
						action="store",
						type=str,
						)
	parser.add_argument("--VMDpath",
						help="path to VMD executable",
						action="store",
						type=str,
						)
	parser.add_argument("--NAMDpath",
						help="path to NAMD executable",
						action="store",
						type=str,
						)
	parser.add_argument("--simProc",
						help="number of processors to run simulation",
						action="store",
						type=int,
						)
	parser.add_argument("--singleRun",
						help="output summary PDB trajectory only",
						action="store_true",
						)
	parser.add_argument("--genStructures",
						help="only generate initial structures",
						action="store_true",
						)
	

	#varScaff options
	parser.add_argument("--scaffID",
						help="ID for the Scaffolding run",
						action="store",
						type=str,
						)
	parser.add_argument("--clustPDBtraj",
						help="cluster results from multiple varScafold singleRun",
						action="store_true",
						)
	parser.add_argument("--loadPDBtraj", nargs='+',
						help="pdb trajectory files to import, list one after another",
						action="store",
						)
	parser.add_argument("--alignmentResidues",
						help="residue numbers used to align trajectory for clustering",
						action="store",
						type=str,
						)
	parser.add_argument("--clusterResidues",
						help="residue numbers to consider when clustering",
						action="store",
						type=str,
						)
	parser.add_argument("--featureTableMethod",
						help="use the feature table method, otherwise RMSD method is used",
						action="store_true",
						)
	parser.add_argument("--PairwiseRMSD",
						help="if RMSD method is used, pass in already created PairwiseRMSD matrix",
						action="store",
						type=str,
						)
	parser.add_argument("--colorTrajectory",
						help="create a full trajectory PDB file where each frame is colored by cluster",
						action="store_true",
						)

	#drugBinding options
	parser.add_argument("--vinaExh",
						help="exhaustiveness parameter of autodock vina",
						action="store",
						default="50",
						type=str,
						)
	parser.add_argument("--VINApath",
						help="path to AutoDock Vina executable",
						action="store",
						type=str,
						)
	parser.add_argument("--PYTHONSHpath",
						help="path to autodock tools python executable",
						action="store",
						type=str,
						)
	parser.add_argument("--ADTpath",
						help="path to autodock utilities directory",
						action="store",
						type=str,
						)
	parser.add_argument("--newBindingConfig",
						help="path to new binding config file",
						action="store",
						type=str,
						)
	parser.add_argument("--flexBinding",
						help="path to file with flex residue numbers",
						action="store",
						type=str,
						)
	parser.add_argument("--drugLibrary",
						help="name of snp2sim drug library",
						action="store",
						type=str,
						)
	parser.add_argument("--singleDrug",
						help="path to single drug PDBQT",
						action="store",
						type=str,
						)
	parser.add_argument("--inputScaff", nargs='+',
						help="pdb scaffold files to import, list one after another",
						action="store",
						)
	parser.add_argument("--bindSingleVar",
						help="only bind single variant scaffolds",
						action="store_true",
						)
	parser.add_argument("--ligandPDB",
						help="path to dir with ligand PDBs",
						action="store",
						type=str,
						)


	parser.add_argument("--config",
						help="use parameters from config.yaml. OVERWRITES command line parameters",
						action="store",
						)
	parser.add_argument("--verbose",
						help="show debugging logs",
						action="store_true",
						)
	#parses command line args
	cmdlineparameters, unknownparams = parser.parse_known_args()
	parameters = copy.copy(cmdlineparameters)

	return parameters




#creates unsolvated PDB and PSF files
#solvates the struct files
def genNAMDstructFiles(parameters):
	if not os.path.exists(parameters.binDIR):
		os.makedirs(parameters.binDIR)

	if os.path.isfile(parameters.templatePDB):
		parameters.logger.debug("Generating solvated/ionized PDB and PSF from %s", parameters.templatePDB)

		#creates the unsolvated PDB and PSF files if not already present
		#also mutates if necessary
		if os.path.isfile(parameters.varPDB) and os.path.isfile(parameters.varPSF):
			parameters.logger.debug("Unsolvated PDB and PSF exist")
			parameters.logger.debug("Using previously created unsolvated PDB and PSF")
		else:
			genStructTCL(parameters)
			genStructCommand = "%s -e %s" % (parameters.VMDpath, parameters.varStructTCL)
			parameters.logger.debug("Running unsolvated structure generation TCL in VMD")
			parameters.logger.debug("VMD Command: %s", genStructCommand)
			try:
				out = subprocess.check_output(genStructCommand, shell = True)
				parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
			except subprocess.CalledProcessError as e:
				parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
				sys.exit(1)

		#Solvates the structure files with a water box for simulation
		genSolvTCL(parameters)
		genSolvCommand = "%s -e %s" % (parameters.VMDpath, parameters.solvTCL)
		if not os.path.exists(parameters.configDIR):
			os.makedirs(parameters.configDIR)
		parameters.logger.debug("Running structure solvation/ionization TCL in VMD")
		parameters.logger.debug("VMD Command: %s", genSolvCommand)
		try:
			out = subprocess.check_output(genSolvCommand, shell = True)
			parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
		except subprocess.CalledProcessError as e:
			parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
			sys.exit(1)
	else:
		parameters.logger.error("No template exists")
		parameters.logger.error("Use --newStruct to specify clean PDB (only cannonical aa) ")
		sys.exit(1)


#Creates a TCL script for VMD
#uses Topology files with template PDB to make PDB and PSF files
#mutates with all variants
def genStructTCL(parameters):
	parameters.logger.debug("Creating unsolvated struct files")
	structFile = open(parameters.varStructTCL,"w+")
	parameters.logger.debug("Writing unsolvated structure generation TCL: %s", parameters.varStructTCL)
	structFile.write("package require psfgen\n")
	for tFile in parameters.simTopology:
		structFile.write("topology %s\n" % tFile)
	structFile.write("pdbalias residue HIS HSE\n")
	structFile.write("segment PROT {pdb %s}\n" % parameters.templatePDB)
	structFile.write("coordpdb %s PROT\n" % parameters.templatePDB)
	structFile.write("guesscoord\n")
	structFile.write("writepdb %s\n" % parameters.wtPDB)
	structFile.write("writepsf %s\n" % parameters.wtPSF)
	if not parameters.variant == "wt":
		longAA = { "G":"GLY","A":"ALA","L":"LEU","M":"MET","F":"PHE",
				   "W":"TRP","K":"LYS","Q":"GLN","E":"GLU","S":"SER",
				   "P":"PRO","V":"VAL","I":"ILE","C":"CYS","Y":"TYR",
				   "H":"HSE","R":"ARG","N":"ASN","D":"ASP","T":"THR"}
		structFile.write("package require mutator\n")
		if (isinstance(parameters.varAA, list) and isinstance(parameters.varResID, list)):
			for x in range(len(parameters.varAA)):
				structFile.write("mutator -psf %s -pdb %s -o %s -ressegname PROT -resid %s -mut %s\n" \
								% (parameters.wtPSF, parameters.wtPDB, \
									parameters.varPrefix, parameters.varResID[x], longAA.get(parameters.varAA[x])))
		else:
			structFile.write("mutator -psf %s -pdb %s -o %s -ressegname PROT -resid %s -mut %s\n" \
							% (parameters.wtPSF, parameters.wtPDB, \
								parameters.varPrefix, parameters.varResID, longAA.get(parameters.varAA)))
	
	structFile.write("quit\n")

#Solvates and ionizes the structure with water box
#measures center and dimensions of structure
def genSolvTCL(parameters):
	#todo - adjustable solvation/ionization parameters

	solvFile = open(parameters.solvTCL,"w+")
	parameters.logger.debug("Writing solvation/ionization TCL: %s", parameters.solvTCL)
	solvFile.write("package require solvate\n")
	solvFile.write("solvate %s %s -t 10 -o %s.wb\n" % \
				   (parameters.varPSF, parameters.varPDB, parameters.structPrefix))
	solvFile.write("package require autoionize\n")
	solvFile.write("autoionize -psf %s.wb.psf -pdb %s.wb.pdb -neutralize -o %s\n" % \
				   (parameters.structPrefix, parameters.structPrefix, parameters.structPrefix))
	solvFile.write("set box [atomselect top all]\n")
	solvFile.write("set output [open %s w]\n" % parameters.solvBoundary)
	solvFile.write("puts $output [measure center $box]\n")
	solvFile.write("puts $output [measure minmax $box]\n")
	solvFile.write("close $output\n")
	solvFile.write("quit")


#Runs NAMD for MD simulation using structures and options
def runNAMD(parameters):

	if not os.path.exists(parameters.configDIR):
		os.makedirs(parameters.configDIR)

	#parses the dimensions of the solvated structure
	if os.path.isfile(parameters.solvBoundary):
		boundaryFile = open(parameters.solvBoundary,"r")        
		bLines = boundaryFile.readlines()

		center = bLines[0]
		center = center.rstrip()
		center = center.split(" ")
		center = [float(i) for i in center]
		
		minmax = bLines[1]
		minmax = minmax.rstrip()
		minmax = minmax.replace("{","")
		minmax = minmax.replace("}","")
		minmax = minmax.split(" ")
		minmax = [float(i) for i in minmax]

		parameters.dimX = minmax[3] - minmax[0] + 0.1
		parameters.dimY = minmax[4] - minmax[1] + 0.1
		parameters.dimZ = minmax[5] - minmax[2] + 0.1

	else:
		parameters.logger.error("Boundary file does not exist")
		parameters.logger.error("Remove solvated/ionized PDB and PSF and resubmit")
		sys.exit(1)

	#creates the NAMD config file and populates it with default options	
	configFile = open(parameters.NAMDconfig,"w+")
	parameters.logger.debug("Writing NAMD config: %s", parameters.NAMDconfig)
	configFile.write("structure          %s\n" % parameters.simPSF)
	configFile.write("coordinates        %s\n" % parameters.simPDB)
	configFile.write("set outputname     %s\n" % parameters.NAMDout)
	configFile.write("paraTypeCharmm on\n")
	for pFile in parameters.simParameters:
		configFile.write("parameters %s\n" % pFile)
	configFile.write("mergeCrossterms yes\n")
	configFile.write("set temperature    310\n")
	configFile.write("firsttimestep      0\n")
	configFile.write("\n")
	configFile.write("temperature         $temperature\n")
	configFile.write("\n")
	configFile.write("\n")
	configFile.write("# Force-Field Parameters\n")
	configFile.write("exclude             scaled1-4\n")
	configFile.write("1-4scaling          1.0\n")
	configFile.write("cutoff              12.0\n")
	configFile.write("switching           on\n")
	configFile.write("switchdist          10.0\n")
	configFile.write("pairlistdist        14.0\n")
	configFile.write("\n")
	configFile.write("# Integrator Parameters\n")
	configFile.write("timestep            2.0  ;# 2fs/step\n")
	configFile.write("rigidBonds          all  ;# needed for 2fs steps\n")
	configFile.write("nonbondedFreq       1\n")
	configFile.write("fullElectFrequency  2\n")
	configFile.write("stepspercycle       10\n")
	configFile.write("\n")
	configFile.write("# Constant Temperature Control\n")
	configFile.write("langevin            on    ;# do langevin dynamics\n")
	configFile.write("langevinDamping     1     ;# damping coefficient (gamma) of 1/ps\n")
	configFile.write("langevinTemp        $temperature\n")
	configFile.write("langevinHydrogen    off    ;# don't couple langevin bath to hydrogens\n\n")
	configFile.write("#Periodic Boundary Conditions\n")
	configFile.write("cellBasisVector1 %.1f 0.0 0.0\n" % parameters.dimX)
	configFile.write("cellBasisVector2 0.0 %.1f 0.0\n" % parameters.dimY)
	configFile.write("cellBasisVector3 0.0 0.0 %.1f\n" % parameters.dimZ)
	configFile.write("cellOrigin %.1f %.1f %.1f\n" % \
					 (center[0], center[1], center[0]))
	configFile.write("margin 1.0\n\n")
	configFile.write("wrapAll on\n\n")
	configFile.write("# PME (for full-system periodic electrostatics)\n")
	configFile.write("PME                 yes\n")
	configFile.write("PMEGridSpacing      1.0\n")
	configFile.write("# Constant Pressure Control (variable volume)\n")
	configFile.write("useGroupPressure      yes ;# needed for rigidBonds\n")
	configFile.write("useFlexibleCell       no\n")
	configFile.write("useConstantArea       no\n")
	configFile.write("langevinPiston        on\n")
	configFile.write("langevinPistonTarget  1.01325 ;#  in bar -> 1 atm\n")
	configFile.write("langevinPistonPeriod  100.0\n")
	configFile.write("langevinPistonDecay   50.0\n")
	configFile.write("langevinPistonTemp    $temperature\n")
	configFile.write("# Output\n")
	configFile.write("outputName          $outputname\n")
	#todo - configure output options i.e. freq of output)
	configFile.write("restartfreq         10000     ;# 500steps = every 1ps\n")
	configFile.write("dcdfreq             5000\n")
	configFile.write("xstFreq             5000\n")
	configFile.write("outputEnergies      5000\n")
	configFile.write("outputPressure      5000\n")
	configFile.write("# Minimization\n")
	configFile.write("minimize            1000\n")
	configFile.write("reinitvels          $temperature\n\n")
	configFile.write("run                 %i\n" % (parameters.simLength*500000)) # using 2 fs step size
	if not os.path.isdir(parameters.trajDIR):
		os.makedirs(parameters.trajDIR)



#Runs NAMD for MD simulation using structures and options
def runNAMD_replicaExchange(parameters):
	if not os.path.exists(parameters.configDIR):
		os.makedirs(parameters.configDIR)

	#parses the dimensions of the solvated structure
	if os.path.isfile(parameters.solvBoundary):
		boundaryFile = open(parameters.solvBoundary,"r")        
		bLines = boundaryFile.readlines()

		center = bLines[0]
		center = center.rstrip()
		center = center.split(" ")
		center = [float(i) for i in center]
		
		minmax = bLines[1]
		minmax = minmax.rstrip()
		minmax = minmax.translate(None,"{}")
		minmax = minmax.split(" ")
		minmax = [float(i) for i in minmax]

		parameters.dimX = minmax[3] - minmax[0] + 0.1
		parameters.dimY = minmax[4] - minmax[1] + 0.1
		parameters.dimZ = minmax[5] - minmax[2] + 0.1

	else:
		parameters.logger.error("Boundary file does not exist")
		parameters.logger.error("Remove solvated/ionized PDB and PSF and resubmit")
		sys.exit(1)

	parameters.repConfig = parameters.configDIR + "replica_s%s.%s.%s.conf" \
								% (parameters.protein, parameters.variant, parameters.simID)
	repConfig = open(parameters.repConfig, "w+")
	parameters.logger.debug("Writing replica exchange config: %s", parameters.repConfig)
	repConfig.write("set num_replicas 8\n")
	repConfig.write("set min_temp 300\n")
	repConfig.write("set max_temp 600\n")
	repConfig.write("set steps_per_run 1000\n")
	repConfig.write("set num_runs %s\n" %(int(parameters.simLength*500)))
	repConfig.write("set runs_per_frame 5\n")
	repConfig.write("set frames_per_restart 2\n")
	repConfig.write("set namd_config_file \"%s\"\n" %parameters.NAMDconfig)
	repConfig.write("set output_root \"%s\"\n" %parameters.NAMDout)

	repConfig.write("set psf_file \"%s\"\n" %parameters.simPSF)
	repConfig.write("set initial_pdb_file \"%s\"\n" %parameters.simPDB)
	repConfig.write("set fit_pdb_file \"%s\"\n" %parameters.templatePDB)

	parameters.jobConfig = parameters.configDIR + "job_s%s.%s.%s.conf" \
								% (parameters.protein, parameters.variant, parameters.simID)
	jobConfig = open(parameters.jobConfig, "w+")

	jobConfig.write("source %s\n" %parameters.repConfig)
	if parameters.NAMDpath:
		jobConfig.write("if { ! [catch numPes] } { cd $(echo $(dirname %s)\"/lib/replica\") \nsource replica.namd }\n")
	else:
		jobConfig.write("if { ! [catch numPes] } { cd $(echo $(dirname $(which namd2))\"/lib/replica\") \nsource replica.namd }\n")




	#creates the NAMD config file and populates it with default options	
	configFile = open(parameters.NAMDconfig,"w+")
	parameters.logger.debug("Writing NAMD config: %s", parameters.NAMDconfig)
	configFile.write("structure          %s\n" % parameters.simPSF)
	configFile.write("coordinates        %s\n" % parameters.simPDB)
	#configFile.write("set outputname     %s\n" % parameters.NAMDout)
	configFile.write("paraTypeCharmm on\n")
	for pFile in parameters.simParameters:
		configFile.write("parameters %s\n" % pFile)
	configFile.write("mergeCrossterms yes\n")
	#configFile.write("set temperature    310\n")
	configFile.write("firsttimestep      0\n")
	configFile.write("\n")
	#configFile.write("temperature         $temperature\n")
	configFile.write("\n")
	configFile.write("\n")
	configFile.write("# Force-Field Parameters\n")
	configFile.write("exclude             scaled1-4\n")
	configFile.write("1-4scaling          1.0\n")
	configFile.write("cutoff              12.0\n")
	configFile.write("switching           on\n")
	configFile.write("switchdist          10.0\n")
	configFile.write("pairlistdist        14.0\n")
	configFile.write("\n")
	configFile.write("# Integrator Parameters\n")
	configFile.write("timestep            2.0  ;# 2fs/step\n")
	configFile.write("rigidBonds          all  ;# needed for 2fs steps\n")
	configFile.write("nonbondedFreq       1\n")
	configFile.write("fullElectFrequency  2\n")
	configFile.write("stepspercycle       10\n")
	configFile.write("\n")
	configFile.write("# Constant Temperature Control\n")
	#configFile.write("langevin            on    ;# do langevin dynamics\n")
	# configFile.write("langevinDamping     1     ;# damping coefficient (gamma) of 1/ps\n")
	#configFile.write("langevinTemp        $temperature\n")
	# configFile.write("langevinHydrogen    off    ;# don't couple langevin bath to hydrogens\n\n")
	configFile.write("#Periodic Boundary Conditions\n")
	configFile.write("cellBasisVector1 %.1f 0.0 0.0\n" % parameters.dimX)
	configFile.write("cellBasisVector2 0.0 %.1f 0.0\n" % parameters.dimY)
	configFile.write("cellBasisVector3 0.0 0.0 %.1f\n" % parameters.dimZ)
	configFile.write("cellOrigin %.1f %.1f %.1f\n" % \
					 (center[0], center[1], center[0]))
	configFile.write("margin 1.0\n\n")
	configFile.write("wrapAll on\n\n")
	configFile.write("# PME (for full-system periodic electrostatics)\n")
	configFile.write("PME                 yes\n")
	configFile.write("PMEGridSpacing      1.0\n")
	configFile.write("# Constant Pressure Control (variable volume)\n")
	configFile.write("useGroupPressure      yes ;# needed for rigidBonds\n")
	configFile.write("useFlexibleCell       no\n")
	configFile.write("useConstantArea       no\n")
	# configFile.write("langevinPiston        on\n")
	# configFile.write("langevinPistonTarget  1.01325 ;#  in bar -> 1 atm\n")
	# configFile.write("langevinPistonPeriod  100.0\n")
	# configFile.write("langevinPistonDecay   50.0\n")
	# configFile.write("langevinPistonTemp    $temperature\n")
	configFile.write("# Output\n")
	#configFile.write("outputName          $outputname\n")
	#todo - configure output options i.e. freq of output)
	configFile.write("restartfreq         10000     ;# 500steps = every 1ps\n")
	#configFile.write("dcdfreq             5000\n")
	configFile.write("xstFreq             5000\n")
	#configFile.write("outputEnergies      5000\n")
	configFile.write("outputPressure      5000\n")
	configFile.write("# Minimization\n")
	configFile.write("minimize            1000\n")
	#configFile.write("reinitvels          $temperature\n\n")
	configFile.write("run                 %i\n" % (parameters.simLength*500000)) # using 2 fs step size
					 
	if not os.path.isdir(parameters.trajDIR):
		os.makedirs(parameters.trajDIR)




				
#if running in CGC, prepares the results of MD sim as PDB file to export
def genSingleRunTCL(parameters):
	parameters.logger.debug("Writing PDB converter TCL: %s", parameters.singleRunTCL)
	pdbTCL = open(parameters.singleRunTCL,"w+")
	pdbTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	variantDIR = parameters.trajDIR

	#loads DCD files from the results into VMD
	parameters.logger.debug("using DCD files in %s", variantDIR)
	for tFile in os.listdir(variantDIR):
		if tFile.endswith(".dcd"):
			dcdFile = variantDIR + tFile
			pdbTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	pdbTCL.write("package require pbctools\n")
	pdbTCL.write("pbc wrap -center com -centersel \"protein\" -compound residue -all\n")

	#writes the trajectory as a PDB file instead of DCD for CGC runs
	allPDBoutput = "%s/%s.%s.%s.pdb" % (variantDIR, parameters.protein,
										parameters.variant, parameters.simID)
	pdbTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
				   % allPDBoutput)
	pdbTCL.write("quit")
	
	genPDBcommand = "%s -e %s" % (parameters.VMDpath, parameters.singleRunTCL)
	
	parameters.logger.debug("Running PDB converter TCL in VMD")
	parameters.logger.debug("VMD Command: %s", genPDBCommand)
	try:
		out = subprocess.check_output(genPDBCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)

#Creates a pairwise RMSD matrix of the trajectory frames for use in multidimensional scaling method
def calcPairwiseRMSD(parameters):

	#only calculated RMSD based on clusterResidues and aligns frames by alignmentResidues
	alignmentRes = parameters.alignmentResidues
	clusterRes = parameters.clusterResidues

	if not os.path.exists(parameters.scaffoldDIR):
		os.makedirs(parameters.scaffoldDIR)

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	if os.path.isfile(scaffLOG):
		parameters.logger.warning("Scaffold log already exists. Remove to recalculate Scaffolds")
		parameters.logger.warning("Creating scaffold PDBs using previous log")
		return
	#loads in the DCD files of the trajectory	
	parameters.logger.debug("Writing RMSD matrix TCL: %s", parameters.scaffoldTCL)
	clustTCL = open(parameters.scaffoldTCL,"w+")
	clustTCL.write("package require csv\n")
	clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	variantDIR = parameters.trajDIR
	parameters.logger.debug("Using DCD files in %s", variantDIR)
	for tFile in sorted(os.listdir(variantDIR)):
		if tFile.endswith(".dcd"):
			dcdFile = variantDIR + tFile
			clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	#aligns the frames to the first frame
	clustTCL.write("set nf [molinfo top get numframes]\n")
	clustTCL.write("set refRes [atomselect top \""+alignmentRes+"\" frame 0]\n")
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top \""+alignmentRes+"\" frame $i]\n")
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")

	clustTCL.write("set output [open %s w]\n" % parameters.featureTable)
	clustTCL.write("set back [atomselect top \""+clusterRes+"\"]\n")
	clustTCL.write("set refclustRes [atomselect top \""+clusterRes+"\" frame 0]\n")

	#Makes the matrix square by padding with 0s to avoid calculating RMSD between pairs twice
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("set row {}\n")
	clustTCL.write("for {set pad 0} {$pad <= $i} {incr pad} {\n")
	clustTCL.write("lappend row 0}\n")

	#calculates RMSD and writes to a file
	clustTCL.write("for {set j [expr $i + 1]} {$j < $nf} {incr j} {\n")
	clustTCL.write("$refclustRes frame $i \n")
	clustTCL.write("$back frame $j \n")
	clustTCL.write("set M [measure rmsd $back $refclustRes] \n")
	clustTCL.write("lappend row $M}\n")
	clustTCL.write("puts $output [::csv::join $row]}\n")


	clustTCL.write("close $output\n")
	clustTCL.write("quit")
	clustTCL.close()

	makeTableCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	parameters.logger.debug("Running RMSD matrix generation TCL in VMD")
	parameters.logger.debug("VMD Command: %s", makeTableCommand)
	try:
		out = subprocess.check_output(makeTableCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)

#Same function as calcPairwaiseRMSD(), but used a PDB format as input trajectory
def calcPairwiseRMSD_PDB(parameters):

	alignmentRes = parameters.alignmentResidues
	clusterRes = parameters.clusterResidues

	if not os.path.exists(parameters.scaffoldDIR):
		os.makedirs(parameters.scaffoldDIR)

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	if os.path.isfile(scaffLOG):
		parameters.logger.warning("Scaffold log already exists. Remove to recalculate Scaffolds")
		parameters.logger.warning("Creating scaffold PDBs using previous log")
		return
	parameters.logger.debug("Writing RMSD matrix TCL: %s", parameters.scaffoldTCL)
	clustTCL = open(parameters.scaffoldTCL,"w+")
	clustTCL.write("package require csv\n")
	variantDIR = parameters.trajDIR
	parameters.logger.debug("Using DCD files in %s", variantDIR)
	for trajFile in sorted(os.listdir(variantDIR)):
		if trajFile.endswith(".pdb"):
			pdbFile = variantDIR + trajFile
			clustTCL.write("mol addfile %s waitfor all\n" % pdbFile)
	clustTCL.write("set nf [molinfo top get numframes]\n")
	clustTCL.write("set refRes [atomselect top \""+alignmentRes+"\" frame 0]\n")
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top \""+alignmentRes+"\" frame $i]\n")
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")

	clustTCL.write("set output [open %s w]\n" % parameters.featureTable)
	clustTCL.write("set back [atomselect top \""+clusterRes+"\"]\n")
	clustTCL.write("set refclustRes [atomselect top \""+clusterRes+"\" frame 0]\n")

	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("set row {}\n")
	clustTCL.write("for {set pad 0} {$pad <= $i} {incr pad} {\n")
	clustTCL.write("lappend row 0}\n")

	clustTCL.write("for {set j [expr $i + 1]} {$j < $nf} {incr j} {\n")
	clustTCL.write("$refclustRes frame $i \n")
	clustTCL.write("$back frame $j \n")
	clustTCL.write("set M [measure rmsd $back $refclustRes] \n")
	clustTCL.write("lappend row $M}\n")
	clustTCL.write("puts $output [::csv::join $row]}\n")


	clustTCL.write("close $output\n")
	clustTCL.write("quit")
	clustTCL.close()

	makeTableCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	parameters.logger.debug("Running RMSD matrix generation TCL in VMD")
	parameters.logger.debug("VMD Command: %s", makeTableCommand)
	try:
		out = subprocess.check_output(makeTableCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)


# Changing the method of clustering - VS
def extractFeatureTable(parameters):
#TODO - currently config must be formated as:
#NOTE:: 3 lines, CASE SPECIFIC, NEEDS IMPROVEMENT
#alignmentRes "(atomselection)"
#clusterRes "(atomselection)"
#rmsdThresh (num)
	

	#old config method, integrated into yaml   - VS
	#no rmsd thresh necessary
	# scaffParameters = open(parameters.scaffParams,"r")
	# scaffLines = scaffParameters.readlines()

	# alignmentRes = scaffLines[0]
	# alignmentRes = alignmentRes.rstrip()
	# alignmentRes = alignmentRes.replace("alignmentRes ","")

	# clusterRes = scaffLines[1]
	# clusterRes = clusterRes.rstrip()
	# clusterRes = clusterRes.replace("clusterRes ","")

	# rmsdThresh = scaffLines[2]
	# rmsdThresh = rmsdThresh.rstrip()
	# rmsdThresh = rmsdThresh.replace("rmsdThresh ","")

	alignmentRes = parameters.alignmentResidues
	clusterRes = parameters.clusterResidues

	#regenerate config even if one exists to ensure new trajecories are included
	
	if not os.path.exists(parameters.scaffoldDIR):
		os.makedirs(parameters.scaffoldDIR)

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	if os.path.isfile(scaffLOG):
		parameters.logger.warning("Scaffold log already exists. Remove to recalculate Scaffolds")
		parameters.logger.warning("Creating scaffold PDBs using previous log")
		return
		
	parameters.logger.debug("Writing feature table TCL: %s", parameters.scaffoldTCL)
	clustTCL = open(parameters.scaffoldTCL,"w+")
	clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	variantDIR = parameters.trajDIR
	parameters.logger.debug("Using DCD files in %s", variantDIR)
	for tFile in sorted(os.listdir(variantDIR)):
		if tFile.endswith(".dcd"):
			dcdFile = variantDIR + tFile
			clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	#clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)
	clustTCL.write("set nf [molinfo top get numframes]\n")	

	clustTCL.write("set refRes [atomselect top \""+alignmentRes+"\" frame 0]\n")
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top \""+alignmentRes+"\" frame $i]\n")
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")

	clustTCL.write("set output [open %s w]\n" % parameters.featureTable)
	clustTCL.write("set back [atomselect top \""+clusterRes+"\"]\n")
	clustTCL.write("set refclustRes [atomselect top \""+clusterRes+"\" frame 0]\n")
	clustTCL.write("set backalpha [atomselect top \""+clusterRes+" and name CA"+"\"]\n")
	#header
	clustTCL.write("puts -nonewline $output \"frame\"\n")
	clustTCL.write("puts -nonewline $output \",rmsd\"\n")

	clustTCL.write("foreach resid [$backalpha get resid] resname [$backalpha get resname] {\n")
	clustTCL.write("puts -nonewline $output \",$resname-$resid-psi\"\n")
	clustTCL.write("puts -nonewline $output \",$resname-$resid-phi\"}\n")

	clustTCL.write("foreach atom [$back get name] resid [$back get resid] resname [$back get resname] {\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-x\"\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-y\"\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-z\"}\n")
	clustTCL.write("puts $output \"\"\n")

	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("$back frame $i \n")
	clustTCL.write("$backalpha frame $i \n")
	clustTCL.write("puts -nonewline $output \"$i\" \n")
	clustTCL.write("set M [measure rmsd $back $refclustRes] \n")
	clustTCL.write("puts -nonewline $output \",$M\" \n")
	clustTCL.write("set row [join [$backalpha get {psi phi}]]\n")
	clustTCL.write("foreach item $row {\n")
	clustTCL.write("puts -nonewline $output \",$item\"}\n")

	clustTCL.write("set row [join [$back get {x y z}]]\n")
	clustTCL.write("foreach item $row {\n")
	clustTCL.write("puts -nonewline $output \",$item\"}\n")
	clustTCL.write("puts $output \"\"}\n")

	clustTCL.write("close $output\n")

	# allPDBoutput = parameters.scaffBASE + ".all.pdb"
	# clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
	# 			   % allPDBoutput)

	clustTCL.write("quit")
	clustTCL.close()

	makeTableCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	parameters.logger.debug("Running feature table generation TCL in VMD")
	parameters.logger.debug("VMD Command: %s", makeTableCommand)
	try:
		out = subprocess.check_output(makeTableCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)


# Old method of clustering using VMD, changed to PCA + K means method - VS

	# clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	# variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
	# print("using DCD files in %s" % variantDIR)
	# for tFile in os.listdir(variantDIR):
	# 	if tFile.endswith(".dcd"):
	# 		dcdFile = variantDIR + tFile
	# 		clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	# clustTCL.write("package require pbctools\n")
	# clustTCL.write("pbc wrap -center com -centersel \"protein\" -compound residue -all\n")
	# clustTCL.write("set nf [molinfo top get numframes]\n")
	# clustTCL.write("set refRes [atomselect top %s frame 0]\n" % alignmentRes)
	# clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	# clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	# clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	# clustTCL.write("  set curRes [atomselect top %s frame $i]\n" % alignmentRes)
	# clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	# clustTCL.write("  $curStruct move $M\n")
	# clustTCL.write("}\n")
	# clustTCL.write("set output [open %s w]\n" % scaffLOG)
	
	# clustTCL.write("set clustRes [atomselect top %s]\n" % clusterRes)
	# clustTCL.write("puts $output \"clusters for RMSD Threshold %s\"\n" % rmsdThresh)
	# clustTCL.write("puts $output [measure cluster $clustRes distfunc rmsd cutoff %s]\n" \
	# 			   % rmsdThresh)
	# clustTCL.write("puts $output \"Only generated summary PDB file\"\n")
	# clustTCL.write("close $output\n")
	# allPDBoutput = parameters.scaffBASE + ".all.pdb"
	# clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
	# 			   % allPDBoutput)
	# clustTCL.write("quit")

#legacy code for surveying the trajectoies
#    for {set i 60} {$i < 80} {incr i} {
#                set thr [expr $i / 100.0]
#                puts $output [puts $thr]
#                puts $output [measure cluster $clustRes distfunc rmsd cutoff $thr]
#            }
	
#	vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
#        print("got here %s" % vmdClustCommand)
#	os.system(vmdClustCommand)

def extractFeatureTable_PDB(parameters):

	#old config method, integrated into yaml   - VS
	#no rmsd thresh necessary
	# scaffParameters = open(parameters.scaffParams,"r")
	# scaffLines = scaffParameters.readlines()

	# alignmentRes = scaffLines[0]
	# alignmentRes = alignmentRes.rstrip()
	# alignmentRes = alignmentRes.replace("alignmentRes ","")

	# clusterRes = scaffLines[1]
	# clusterRes = clusterRes.rstrip()
	# clusterRes = clusterRes.replace("clusterRes ","")

	# rmsdThresh = scaffLines[2]
	# rmsdThresh = rmsdThresh.rstrip()
	# rmsdThresh = rmsdThresh.replace("rmsdThresh ","")

	alignmentRes = parameters.alignmentResidues
	clusterRes = parameters.clusterResidues

	#regenerate config even if one exists to ensure new trajecories are included
	if not os.path.exists(parameters.scaffoldDIR):
		os.makedirs(parameters.scaffoldDIR)

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	if os.path.isfile(scaffLOG):
		parameters.logger.warning("Scaffold log already exists. Remove to recalculate Scaffolds")
		parameters.logger.warning("Creating scaffold PDBs using previous log")
		return

	parameters.logger.debug("Writing feature table TCL: %s", parameters.scaffoldTCL)

	if not os.path.exists(parameters.binDIR):
		os.makedirs(parameters.binDIR)

	clustTCL = open(parameters.scaffoldTCL,"w+")
	variantDIR = parameters.trajDIR
	parameters.logger.debug("Using PDB files in %s", variantDIR)
	for trajFile in sorted(os.listdir(variantDIR)):
		if trajFile.endswith(".pdb"):
			pdbFile = variantDIR + trajFile
			clustTCL.write("mol addfile %s waitfor all\n" % pdbFile)

	clustTCL.write("set nf [molinfo top get numframes]\n")	

	clustTCL.write("set refRes [atomselect top \""+alignmentRes+"\" frame 0]\n")
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top \""+alignmentRes+"\" frame $i]\n")
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")

	clustTCL.write("set output [open %s w]\n" % parameters.featureTable)
	clustTCL.write("set back [atomselect top \""+clusterRes+"\"]\n")
	clustTCL.write("set refclustRes [atomselect top \""+clusterRes+"\" frame 0]\n")
	clustTCL.write("set backalpha [atomselect top \""+clusterRes+" and name CA"+"\"]\n")
	#header
	clustTCL.write("puts -nonewline $output \"frame\"\n")
	clustTCL.write("puts -nonewline $output \",rmsd\"\n")

	clustTCL.write("foreach resid [$backalpha get resid] resname [$backalpha get resname] {\n")
	clustTCL.write("puts -nonewline $output \",$resname-$resid-psi\"\n")
	clustTCL.write("puts -nonewline $output \",$resname-$resid-phi\"}\n")

	clustTCL.write("foreach atom [$back get name] resid [$back get resid] resname [$back get resname] {\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-x\"\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-y\"\n")
	clustTCL.write("puts -nonewline $output \",$atom-$resname-$resid-z\"}\n")
	clustTCL.write("puts $output \"\"\n")

	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("$back frame $i \n")
	clustTCL.write("$backalpha frame $i \n")
	clustTCL.write("puts -nonewline $output \"$i\" \n")
	clustTCL.write("set M [measure rmsd $back $refclustRes] \n")
	clustTCL.write("puts -nonewline $output \",$M\" \n")
	clustTCL.write("set row [join [$backalpha get {psi phi}]]\n")
	clustTCL.write("foreach item $row {\n")
	clustTCL.write("puts -nonewline $output \",$item\"}\n")

	clustTCL.write("set row [join [$back get {x y z}]]\n")
	clustTCL.write("foreach item $row {\n")
	clustTCL.write("puts -nonewline $output \",$item\"}\n")
	clustTCL.write("puts $output \"\"}\n")

	clustTCL.write("close $output\n")

	# allPDBoutput = parameters.scaffBASE + ".all.pdb"
	# clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
	# 			   % allPDBoutput)

	clustTCL.write("quit")
	clustTCL.close()

	makeTableCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	parameters.logger.debug("Running feature table generation TCL in VMD")
	parameters.logger.debug("VMD Command: %s", makeTableCommand)
	try:
		out = subprocess.check_output(makeTableCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)

#Old method of clustering, switched to PCA + kmeans - VS
	# clustTCL = open(parameters.scaffoldTCL,"w+")
	# variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory"
	# if not os.path.exists(variantDIR):
	# 	os.makedirs(variantDIR)
		
	# for trajFile in os.listdir(variantDIR):
	# 	if trajFile.endswith(".pdb"):
	# 		pdbFile = variantDIR + "/" + trajFile
	# 		clustTCL.write("mol addfile %s waitfor all\n" % pdbFile)
	# clustTCL.write("package require csv\n")
	# clustTCL.write("set nf [molinfo top get numframes]\n")
	# clustTCL.write("set refRes [atomselect top %s frame 0]\n" % alignmentRes)
	# clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	# clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	# clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	# clustTCL.write("  set curRes [atomselect top %s frame $i]\n" % alignmentRes)
	# clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	# clustTCL.write("  $curStruct move $M\n")
	# clustTCL.write("}\n")
	# clustTCL.write("set output [open %s w]\n" % scaffLOG)
	# clustTCL.write("set clustRes [atomselect top %s]\n" % clusterRes)
	# #clustTCL.write("puts $output \"clusters for RMSD Threshold %s\"\n" % rmsdThresh)
	# clustTCL.write("set clust [measure cluster $clustRes distfunc rmsd cutoff %s]\n" \
	# 			   % rmsdThresh)
	# clustTCL.write("puts $output [::csv::joinlist $clust]\n")
	# clustTCL.write("close $output\n")


	# allPDBoutput = parameters.scaffBASE + ".all.pdb"
	# clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
	# 			   % allPDBoutput)
	# clustTCL.write("quit")

#legacy code for surveying the trajectoies
#    for {set i 60} {$i < 80} {incr i} {
#                set thr [expr $i / 100.0]
#                puts $output [puts $thr]
#                puts $output [measure cluster $clustRes distfunc rmsd cutoff $thr]
#            }

# moving this to runVarScaffold
#	vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
#	os.system(vmdClustCommand) 


def clusterTrajectory(parameters):
	scaffLOG = parameters.scaffBASE + ".log"
	clusterFigures = parameters.scaffoldDIR + "cluster_figures"
	if not os.path.isdir(clusterFigures):
		os.makedirs(clusterFigures)
	if parameters.featureTableMethod:
		clustCommand = "Rscript %s/snp2sim_analysis/scaffoldAnalysis/make_pca_plot.R %s %s %s" % \
											(parameters.programDIR, parameters.featureTable,
											clusterFigures, scaffLOG)
	else:
		clustCommand = "Rscript %s/snp2sim_analysis/scaffoldAnalysis/rmsd_mds_clustering.R %s %s %s" % \
											(parameters.programDIR, parameters.featureTable,
											clusterFigures, scaffLOG)

	parameters.logger.info("Running clustering script")
	parameters.logger.debug("Executing command: %s", clustCommand)
	try:
		r_out = subprocess.check_output(clustCommand, shell = True, stderr = subprocess.STDOUT)
		parameters.logger.debug("Clustering output: \n%s", r_out.decode("utf-8"))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("Clustering output: \n%s", e.output.decode('utf-8'))
		sys.exit(1)



def sortPDBclusters(parameters):

	def colorTrajectory():
		parameters.colorTraj = parameters.scaffoldDIR + "colorTraj.tcl"
		clustTCL = open(parameters.colorTraj, "w")
		allPDBoutput = parameters.scaffBASE + ".all.pdb"
		clustTCL.write("mol new %s waitfor all\n" %allPDBoutput)
		clustTCL.write("set cur [atomselect top \"%s\"]\n" %parameters.clusterResidues)
		clustTCL.write("set colorinput [open %s/cluster_by_frame.csv r]\n" % parameters.scaffoldDIR)
		clustTCL.write("gets $colorinput line\n")
		clustTCL.write("set clustMax 0\n")
		clustTCL.write("while {[gets $colorinput line]>=0} {\n\
						set fields [split $line \",\"]\n\
						animate goto [lindex $fields 0]\n\
						$cur frame [lindex $fields 0]\n\
						$cur set user [lindex $fields 1]\n\
						set clustMax [expr max($clustMax, [lindex $fields 1])]\n\
						}\n")
		

		clustTCL.write("mol modcolor 0 0 User\n")
		clustTCL.write("mol colupdate 0 0 1\n")
		clustTCL.write("mol scaleminmax 0 0 0.0 $clustMax\n")

	clustTCL = open(parameters.createClustPDB,"w+")
	parameters.logger.debug("Writing scaffold extraction TCL: %s", parameters.createClustPDB)
	scaffLOG = parameters.scaffBASE + ".log"
	variantDIR = parameters.trajDIR
	if parameters.clustPDBtraj:
		parameters.logger.debug("Loading PDB files in %s", variantDIR)
		for trajFile in sorted(os.listdir(variantDIR)):
			if trajFile.endswith(".pdb"):
				pdbFile = variantDIR + trajFile
				clustTCL.write("mol addfile %s waitfor all\n" % pdbFile)
	else:
		parameters.logger.debug("Loading DCD files in %s" % variantDIR)
		clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
		for trajFile in sorted(os.listdir(variantDIR)):
			if trajFile.endswith(".dcd"):
				dcdFile = variantDIR + trajFile
				clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	clustLogfile = open(scaffLOG, "r")
	clusterMembership = [x.split(",") for x in clustLogfile.readlines()]
	clusterMembership = [x for x in clusterMembership if x]
	scaffNum = 1
	for indCluster in clusterMembership:
		repStructIndex = int(indCluster[0])
		clustTCL.write("set cur [atomselect top \"protein not element H\" frame " + str(repStructIndex) + "]\n")
		clustTCL.write("$cur writepdb " + parameters.scaffBASE + "_cl" + str(scaffNum) + ".scaffold.pdb\n" )
		scaffNum += 1

	if parameters.colorTrajectory:
		allPDBoutput = parameters.scaffBASE + ".all.pdb"
		clustTCL.write("animate write pdb %s beg 0 end -1 sel $cur\n" \
					   % allPDBoutput)
		colorTrajectory()
	clustTCL.write("quit")
	clustTCL.close()

	genScaffolds = "%s -e %s" % (parameters.VMDpath, parameters.createClustPDB)
	parameters.logger.debug("Running Scaffold extraction TCL in VMD")
	parameters.logger.debug("VMD Command: %s", genScaffolds)
	try:
		out = subprocess.check_output(genScaffolds, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)


	#old method of extracting scaffolds to PDB
	#Doesnt require saving of concatenated PDB, saving space on drive
	# allPDBoutput = parameters.scaffBASE + ".all.pdb"
	# #traj = md.load(allPDBoutput)

	# allPDBstruct = open(allPDBoutput, "r")
	# allPDBlines = allPDBstruct.readlines()
	# pdbHeader = allPDBlines.pop(0)
	# indPDB = []
	# currPDB = ""
	# for pdbLine in allPDBlines:
	# 	currPDB = currPDB + pdbLine
	# 	if "END" in pdbLine:
	# 		indPDB.append(currPDB)
	# 		currPDB = ""
	
	# scaffLOG = parameters.scaffBASE + ".log"    
	# clustLogfile = open(scaffLOG, "r")
	# clusterMembership = [x.split(",") for x in clustLogfile.readlines()]
	# clusterMembership = [x for x in clusterMembership if x]
	
	# #Not necessary, no catch-all cluster in new method
	# #del clusterMembership[-1]
	# scaffNum = 1;
	# #print(traj.n_frames)
	# print(len(indPDB))
	# for indCluster in clusterMembership:
	# 	repStructIndex = int(indCluster[0])
	# 	#if len(indCluster.split(" ")) > parameters.clustThresh*traj.n_frames:
	# 	#    repStruct = traj[repStructIndex]
	# 	#    scaffFileName = parameters.scaffBASE + ".cl" + str(scaffNum) + ".scaffold.pdb"
	# 	#    repStruct.save(scaffFileName)
	# 	#    scaffNum += 1 
	# 	#print str(len(indCluster)) + " cluster " + str(scaffNum)
	# 	#if len(indCluster) > parameters.clustThresh*len(indPDB):
	# 		#scaffFileName = parameters.scaffBASE + ".cluster" + str(scaffNum) + ".pdb"

	# 	#no longer checks for clusterthresh in new method - VS
	# 	scaffFileName = parameters.scaffBASE + ".cl" + str(scaffNum) + ".scaffold.pdb"
	# 	scaffFile = open(scaffFileName,"w+")
	# 	scaffFile.write(pdbHeader)
	# 	scaffFile.write(indPDB[int(repStructIndex)])
	# 	scaffFile.close()
	# 	scaffNum += 1
	# #        for structID in indCluster:
	# #            print structID
	# #            if structID:
	# #                scaffFile.write(indPDB[int(structID)])
	# #    scaffNum += 1 

# def genScaffoldTCL(parameters):
	
# 	scaffParameters = open(parameters.scaffParams,"r")
# 	scaffLines = scaffParameters.readlines()

# 	alignmentRes = scaffLines[0]
# 	alignmentRes = alignmentRes.rstrip()
# 	alignmentRes = alignmentRes.replace("alignmentRes ","")

# 	scaffDIR = parameters.resultsDIR + "/" + parameters.variant + "/scaffold/"
# 	genScaff = open(parameters.clusterTCL, "w+")
# 	for tFile in os.listdir(scaffDIR):
# 		if tFile.endswith(".pdb"):
# 			if "cluster" in tFile:
# 				pdbClustFile = scaffDIR + tFile
# 				pdbScaffFile = pdbClustFile
# 				pdbScaffFile = pdbScaffFile.replace(".pdb",".scaffold.pdb")
# 				genScaff.write("mol new %s waitfor all\n" % pdbClustFile)
# #                genScaff.write("set domain [atomselect top %s]\n" % alignmentRes)
# 				genScaff.write("set domain [atomselect top all]\n")
# 				genScaff.write("set avePos [measure avpos $domain]\n")
# 				genScaff.write("$domain set {x y z} $avePos\n")
# 				genScaff.write("$domain writepdb %s\n" % pdbScaffFile)
				
# 	genScaff.write("quit\n")

#now defunct because no longer using mdtraj

# def genScaffoldMDTRAJ(parameters):
# 	#TODO generate rep structure using cluster res
# 	scaffParameters = open(parameters.scaffParams,"r")
# 	scaffLines = scaffParameters.readlines()

# 	clusterRes = scaffLines[1]
# 	clusterRes = clusterRes.rstrip()
# 	clusterRes = clusterRes.replace("clusterRes ","")

# 	scaffDIR = parameters.resultsDIR + "/" + parameters.variant + "/scaffold/"
# 	genScaff = open(parameters.clusterTCL, "w+")
# 	for tFile in os.listdir(scaffDIR):
# 		if tFile.endswith(".pdb"):
# 			if "cluster" in tFile:
# 				pdbClustFile = scaffDIR + tFile
# 				pdbScaffFile = pdbClustFile
# 				pdbScaffFile = pdbScaffFile.replace(".pdb",".scaffold.pdb")
# 				traj = md.load(pdbClustFile)
# 				atom_indices = [a.index for a in traj.topology.atoms if a.element.symbol != 'H']
# 				numStruct = traj.n_frames
# 				structLimit = 1500
# 				if numStruct > structLimit:
# 					numStruct = structLimit
# 					repSet = traj[random.sample(range(traj.n_frames),structLimit)]
# 				else:
# 					repSet = traj

# 				distances = np.empty((repSet.n_frames,repSet.n_frames))
# 				for i in range(repSet.n_frames):
# 					#print i
# 					distances[i] = md.rmsd(repSet,repSet,i,atom_indices = atom_indices)

# 				index = np.exp(-distances/distances.std()).sum(axis=1).argmax()
# 				print(index)
# 				#centroid = repSet[index]
# 				#centroid.save(pdbScaffFile)
		
def calcSearchSpace(parameters):
	out = parameters.drugBindConfig
	scaffRes = parameters.searchResidues
	parameters.boxTCL = parameters.binDIR + "%s_%s_%s_generateSearchSpace.tcl" % \
								   (parameters.protein, parameters.variant, parameters.bindingID)
	parameters.logger.debug("Writing TCL script to calculate search space bounds: %s", parameters.boxTCL)
	templatePDB = parameters.templatePDB
	clustTCL = open(parameters.boxTCL,"w+")
	clustTCL.write("set output [open %s w]\n" % out)

	clustTCL.write("mol new %s\n" % templatePDB)
	clustTCL.write("set pocket [atomselect top \"resid "+scaffRes+"\"]\n")
	clustTCL.write("set cen [measure center $pocket]\n")
	clustTCL.write("set bound [measure minmax $pocket]\n")

	clustTCL.write("puts -nonewline $output \"center_x = \"\n")
	clustTCL.write("puts $output [lindex $cen 0]\n")
	clustTCL.write("puts -nonewline $output \"center_y = \"\n")
	clustTCL.write("puts $output [lindex $cen 1]\n")
	clustTCL.write("puts -nonewline $output \"center_z = \"\n")
	clustTCL.write("puts $output [lindex $cen 2]\n")

	clustTCL.write("puts -nonewline $output \"size_x = \"\n")
	clustTCL.write("puts $output [expr 1.1 * [expr abs([lindex [lindex $bound 1] 0] - [lindex [lindex $bound 0] 0])]]\n")
	clustTCL.write("puts -nonewline $output \"size_y = \"\n")
	clustTCL.write("puts $output [expr 1.1 * [expr abs([lindex [lindex $bound 1] 1] - [lindex [lindex $bound 0] 1])]]\n")
	clustTCL.write("puts -nonewline $output \"size_z = \"\n")
	clustTCL.write("puts $output [expr 1.1 * [expr abs([lindex [lindex $bound 1] 2] - [lindex [lindex $bound 0] 2])]]\n")

	clustTCL.write("close $output\n")
	clustTCL.write("quit\n")
	clustTCL.close()

	vmdCommand = "%s -e %s" % (parameters.VMDpath, parameters.boxTCL)
	parameters.logger.debug("Running search space calculator TCL in VMD")
	parameters.logger.debug("VMD Command: %s", vmdCommand)
	try:
		out = subprocess.check_output(vmdCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)

def parseADconfig(parameters):
	paramFile = parameters.configDIR + "%s.autodock" %(parameters.bindingID)
	autodockParams = open(paramFile, "r")
	paramData = autodockParams.readlines()
	parameters.ADsearchSpace = [paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip()]

def getFlexRes(pdbFile,flexResIn):	
	currScaff = open(pdbFile,"r").readlines()
	currScaff.pop(0)
	flexRes = flexResIn[:]
	flexRes.sort()
	flexResNum = flexRes.pop(0)
	flexResID = ""
	for line in currScaff:
		resNum = line[22:26]
		if len(resNum) > 1:
			resNum = int(resNum)
		if flexResNum == resNum:
#            print resNum
			resName = line[17:20]
			flexResID = "%s_%s%s" % (flexResID,resName,str(flexResNum))
			if len(flexRes) > 0:

				flexResNum = flexRes.pop(0)
			else:
				flexResNum = 0
	return(flexResID[1:])

def alignScaff(parameters, currScaff):
	
	clusterRes = parameters.clusterResidues

	if not os.path.isdir(parameters.binDIR):
		os.makedirs(parameters.binDIR)

	
	alignmentConfig = parameters.binDIR + "%s.%s.alignStruct.TCL" \
							% (parameters.protein, parameters.variant)
	parameters.logger.debug("Writing scaffold alignment TCL: %s", alignmentConfig)

	alignmentTCL = open(alignmentConfig, "w+")
	alignmentTCL.write("mol new %s\n" % parameters.templatePDB)
	alignmentTCL.write("set ref [atomselect top \"%s\"]\n" % clusterRes)
	alignmentTCL.write("mol new %s\n" % currScaff)
	alignmentTCL.write("set target [atomselect top \"%s\"]\n" % clusterRes)
	alignmentTCL.write("set allP [atomselect top \"all\"]\n")
	alignmentTCL.write("set M [measure fit $target $ref]\n")
	alignmentTCL.write("$allP move $M\n")
	alignmentTCL.write("$allP writepdb %s\n" % currScaff)
	alignmentTCL.write("quit\n")
	alignmentTCL.close()
	alignmentCommand = "%s -e %s" % (parameters.VMDpath, alignmentConfig)

	parameters.logger.debug("Running scaffold alignment TCL in VMD")
	parameters.logger.debug("VMD Command: %s", alignmentCommand)
	try:
		out = subprocess.check_output(alignmentCommand, shell = True)
		parameters.logger.debug("VMD Output: \n\n%s\n\n", out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		parameters.logger.error("VMD output: \n%s", e.output.decode('ascii'))
		sys.exit(1)

def genVinaConfig(parameters):
	if hasattr(parameters, "numTrials") and parameters.numTrials > 1:
		for x in range(parameters.numTrials):
			vinaConfig = open(parameters.vinaConfigFolder+"/run%d.vina" %x,"w+")
			if not parameters.flexBinding:
				vinaConfig.write("receptor = %s\n" % parameters.scaff1out)
			else:
				vinaConfig.write("receptor = %s\n" % parameters.scaffRigid)
				vinaConfig.write("flex = %s\n" % parameters.scaffFlex)
			vinaConfig.write("ligand = %s\n" % parameters.currDrugPath)
			vinaConfig.write("out = %s/%s.%d.pdbqt\n" % (parameters.vinaOutDir, parameters.vinaBase, x))
			vinaConfig.write("log = %s/%s.%d.log\n" % (parameters.vinaOutDir, parameters.vinaBase, x))
			vinaConfig.write("exhaustiveness = %s\n" % parameters.vinaExh)
			for searchParam in parameters.ADsearchSpace:
				vinaConfig.write("%s\n" % searchParam)
			vinaConfig.close()
	else:
		vinaConfig = open(parameters.vinaConfig,"w+")
		if not parameters.flexBinding:
			vinaConfig.write("receptor = %s\n" % parameters.scaff1out)
		else:
			vinaConfig.write("receptor = %s\n" % parameters.scaffRigid)
			vinaConfig.write("flex = %s\n" % parameters.scaffFlex)
		vinaConfig.write("ligand = %s\n" % parameters.currDrugPath)
		vinaConfig.write("out = %s/%s.pdbqt\n" % (parameters.vinaOutDir, parameters.vinaBase))
		vinaConfig.write("log = %s/%s.log\n" % (parameters.vinaOutDir, parameters.vinaBase))
		vinaConfig.write("exhaustiveness = %s\n" % parameters.vinaExh)
		if hasattr(parameters, "vinaCPU") and isinstance(parameters.vinaCPU, int):
			vinaConfig.write("cpu = %d\n" % parameters.vinaCPU)
		for searchParam in parameters.ADsearchSpace:
			vinaConfig.write("%s\n" % searchParam)

#moving to runDrugSearch                
#	vinaCommand = "%s --config %s" % (parameters.VINApath, parameters.vinaConfig)
#	os.system(vinaCommand)





#### start of program move to main()
def runVarMDsim(parameters):
	if parameters.newStruct:
		if os.path.isdir(parameters.trajDIR) and not parameters.clean and not parameters.genStructures:
			parameters.logger.error("ERROR - %s-%s is already created... resubmit with new protein/variant combination", parameters.protein, parameters.variant)
			sys.exit(1)
		if hasattr(parameters, "clean") and parameters.clean:
			if os.path.isdir(parameters.trajDIR):
				shutil.rmtree(parameters.trajDIR)
			os.makedirs(parameters.trajDIR)
		if not os.path.isdir(parameters.structDIR):
			os.makedirs(parameters.structDIR)
		os.system("cp " + parameters.newStruct + " " + parameters.structDIR + "%s.template.pdb" %parameters.protein)
	if not os.path.isdir(parameters.trajDIR):
		os.makedirs(parameters.trajDIR)
	if os.path.isfile(parameters.simPDB):
		if os.path.isfile(parameters.simPSF):
			parameters.logger.debug("Using %s as simulation PDB", parameters.simPDB)
			parameters.logger.debug("Using %s as simulation PSF", parameters.simPSF)

			#if no variant structure exists, use clean template to create
		else:
			parameters.logger.debug("%s does not exist" % parameters.simPSF)
			parameters.logger.info("Generating pdb and psf for simulation")
			genNAMDstructFiles(parameters)
	else:
		parameters.logger.debug("%s does not exist" % parameters.simPDB)
		parameters.logger.info("Generating pdb and psf for simulation")
		genNAMDstructFiles(parameters)
	if parameters.genStructures:
		return
	if parameters.simLength:
		parameters.logger.info("Performing Variant %.3f ns Simulation", parameters.simLength)
		if hasattr(parameters, "replica") and parameters.replica:
			runNAMD_replicaExchange(parameters)
			os.makedirs("%s" %parameters.NAMDout)
			os.system("cd %s; mkdir 0 1 2 3 4 5 6 7" %parameters.NAMDout)
			runNAMDcommand = "mpirun --allow-run-as-root %s +p%i +replicas 8 %s > %s.log" % \
					 (parameters.NAMDpath, parameters.simProc,
					  parameters.jobConfig, parameters.NAMDout)
			parameters.logger.info("Running replica exchange NAMD")
			parameters.logger.debug("NAMD Command: %s", runNAMDcommand)
			try:
				out = subprocess.check_output(runNAMDcommand, shell = True)
				parameters.logger.debug("NAMD Output: \n\n%s\n\n", out.decode('ascii'))
			except subprocess.CalledProcessError as e:
				parameters.logger.error("NAMD output: \n%s", e.output.decode('ascii'))
				sys.exit(1)
			if not os.path.isfile("%s.dcd" % parameters.NAMDout):
				parameters.logger.error("NAMD run failed")
				sys.exit(1)
		else:
			runNAMD(parameters)
			runNAMDcommand = "%s +p%i %s > %s.log" % \
						 (parameters.NAMDpath, parameters.simProc,
						  parameters.NAMDconfig, parameters.NAMDout)
			parameters.logger.info("Running NAMD")
			parameters.logger.debug("NAMD Command: %s", runNAMDcommand)
			try:
				out = subprocess.check_output(runNAMDcommand, shell = True)
				parameters.logger.debug("NAMD Output: \n\n%s\n\n", out.decode('ascii'))
			except subprocess.CalledProcessError as e:
				parameters.logger.error("NAMD output: \n%s", e.output.decode('ascii'))
				sys.exit(1)
			if not os.path.isfile("%s.dcd" % parameters.NAMDout):
				parameters.logger.error("NAMD run failed")
				sys.exit(1)
		if parameters.singleRun:
			genSingleRunTCL(parameters)
			if parameters.cgcRun:
				CGCpdb = parameters.NAMDout + ".pdb"
				CGClog = parameters.NAMDout + ".log"
				cwd = os.getcwd()
				os.system("mv %s %s" % (CGCpdb, cwd))
				os.system("mv %s %s" % (CGClog, cwd))

	else:
		parameters.logger.error("Simulation length not specified")        
		sys.exit(1)

def runVarScaffold(parameters):
	#moved to start    
#    if parameters.newScaff:
#        if not os.path.isdir("%s/variantSimulations/%s/config" % \
#                             (parameters.runDIR,parameters.protein)):
#            os.makedirs("%s/variantSimulations/%s/config" % \
#                        (parameters.runDIR,parameters.protein))
#        if not os.path.isfile(parameters.scaffParams):
#            os.system("cp %s %s/variantSimulations/%s/config/%s.scaff" \
#                      % (parameters.newScaff,parameters.runDIR,
#                         parameters.protein, parameters.scaffID))
#
#        else:
#            print "Scaff config for %s exists." % (parameters.scaffParams)
#            print "Select new scaffID or remove existing scaff config"
#            sys.exit()

	if hasattr(parameters, "alignmentResidues") and hasattr(parameters, "clusterResidues"):
		#check if variant specified, otherwise perform clustering for all variants
		if parameters.variant and not isinstance(parameters.varAA, list) and not isinstance(parameters.varResID, list):
			variantList = (parameters.variant,)
		else:
			parameters.logger.info("Generating all variant scaffolds")
			variantList = ["".join([str(parameters.varResID[x]), str(parameters.varAA[x])]) for x in range(len(parameters.varAA))]

## TODO include option to analyze trajectory clusters to determine rmsd threshold
		for varSimResult in variantList:
			if hasattr(parameters, "clean") and parameters.clean:
				shutil.rmtree(parameters.scaffoldDIR)
				os.makedirs(parameters.scaffoldDIR)
			parameters.variant = varSimResult
			parameters.simPSF = parameters.structDIR + "%s.%s.psf" \
								% (parameters.protein,parameters.variant)
			parameters.scaffoldTCL =  parameters.binDIR + "%s.%s.%s.genFeatureTable.tcl" % \
									  (parameters.protein, parameters.variant, parameters.scaffID)
			parameters.scaffBASE = parameters.scaffoldDIR + "%s.%s.%s" % \
								   (parameters.protein, parameters.variant, parameters.scaffID)
			parameters.featureTable = parameters.configDIR + "%s.%s.%s.featureTable.csv" % \
								   (parameters.protein, parameters.variant, parameters.scaffID)
			parameters.createClustPDB = parameters.configDIR + "%s.%s.%s.createClustPDB.tcl" % \
								   (parameters.protein, parameters.variant, parameters.scaffID)

			parameters.logger.info("Generating scaffolds for %s", parameters.variant)
			#will not generate new TCL if previous scaffID log exists
			scaffLOG = parameters.scaffBASE + ".log"
			parameters.logger.debug("Scaffold log: %s", scaffLOG)
			if not os.path.isfile(scaffLOG):
				if parameters.featureTableMethod:
					if not parameters.clustPDBtraj:
						extractFeatureTable(parameters)
					else:
						extractFeatureTable_PDB(parameters)	

					clusterTrajectory(parameters)

					sortPDBclusters(parameters)
					#old (wrong) way to gen rep structure
					#genScaffoldTCL(parameters)
					#vmdScaffCommand = "%s -e %s" % (parameters.VMDpath, parameters.clusterTCL)
					#os.system(vmdScaffCommand)
					#genScaffoldMDTRAJ(parameters)

					if parameters.cgcRun:
						cwd = os.getcwd()
						scaffPDB = parameters.scaffBASE + "*.scaffold.pdb"
						os.system("mv %s %s" % (scaffLOG, cwd))
						os.system("mv %s %s" % (scaffPDB, cwd))

				else:
					if not hasattr(parameters, "PairwiseRMSD"):
						if not parameters.clustPDBtraj:
							calcPairwiseRMSD(parameters)
						else:
							parameters.scaffoldTCL =  parameters.binDIR + "%s.%s.%s.genFeatureTable.tcl" % \
													  (parameters.protein, parameters.variant, parameters.scaffID)
							calcPairwiseRMSD_PDB(parameters)	
					else:
						parameters.featureTable = parameters.PairwiseRMSD
					clusterTrajectory(parameters)

					sortPDBclusters(parameters)

					if parameters.cgcRun:
						cwd = os.getcwd()
						scaffPDB = parameters.scaffBASE + "*.scaffold.pdb"
						os.system("mv %s %s" % (scaffLOG, cwd))
						os.system("mv %s %s" % (scaffPDB, cwd))
			else:
					parameters.logger.warning("Scaffold log already exists. Remove to recalculate Scaffolds")
					parameters.logger.warning("Creating scaffold PDBs using previous log")
					sortPDBclusters(parameters)

					if parameters.cgcRun:
						cwd = os.getcwd()
						scaffPDB = parameters.scaffBASE + "*.scaffold.pdb"
						os.system("mv %s %s" % (scaffLOG, cwd))
						os.system("mv %s %s" % (scaffPDB, cwd))
	else:
		parameters.logger.error("Scaffold Parameters (alignmentResidues and clusterResidues) do not exist")
		sys.exit(1)


def runDrugSearch(parameters):

	#copy over the binding config with the search coords
	#drugBindConfig is the location of the bidning config
	parameters.vinaOutDir = parameters.drugBindingDIR
	if not os.path.isdir(parameters.vinaOutDir):
		os.makedirs(parameters.vinaOutDir)

	#Either determine search space automatically or take file input
	if hasattr(parameters, "autoSearchSpace") and parameters.autoSearchSpace and \
		hasattr(parameters, "searchResidues") and parameters.searchResidues:
		calcSearchSpace(parameters)
	else:
		if parameters.newBindingConfig:
			configDIR = parameters.configDIR
			if not os.path.isdir(configDIR):
				os.makedirs(configDIR)
				
			shutil.copyfile(parameters.newBindingConfig, os.path.join(configDIR, "%s.autodock" %parameters.bindingID))

	if os.path.isfile(parameters.drugBindConfig):
		parameters.logger.debug("Using search space config %s", parameters.newBindingConfig)
		parseADconfig(parameters)
	else:
		parameters.logger.error("Search space config file does not exist at path: %s", parameters.drugBindConfig)
		sys.exit(1)

	#copying over the flex config with residue numbers	

	if parameters.flexBinding:
		parameters.flexBinding = parameters.flexBinding.rstrip().split(" ")
		parameters.flexBinding = [int(x) for x in parameters.flexBinding]
		parameters.logger.debug("Flex residues: %s", ", ".join(str(x) for x in parameters.flexBinding))

	#copies input scaffolds into a scaffold folder to be read in later
	if parameters.inputScaff:
		parameters.logger.debug("Using input scaffolds")
		variantDIR = parameters.scaffoldDIR
		#if not os.path.exists(variantDIR):
		#    os.makedirs(variantDIR)
		os.makedirs(variantDIR)
		inputCount = 1
		for pdbFile in os.listdir(parameters.inputScaff):
#            pdbNewLoc = variantDIR + os.path.basename(pdbFile)
#            pdbNewLoc = os.path.splitext(pdbNewLoc)[0] + ".scaffold.pdb"
			inputID = "cl%s" % str(inputCount)
			pdbNewLoc = "%s/%s.%s.input_%s.scaffold.pdb" % (variantDIR, parameters.protein,
													  parameters.variant, inputID)
			inputCount += 1
			shutil.copyfile(pdbFile, pdbNewLoc)

	#if pdb files given as ligands, converts to pdbqt and places in drugLibrary directory
	if hasattr(parameters, "ligandPDB") and parameters.ligandPDB:
		library = parameters.programDIR + "/drugLibraries/" + os.path.basename(parameters.ligandPDB)
		if not os.path.isdir(library):
			os.makedirs(library)
		for file in os.listdir(parameters.ligandPDB):
			prepLigand = "%s %s/prepare_ligand4.py -l %s" \
										 % (parameters.PYTHONSHpath, parameters.ADTpath,
											file)
			os.system(prepLigand)
			os.system("cp %sqt %s/%sqt" % (file, library, file))
		parameters.drugLibPath = library

	#sets the drugLibrary, and creates a new library if a single drug is given
	if parameters.drugLibrary:
		if parameters.singleDrug:
			#todo - add single drug binding
			parameters.logger.debug("Binding a single drug: %s", parameters.singleDrug)
			if not os.path.isdir(parameters.drugLibPath):
				parameters.logger.debug("Drug library does not exist! Creating it!")
				os.makedirs(parameters.drugLibPath)
				shutil.copyfile(parameters.singleDrug, os.path.join(parameters.drugLibPath, os.path.basename(parameters.singleDrug)))
		else:
			if not os.path.isdir(parameters.drugLibPath):
				parameters.logger.error("Drug library does not exist")
				sys.exit(1)
		parameters.logger.debug("Using small molecules in %s" % parameters.drugLibPath)
			
										
	if parameters.bindingID:
		if not os.path.isfile(parameters.templatePDB):
			if parameters.bindingTemplate:
				if not os.path.isdir(parameters.structDIR):
					os.makedirs(parameters.structDIR)
				shutil.copyfile(parameters.bindingTemplate, parameters.structDIR + "%s.template.pdb" \
						  % (parameters.protein))
			else:
				parameters.logger.error("No binding template specified or template PDB available")
				sys.exit(1)
	   
				
		if parameters.bindSingleVar:
			bindingVar = [parameters.variant,]
		else:
			bindingVar = os.listdir(parameters.varsDIR)
			if "analysis" in bindingVar:
				bindingVar.remove("analysis")
			
		for var in bindingVar:
			#regenerate pdbqt for all scaffold.pdb files
			parameters.logger.info("Binding to variant %s", var)
			#todo check if pdbqt exist already
			scaffDIR = parameters.varsDIR + var
			for scaffPDB in os.listdir(scaffDIR):
				if scaffPDB.endswith("scaffold.pdb"):
					parameters.logger.debug("Binding to scaffold: %s", scaffPDB)
					currScaffPath = scaffDIR + scaffPDB
					#align scaff to template
					alignScaff(parameters,currScaffPath)
					scaffBase = scaffDIR + os.path.splitext(scaffPDB)[0]
					scaffBase = os.path.splitext(scaffBase)[0]

					parameters.scaff1out = scaffBase + ".pdbqt"
					#prepBaseScaff =  "%s %s/prepare_receptor4.py -U nphs -r %s -o %s" \
					prepBaseScaff =  "%s %s/prepare_receptor4.py -r %s -o %s" \
														 % (parameters.PYTHONSHpath,
															parameters.ADTpath,
															currScaffPath,
															parameters.scaff1out)
					parameters.logger.debug("Preparing scaffold as receptor")
					parameters.logger.debug("Prepare receptor command: %s", prepBaseScaff)
					try:
						out = subprocess.check_output(prepBaseScaff, shell = True)
						parameters.logger.debug("Prepare Receptor Script Output: %s", out.decode('ascii'))
					except subprocess.CalledProcessError as e:
						parameters.logger.error("Prepare Receptor Script Output: \n%s", e.output.decode('ascii'))
						sys.exit(1)

					if parameters.flexBinding:
						parameters.logger.debug("Found Flex Config: Performing flexable residue binding")
						
						scaffFlexRes = getFlexRes(currScaffPath,parameters.flexBinding)
						parameters.logger.debug("Residues in Flex Binding region: %s", scaffFlexRes)
						parameters.scaffFlex = scaffBase + ".flex.pdbqt"
						parameters.scaffRigid = scaffBase + ".rigid.pdbqt"
						prepFlexScaff = "%s %s/prepare_flexreceptor4.py -r %s -s %s -g %s -x %s" \
									 % (parameters.PYTHONSHpath, parameters.ADTpath,
										parameters.scaff1out, scaffFlexRes,
										parameters.scaffRigid, parameters.scaffFlex)
						parameters.logger.debug("Preparing scaffold as flex receptor")
						parameters.logger.debug("Prepare flex receptor command: %s", prepFlexScaff)
						try:
							out = subprocess.check_output(prepFlexScaff, shell = True)
						except subprocess.CalledProcessError as e:
							parameters.logger.error("Prepare Receptor Script Output: \n%s", e.output.decode('ascii'))
							sys.exit(1)

					for dFile in os.listdir(parameters.drugLibPath):
						if dFile.endswith(".pdbqt"):
							parameters.currDrugPath = parameters.drugLibPath + dFile
							parameters.drugBase = os.path.splitext(dFile)[0]
							parameters.vinaBase = os.path.basename("%s.%s.%s" % \
																   (scaffBase,parameters.bindingID,
																	parameters.drugBase))

							if hasattr(parameters, "numTrials") and parameters.numTrials > 1:
								parameters.vinaConfigFolder = parameters.configDIR + parameters.vinaBase
								if not os.path.exists(parameters.vinaConfigFolder):
									os.makedirs(parameters.vinaConfigFolder)
								genVinaConfig(parameters)
								for config in os.listdir(parameters.vinaConfigFolder):
									vinaCommand = "%s --config %s" % (parameters.VINApath, parameters.vinaConfigFolder + "/" + config)
									parameters.logger.debug("Running Vina")
									parameters.logger.debug("Vina command: %s", vinaCommand)
									os.system(vinaCommand)
							else:		
								parameters.vinaConfig = parameters.configDIR + "%s.vina"

								genVinaConfig(parameters)
								vinaCommand = "%s --config %s" % (parameters.VINApath, parameters.vinaConfig)
								parameters.logger.debug("Running Vina")
								parameters.logger.debug("Vina command: %s", vinaCommand)
								os.system(vinaCommand)

						
					if parameters.cgcRun:
						cwd = os.getcwd()
						vinaOUT = parameters.vinaOutDir + "*.pdbqt"
						vinaLOG = parameters.vinaOutDir + "*.log"
						os.system("cp %s %s" % (vinaLOG, cwd))
						os.system("cp %s %s" % (vinaOUT, cwd))

			
	# - - align to reference pdb
	# - - create pdbqt from file
	# - - create flex/rigid pdbqt from pdbqt
	# - Check for proper flex/rigid files for given
	# - - if single run, set variable as file path
	# - - else set variable as multiple file paths

					
	else:
		parameters.logger.error("Specify bindingID")
		sys.exit(1)
		

	#import new scaffolds - need config defining flex residues
	# - if new pdbqt


	#import new drugs
	# - if new pdbqt
	# - - create pdbqt from file
	# - Check for drug pdbqt from parameter name
	# - - if single run, set variable as file path
	# - - else set variable as multiple file paths
	
	#for each scaffold
	# - for each drug
	# - - bind drug to single scaffold
	# - - - build config - need search space 

#Runs analysis scripts - VS
def runAnalysis(parameters):

	# if parameters.variant == "wt" and not hasattr(parameters, "varAA") and not hasattr(parameters, "varResID"):
	# 	print("generating figures for all Variant Scaffolds")
	# 	variantList = os.listdir(parameters.resultsDIR)
	# elif (isinstance(parameters.varAA, list) and isinstance(parameters.varResID, list)):
	# 	print("generating figures for variants: %s" % (parameters.variant))
	# 	variantList = [str(self.varResID[x]) + self.varAA[x] for x in range(len(self.varAA))]
	# 	print("generating figures for variants: %s" % (", ".join(variantList)))
	# else:
	# 	print("generating figures for %s ONLY" %parameters.variant)
	# 	variantList = (parameters.variant, )


	variantList = parameters.analysisVariants
	parameters.logger.info("Running analysis on variants: %s", ", ".join(variantList))

	for v in variantList:
		if v not in os.listdir(parameters.varsDIR):
			parameters.logger.error("Variant " + v + " not in results directory!")
			sys.exit(1)

	

	analysisVars = "_".join(variantList)
	if hasattr(parameters, "analysisDrugLibrary"):
		analysisVars = parameters.analysisDrugLibrary + "_" + analysisVars
	elif hasattr(parameters, "analysisDrug"):
		analysisVars = parameters.analysisDrug + "_" + analysisVars
	curAnalysisDir = "%s/%s/binding_structures" % \
											(parameters.analysisDIR, analysisVars)
	if not os.path.isdir(curAnalysisDir):
		os.makedirs(curAnalysisDir)

	if hasattr(parameters, "analysisDrugLibrary"):
		if not isinstance(parameters.analysisDrugLibrary, list):
			parameters.analysisDrugLibrary = [parameters.analysisDrugLibrary]
		for v in variantList:
			for scaff in os.listdir(os.path.join(parameters.varsDIR, v, "results", "drugBinding")):
				if scaff.endswith(".pdbqt"):
					if scaff.split(".")[3] in parameters.analysisDrugLibrary:
						path = os.path.join(parameters.varsDIR, v, "results", "drugBinding", scaff)
						shutil.copy(path, curAnalysisDir)
	elif hasattr(parameters, "analysisDrug"):
		if not isinstance(parameters.analysisDrug, list):
			parameters.analysisDrug = [parameters.analysisDrug]
		for v in variantList:
			for scaff in os.listdir(os.path.join(parameters.varsDIR, v, "results", "drugBinding")):
				if scaff.endswith(".pdbqt"):
					if scaff.split(".")[4] in parameters.analysisDrug:
						path = os.path.join(parameters.varsDIR, v, "results", "drugBinding", scaff)
						shutil.copy(path, curAnalysisDir)
	else:	
		for v in variantList:
			for scaff in os.listdir(os.path.join(parameters.varsDIR, v, "results", "drugBinding")):
				if scaff.endswith(".pdbqt"):
					path = os.path.join(parameters.varsDIR, v, "results", "drugBinding", scaff)
					shutil.copy(path, curAnalysisDir)

	files = [len(file.split(".")) for file in os.listdir(curAnalysisDir)]
	if not len(set(files)) <= 1:
		parameters.logger.error("Some files are from multiple trial runs and others from single trial runs.")
		parameters.logger.error("Run varAnalysis on only a single type of run results")
		sys.exit(1)
	if files[0] == 7:
		parameters.mulTrials = True
		parameters.logger.debug("Creating summary file")
		os.system("%s/snp2sim_analysis/vinaAnalysis/collectResults_mulTrials.sh %s %s" \
				%(parameters.programDIR, curAnalysisDir, parameters.analysisDIR + analysisVars))
	else:
		parameters.mulTrials = False
		parameters.logger.debug("Creating summary file")
		os.system("%s/snp2sim_analysis/vinaAnalysis/collectResults.sh %s %s" \
				%(parameters.programDIR, curAnalysisDir, parameters.analysisDIR + analysisVars))


	#r_out = subprocess.check_output("Rscript %s/snp2sim_analysis/vinaAnalysis/visualizations.R %s/analysis/%s/vinaSummary_%s.txt" % (parameters.programDIR, 
	#											parameters.resultsDIR, analysisVars, parameters.protein), shell = True)
	#print(r_out)

	if not hasattr(parameters, "customScaff") or not parameters.customScaff:
		summary = "%s/vinaSummary_%s.txt" % (parameters.analysisDIR + analysisVars, parameters.protein)
		summary = open(summary, "r")
		summary = [x.split() for x in summary.readlines()]


		for v in variantList:
			if len(glob.glob("%s/%s/results/old_scaffold/%s.%s.*.log" % (parameters.varsDIR, v, parameters.protein, v))) == 1:
				scaffLOG = glob.glob("%s/%s/results/old_scaffold/%s.%s.*.log" % (parameters.varsDIR, v, parameters.protein, v))[0]
				shutil.copy(scaffLOG, parameters.analysisDIR + analysisVars)

				#if hasattr(parameters, "legacyScaff") and parameters.legacyScaff:
				#	os.system("%s/snp2sim_analysis/scaffoldAnalysis/legacy_scaffold.sh %s" %(parameters.programDIR, scaffLOG))
				clustLogfile = open(scaffLOG, "r")
				clusterMembership = [x.split(",") for x in clustLogfile.read().splitlines()]
				total_frames = sum(len(x) for x in clusterMembership)
				if hasattr(parameters, "legacyScaff") and parameters.legacyScaff:
					del clusterMembership[-1]
					clusterMembership = [x for x in clusterMembership if x and len(x) > .09 * total_frames]
				total_frames = sum(len(x) for x in clusterMembership)
				propMem = [len(x)*1.0/total_frames for x in clusterMembership]

				for row in summary[1:]:
					if row[4] == v:
						clust = int(re.search('(\d+)$', row[5]).group(0)) - 1
						row.append(propMem[clust])
			elif len(glob.glob("%s/%s/results/scaffold/%s.%s.*.log" % (parameters.varsDIR, v, parameters.protein, v))) > 1:
				parameters.logger.error("More than one scaffold log in scaffold directory")
				sys.exit(1)
			else:
				parameters.logger.error("Scaffold log does not exist in %s/%s/results/scaffold/ does not exist" % (parameters.varsDIR, v))
				sys.exit(1)
		summary[0].append("weight")

		with open("%s/weighted_vinaSummary_%s.txt" % (parameters.analysisDIR + analysisVars, parameters.protein),"w") as tsv:
			csvWriter = csv.writer(tsv, delimiter='\t')
			csvWriter.writerows(summary)

		curAnalysisDir = "%s/figures" % \
												(parameters.analysisDIR + analysisVars)
		if not os.path.isdir(curAnalysisDir):
			os.makedirs(curAnalysisDir)

		analysisCommand = "Rscript %s/snp2sim_analysis/vinaAnalysis/visualizations_weighted.R %s/weighted_vinaSummary_%s.txt %s" % (parameters.programDIR, 
													parameters.analysisDIR + analysisVars, parameters.protein, parameters.mulTrials)
		parameters.logger.info("Running analysis script")
		parameters.logger.debug("Executing command: %s", analysisCommand)
		try:
			r_out = subprocess.check_output(analysisCommand, shell = True, stderr = subprocess.STDOUT)
			parameters.logger.debug("Analysis output: \n%s", r_out.decode("utf-8"))
		except subprocess.CalledProcessError as e:
			parameters.logger.error("Clustering output: \n%s", e.output.decode('utf-8'))
			sys.exit(1)

		with open("%s/interactive_visualize.sh" %(parameters.analysisDIR + analysisVars), 'w') as shiny:
			shiny.write(" Rscript %s/snp2sim_analysis/vinaAnalysis/weighted_app.R %s/weighted_vinaSummary_%s.txt %s" %(parameters.programDIR, parameters.analysisDIR + analysisVars, parameters.protein, parameters.mulTrials))
		os.system("chmod +x %s/interactive_visualize.sh" %(parameters.analysisDIR + analysisVars))
		parameters.logger.info("Interactive vizualization tool: %s/interactive_visualize.sh" %(parameters.analysisDIR, analysisVars))
	else:
		curAnalysisDir = "%s/figures" % \
												(parameters.analysisDIR + analysisVars)
		if not os.path.isdir(curAnalysisDir):
			os.makedirs(curAnalysisDir)
		analysisCommand = "Rscript %s/snp2sim_analysis/vinaAnalysis/visualizations.R %s/vinaSummary_%s.txt %s" % (parameters.programDIR, 
													parameters.analysisDIR + analysisVars, parameters.protein, parameters.mulTrials)
		parameters.logger.info("Running analysis script")
		parameters.logger.debug("Executing command: %s", analysisCommand)
		try:
			r_out = subprocess.check_output(analysisCommand, shell = True, stderr = subprocess.STDOUT)
			parameters.logger.debug("Analysis output: \n%s", r_out.decode("utf-8"))
		except subprocess.CalledProcessError as e:
			parameters.logger.error("Clustering output: \n%s", e.output.decode('utf-8'))
			sys.exit(1)

		with open("%s/interactive_visualize.sh" %(parameters.analysisDIR + analysisVars), 'w') as shiny:
			shiny.write(" Rscript %s/snp2sim_analysis/vinaAnalysis/unweighted_app.R %s/vinaSummary_%s.txt %s" %(parameters.programDIR, parameters.analysisDIR + analysisVars, parameters.protein, parameters.mulTrials))
		os.system("chmod +x %s/interactive_visualize.sh" %(parameters.analysisDIR + analysisVars))
		parameters.logger.info("Interactive vizualization tool: %s/interactive_visualize.sh" %(parameters.analysisDIR + analysisVars))

def main():

	parameters = argParse()
	parameters.checkRequiredArgs()
	parameters.setDefault()

	

	if parameters.mode == "varMDsim":
		start = time.time()
		parameters.logger = logging.getLogger('varMDsim')
		parameters.logger.info("Starting varMDsim module")
		runVarMDsim(parameters)
		end = time.time()
		parameters.logger.info("MD simulation finished in %.3f sec!" %(end-start))
	elif parameters.mode == "varScaffold":
		start = time.time()
		parameters.logger = logging.getLogger('varScaffold')
		parameters.logger.info("Starting varScaffold module")
		runVarScaffold(parameters)    
		end = time.time()   
		parameters.logger.info("Scaffold clustering finished in %.3f sec!" %(end-start))
	elif parameters.mode == "drugSearch":
		start = time.time()
		parameters.logger = logging.getLogger('drugSearch')
		parameters.logger.info("Starting drugSearch module")
		runDrugSearch(parameters)
		end = time.time()
		parameters.logger.info("Drug Binding simulation finished in %.3f sec!" %(end-start))
	elif parameters.mode == "varAnalysis":
		start = time.time()
		parameters.logger = logging.getLogger('varAnalysis')
		parameters.logger.info("Starting varAnalysis module")
		runAnalysis(parameters)
		end = time.time()
		parameters.logger.info("Analysis completed in %.3f sec!" %(end-start))
	else:
		parameters.logger.error("Invalid mode selected")
	
main()
#--protein
#  name used to refer to group of variants
#
#--variant
#  variant to simulate
#
#Mode specific options:
#varTraj
#  --simulation time
#
#mdScaffold
#  --alignmentResID
#  --clusterResID
#  --clusterThresh
#
#drugSearch
#  --searchSpace
#  --library
#
#varAnalysis
#
#--------
#additional tools (in development)
#  buildSearchSpace - defines protein orientation and search parameters
#  buildDrugLibrary - creates new library from set of small molecules pdb files


