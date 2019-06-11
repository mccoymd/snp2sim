#!/usr/bin/env python

import os
import sys
import argparse
import copy
import random
import multiprocessing
import yaml
import shutil

#now unneccesary because mdtraj not used
#import mdtraj as md
#import numpy as np

class argParse():
	def __init__(self):
		self.requiredArgs = ['protein', 'mode']
		self.commandargs = _parseCommandLine()
		self.__dict__.update(self.commandargs.__dict__)
		#self.args = yaml.load(open('config.yaml'))
		#self.__dict__.update(self.args)
	def checkRequiredArgs(self):
		for arg in self.requiredArgs:
			assert (hasattr(self, arg) and getattr(self, arg)), arg + " not specified! " + arg + " is required."

		if self.mode == "varMDsim":
			if not self.simID:
				self.simID = str(random.randint(1,1000))
		if self.mode == "varScaffold":
			assert (hasattr(self, "scaffID") and getattr(self, "scaffID")), "no scaff ID - exiting"

	def setDefault(self):
		#self.runDIR = os.path.abspath(__file__)
		#self.runDIR = os.path.dirname(self.runDIR)
		self.runDIR = "/opt/snp2sim_results"
		if not os.path.isdir(self.runDIR):
				os.makedirs(self.runDIR)
		print(self.runDIR)

		print(self.protein)
		if not self.VMDpath:
			self.VMDpath = "vmd"
		if not self.NAMDpath:
			self.NAMDpath = "namd2"
		if not self.PYTHONSHpath:
		   #make sure pythonsh (from AutoDockTools) has been added to path
			self.PYTHONSHpath = "pythonsh"
		if not self.ADTpath:
			#must have path for "prepare_xxx.py" scripts from autodock tools
			self.ADTpath = "/opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
		if not self.VINApath:
			self.VINApath = "vina"

		if self.varResID and self.varAA:
			self.variant = self.varResID + self.varAA
		else:
			print("varResID and varAA not specified")
			print("Using WT structure for simulation")
			self.variant = "wt"



			
		self.simTopology = ("%s/simParameters/top_all36_prot.rtf" % self.runDIR,
								  "%s/simParameters/toppar_water_ions_namd.str" % self.runDIR)
		self.simParameters = ("%s/simParameters/par_all36_prot.prm" % self.runDIR,
									"%s/simParameters/toppar_water_ions_namd.str" %self.runDIR)

		if not self.simProc:
			self.simProc = multiprocessing.cpu_count()    
		
		self.simPDB = "%s/variantSimulations/%s/structures/%s.%s.pdb" \
							% (self.runDIR,self.protein,
							   self.protein,self.variant)
		self.simPSF = "%s/variantSimulations/%s/structures/%s.%s.psf" \
							% (self.runDIR,self.protein,
							   self.protein,self.variant)
		self.templatePDB = "%s/variantSimulations/%s/structures/%s.template.pdb" \
								 % (self.runDIR,self.protein,
									self.protein)
		self.wtPDB = "%s/variantSimulations/%s/structures/%s.wt.UNSOLVATED.pdb" \
						   % (self.runDIR,self.protein,
							  self.protein)
		self.wtPSF = "%s/variantSimulations/%s/structures/%s.wt.UNSOLVATED.psf" \
						   % (self.runDIR,self.protein,
							  self.protein)
		self.varPrefix = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED" \
							   % (self.runDIR,self.protein,
								  self.protein, self.variant)
		self.varPDB = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED.pdb" \
							% (self.runDIR, self.protein,
							   self.protein, self.variant)
		self.varPSF = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED.psf" \
							% (self.runDIR,self.protein,
							   self.protein, self.variant)
		self.wtStructTCL = "%s/variantSimulations/%s/bin/%s.wt.genStructFiles.tcl" \
								 % (self.runDIR,self.protein, self.protein)
		self.varStructTCL = "%s/variantSimulations/%s/bin/%s.%s.genStructFiles.tcl" \
								  % (self.runDIR,self.protein,
									 self.protein, self.variant)
		self.singleRunTCL = "%s/variantSimulations/%s/bin/%s.%s.genSingleRunPDB.tcl" \
								  % (self.runDIR,self.protein,
									 self.protein, self.variant)
		self.solvTCL = "%s/variantSimulations/%s/bin/%s.%s.genSolvStruct.tcl" \
							 % (self.runDIR,self.protein,
								self.protein, self.variant)
		self.solvBoundary = "%s/variantSimulations/%s/config/%s.%s.solvBoundary.txt" \
								  % (self.runDIR,self.protein,
									 self.protein, self.variant)
		self.structPrefix = "%s/variantSimulations/%s/structures/%s.%s" \
								  % (self.runDIR,self.protein,
									 self.protein, self.variant)
		self.NAMDconfig = "%s/variantSimulations/%s/config/%s.%s.%s.NAMD" \
								% (self.runDIR, self.protein, self.protein,
								   self.variant, self.simID)
		self.NAMDout = "%s/variantSimulations/%s/results/%s/trajectory/%s.%s.%s" \
							 % (self.runDIR, self.protein, self.variant,
								self.protein, self.variant, self.simID)
		self.scaffParams = "%s/variantSimulations/%s/config/%s.scaff" \
							 % (self.runDIR, self.protein, self.scaffID)
		self.resultsDIR =  "%s/variantSimulations/%s/results" % \
								 (self.runDIR, self.protein)
		self.trajDIR = "%s/variantSimulations/%s/results/%s/trajectory/" \
							 % (self.runDIR, self.protein, self.variant)

		self.drugLibPath = "%s/drugLibraries/%s/" % \
				   (self.runDIR, self.drugLibrary)

                
		#hardcoding bindingID as drugLibrary
		#TODO - refactor to remove "bindingID"
		self.bindingID = self.drugLibrary
		
		if self.bindingID:
			self.drugBindConfig = "%s/variantSimulations/%s/config/%s.autodock" \
										% (self.runDIR, self.protein, self.bindingID)
		

		if self.loadPDBtraj:
			variantDIR = self.resultsDIR + "/" + self.variant + "/trajectory/"
			if not os.path.exists(variantDIR):
				os.makedirs(variantDIR)

			for pdbFile in self.loadPDBtraj:
				pdbNewLoc = self.trajDIR + os.path.basename(pdbFile)
				print(pdbNewLoc)
				os.system("cp %s %s" % (pdbFile,pdbNewLoc))
		
		if self.newScaff:
			if not os.path.isdir("%s/variantSimulations/%s/config" % \
								 (self.runDIR,self.protein)):
				os.makedirs("%s/variantSimulations/%s/config" % \
							(self.runDIR,self.protein))
			if not os.path.isfile(self.scaffParams):
				os.system("cp %s %s/variantSimulations/%s/config/%s.scaff" \
						  % (self.newScaff,self.runDIR,
							 self.protein, self.scaffID))

			else:
				print("Scaff config for %s exists." % (self.scaffParams))
				print("Select new scaffID or remove existing scaff config")
				sys.exit()


def _parseCommandLine():

	parser = argparse.ArgumentParser(
		#version="0.1",
		description="snp2sim - Molecular Simulation of Somatic Variation",
		epilog="written by Matthew McCoy, mdm299@georgetown.edu"
	)

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
	parser.add_argument("--newScaff",
						help="path to scaffold config (search RMSD and # iterations, alignment aa, clusteraa",
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
	parser.add_argument("--vinaExh",
						help="exhaustiveness parameter of autodock vina",
						action="store",
						default="50",
						type=str,
						)
	parser.add_argument("--clustThresh",
						help="exhaustiveness parameter of autodock vina",
						action="store",
						default="0.09",
						type=float,
						)
	parser.add_argument("--simID",
						help="amino acid to change",
						action="store",
						type=str,
						)
	parser.add_argument("--bindingID",
						help="name for autodock template config (search space and alignment res)",
						action="store",
						type=str,
						)
	parser.add_argument("--scaffID",
						help="name of scaffolding parameters set",
						action="store",
						type=str,
						)
	parser.add_argument("--VMDpath",
						help="path to VMD executable",
						action="store",
						type=str,
						)
	parser.add_argument("--VINApath",
						help="path to AutoDock Vina executable",
						action="store",
						type=str,
						)
	parser.add_argument("--NAMDpath",
						help="path to NAMD executable",
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

	parser.add_argument("--simProc",
						help="number of processors to run simulation",
						action="store",
						type=int,
						)
	parser.add_argument("--singleRun",
						help="output summary PDB trajectory only",
						action="store_true",
						)
	parser.add_argument("--clustPDBtraj",
						help="cluster results from multiple varScafold singleRun",
						action="store_true",
						)
	parser.add_argument("--loadPDBtraj", nargs='+',
						help="pdb trajectory files to import, list one after another",
						action="store",
						)
	parser.add_argument("--inputScaff", nargs='+',
						help="pdb scaffold files to import, list one after another",
						action="store",
						)
	parser.add_argument("--cgcRun",
                            help="if Run is on CGC platform, move pdb and log to snp2sim root for processing",
			    action="store_true",
        )
	parser.add_argument("--bindSingleVar",
						help="only bind single variant scaffolds",
						action="store_true",
						)
	parser.add_argument("--clean",
						help="remove all files from previous run with same protein name",
						action="store_true",
						)
	

	cmdlineparameters, unknownparams = parser.parse_known_args()
	parameters = copy.copy(cmdlineparameters)

	return parameters

def genNAMDstructFiles(parameters):
	if not os.path.exists("%s/variantSimulations/%s/bin/" % (parameters.runDIR,parameters.protein)):
		os.makedirs("%s/variantSimulations/%s/bin/" % (parameters.runDIR,parameters.protein))
	if os.path.isfile(parameters.templatePDB):
		print("generating solvated/ionized PDB and PSF from %s" \
			% parameters.templatePDB)

		if os.path.isfile(parameters.varPDB) and os.path.isfile(parameters.varPSF):
			print("unsolvated PDB and PSF exist")
		else:
			genStructTCL(parameters)
			genStructCommand = "%s -e %s" % (parameters.VMDpath, parameters.varStructTCL)
			os.system(genStructCommand)

		genSolvTCL(parameters)
		genSolvCommand = "%s -e %s" % (parameters.VMDpath, parameters.solvTCL)
		if not os.path.exists("%s/variantSimulations/%s/config/" % \
							  (parameters.runDIR,parameters.protein)):
			os.makedirs("%s/variantSimulations/%s/config/" % (parameters.runDIR,parameters.protein))
		os.system(genSolvCommand)
	else:
		print("no template exists")
		print("use --new to specify clean PDB (only cannonical aa) ")
		sys.exit()


def genStructTCL(parameters):
	print("generating struct files")
	structFile = open(parameters.varStructTCL,"w+")
	
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
		structFile.write("mutator -psf %s -pdb %s -o %s -ressegname PROT -resid %s -mut %s\n" \
						 % (parameters.wtPSF, parameters.wtPDB, \
							parameters.varPrefix, parameters.varResID, longAA.get(parameters.varAA)))
	
	structFile.write("quit\n")


def genSolvTCL(parameters):
	#todo - adjustable solvation/ionization parameters
	print("generating solvated/ionized pdb and psf files")

	solvFile = open(parameters.solvTCL,"w+")
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

def runNAMD(parameters):
	if not os.path.exists("%s/variantSimulations/%s/config/" % (parameters.runDIR,parameters.protein)):
		os.makedirs("./variantSimulations/%s/config/" % (parameters.runDIR,parameters.protein))


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
		print("boundary file does not exist")
		print("remove solvated/ionized PDB and PSF and resubmit")
		sys.exit()
		
	configFile = open(parameters.NAMDconfig,"w+")
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
					 
	if not os.path.isdir("%s/variantSimulations/%s/results/%s/trajectory/" \
						 % (parameters.runDIR, parameters.protein, parameters.variant)):
		os.makedirs("%s/variantSimulations/%s/results/%s/trajectory/" \
						 % (parameters.runDIR, parameters.protein, parameters.variant))

	print("running NAMD with %i processors" % parameters.simProc)

#moving to runVarMDsim
#        runNAMDcommand = "%s +p%i %s > %s.log" % \
#			 (parameters.NAMDpath, parameters.simProc,
#			  parameters.NAMDconfig, parameters.NAMDout)
#	os.system(runNAMDcommand)
                

def genSingleRunTCL(parameters):
	print("generating %s" % parameters.singleRunTCL)
	pdbTCL = open(parameters.singleRunTCL,"w+")
	pdbTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
	print("using DCD files in %s" % variantDIR)
	for tFile in os.listdir(variantDIR):
		if tFile.endswith(".dcd"):
			dcdFile = variantDIR + tFile
			pdbTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	pdbTCL.write("package require pbctools\n")
	pdbTCL.write("pbc wrap -center com -centersel \"protein\" -compound residue -all\n")

	allPDBoutput = "%s/%s.%s.%s.pdb" % (variantDIR, parameters.protein,
										parameters.variant, parameters.simID)
	pdbTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
				   % allPDBoutput)
	pdbTCL.write("quit")
	
	genPDBcommand = "%s -e %s" % (parameters.VMDpath, parameters.singleRunTCL)
	
	os.system(genPDBcommand)

def runClusterTCL(parameters):
#TODO - currently config must be formated as:
#NOTE:: 3 lines, CASE SPECIFIC, NEEDS IMPROVEMENT
#alignmentRes "(atomselection)"
#clusterRes "(atomselection)"
#rmsdThresh (num)
	
	scaffParameters = open(parameters.scaffParams,"r")
	scaffLines = scaffParameters.readlines()

	alignmentRes = scaffLines[0]
	alignmentRes = alignmentRes.rstrip()
	alignmentRes = alignmentRes.replace("alignmentRes ","")

	clusterRes = scaffLines[1]
	clusterRes = clusterRes.rstrip()
	clusterRes = clusterRes.replace("clusterRes ","")

	rmsdThresh = scaffLines[2]
	rmsdThresh = rmsdThresh.rstrip()
	rmsdThresh = rmsdThresh.replace("rmsdThresh ","")


	#regenerate config even if one exists to ensure new trajecories are included
	
	if not os.path.exists("%s/variantSimulations/%s/results/%s/scaffold" % \
						  (parameters.runDIR,parameters.protein,parameters.variant)):
		os.makedirs("%s/variantSimulations/%s/results/%s/scaffold" % \
					(parameters.runDIR,parameters.protein,parameters.variant))

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	print(scaffLOG)
	if os.path.isfile(scaffLOG):
		print("%s exists. Remove to recalculate Scaffolds" % scaffLOG)
		return
		
	print("generating %s" % parameters.scaffoldTCL)
	clustTCL = open(parameters.scaffoldTCL,"w+")
	clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
	variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
	print("using DCD files in %s" % variantDIR)
	for tFile in os.listdir(variantDIR):
		if tFile.endswith(".dcd"):
			dcdFile = variantDIR + tFile
			clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

	clustTCL.write("package require pbctools\n")
	clustTCL.write("pbc wrap -center com -centersel \"protein\" -compound residue -all\n")
	clustTCL.write("set nf [molinfo top get numframes]\n")
	clustTCL.write("set refRes [atomselect top %s frame 0]\n" % alignmentRes)
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top %s frame $i]\n" % alignmentRes)
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")
	clustTCL.write("set output [open %s w]\n" % scaffLOG)
	
	clustTCL.write("set clustRes [atomselect top %s]\n" % clusterRes)
	clustTCL.write("puts $output \"clusters for RMSD Threshold %s\"\n" % rmsdThresh)
	clustTCL.write("puts $output [measure cluster $clustRes distfunc rmsd cutoff %s]\n" \
				   % rmsdThresh)
	clustTCL.write("puts $output \"Only generated summary PDB file\"\n")
	clustTCL.write("close $output\n")
	allPDBoutput = parameters.scaffBASE + ".all.pdb"
	clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
				   % allPDBoutput)
	clustTCL.write("quit")

#legacy code for surveying the trajectoies
#    for {set i 60} {$i < 80} {incr i} {
#                set thr [expr $i / 100.0]
#                puts $output [puts $thr]
#                puts $output [measure cluster $clustRes distfunc rmsd cutoff $thr]
#            }
	
#	vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
#        print("got here %s" % vmdClustCommand)
#	os.system(vmdClustCommand)

def runPDBclustTCL(parameters):
	scaffParameters = open(parameters.scaffParams,"r")
	scaffLines = scaffParameters.readlines()

	alignmentRes = scaffLines[0]
	alignmentRes = alignmentRes.rstrip()
	alignmentRes = alignmentRes.replace("alignmentRes ","")

	clusterRes = scaffLines[1]
	clusterRes = clusterRes.rstrip()
	clusterRes = clusterRes.replace("clusterRes ","")

	rmsdThresh = scaffLines[2]
	rmsdThresh = rmsdThresh.rstrip()
	rmsdThresh = rmsdThresh.replace("rmsdThresh ","")


	#regenerate config even if one exists to ensure new trajecories are included
	scaffPath = "%s/variantSimulations/%s/results/%s/scaffold" % \
						  (parameters.runDIR,parameters.protein,parameters.variant)
	if not os.path.exists(scaffPath):
		os.makedirs(scaffPath)

	#will not generate new TCL if previous scaffID log exists
	scaffLOG = parameters.scaffBASE + ".log"
	print(scaffLOG)
	if os.path.isfile(scaffLOG):
		print("%s exists. Remove to recalculate Scaffolds" % scaffLOG)
		return

	print("generating %s" % parameters.scaffoldTCL)

	binPath = "%s/variantSimulations/%s/bin/" % \
						  (parameters.runDIR,parameters.protein)
	if not os.path.exists(binPath):
		os.makedirs(binPath)

			
	clustTCL = open(parameters.scaffoldTCL,"w+")
	variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
	if not os.path.exists(variantDIR):
		os.makedirs(variantDIR)
		
	for trajFile in os.listdir(variantDIR):
		if trajFile.endswith(".pdb"):
			pdbFile = variantDIR + "/" + trajFile
			clustTCL.write("mol addfile %s waitfor all\n" % pdbFile)
	clustTCL.write("set nf [molinfo top get numframes]\n")
	clustTCL.write("set refRes [atomselect top %s frame 0]\n" % alignmentRes)
	clustTCL.write("set refStruct [atomselect top all frame 0]\n")
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("  set curStruct [atomselect top all frame $i]\n")
	clustTCL.write("  set curRes [atomselect top %s frame $i]\n" % alignmentRes)
	clustTCL.write("  set M [measure fit $curRes $refRes]\n")
	clustTCL.write("  $curStruct move $M\n")
	clustTCL.write("}\n")
	clustTCL.write("set output [open %s w]\n" % scaffLOG)
	clustTCL.write("set clustRes [atomselect top %s]\n" % clusterRes)
	clustTCL.write("puts $output \"clusters for RMSD Threshold %s\"\n" % rmsdThresh)
	clustTCL.write("puts $output [measure cluster $clustRes distfunc rmsd cutoff %s]\n" \
				   % rmsdThresh)
	clustTCL.write("close $output\n")


	allPDBoutput = parameters.scaffBASE + ".all.pdb"
	clustTCL.write("animate write pdb %s beg 0 end -1 sel [atomselect top \"protein\"]\n" \
				   % allPDBoutput)
	clustTCL.write("quit")

#legacy code for surveying the trajectoies
#    for {set i 60} {$i < 80} {incr i} {
#                set thr [expr $i / 100.0]
#                puts $output [puts $thr]
#                puts $output [measure cluster $clustRes distfunc rmsd cutoff $thr]
#            }

# moving this to runVarScaffold
#	vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
#	os.system(vmdClustCommand)  

def sortPDBclusters(parameters):
	allPDBoutput = parameters.scaffBASE + ".all.pdb"
	#traj = md.load(allPDBoutput)

	allPDBstruct = open(allPDBoutput, "r")
	allPDBlines = allPDBstruct.readlines()
	pdbHeader = allPDBlines.pop(0)
	indPDB = []
	currPDB = ""
	for pdbLine in allPDBlines:
		currPDB = currPDB + pdbLine
		if "END" in pdbLine:
			indPDB.append(currPDB)
			currPDB = ""
	
	scaffLOG = parameters.scaffBASE + ".log"    
	clustLogfile = open(scaffLOG, "r")
	clustLogLines = clustLogfile.readlines()
	clusterMembership = clustLogLines[1]
	clusterMembership = clusterMembership.rstrip()
	clusterMembership = clusterMembership.replace("} {","_")
	clusterMembership = clusterMembership.replace("{","")
	clusterMembership = clusterMembership.replace("}","")
	clusterMembership = clusterMembership.split("_")
	#remove "unclustered struct"
	del clusterMembership[-1]
	scaffNum = 1;
	#print(traj.n_frames)
	print(len(indPDB))
	for indCluster in clusterMembership:
		repStructIndex = int(indCluster.split(" ")[0])
		#if len(indCluster.split(" ")) > parameters.clustThresh*traj.n_frames:
		#    repStruct = traj[repStructIndex]
		#    scaffFileName = parameters.scaffBASE + ".cl" + str(scaffNum) + ".scaffold.pdb"
		#    repStruct.save(scaffFileName)
		#    scaffNum += 1 
		#print str(len(indCluster)) + " cluster " + str(scaffNum)
		if len(indCluster.split(" ")) > parameters.clustThresh*len(indPDB):
			#scaffFileName = parameters.scaffBASE + ".cluster" + str(scaffNum) + ".pdb"
			scaffFileName = parameters.scaffBASE + ".cl" + str(scaffNum) + ".scaffold.pdb"
			scaffFile = open(scaffFileName,"w+")
			scaffFile.write(pdbHeader)
			scaffFile.write(indPDB[int(repStructIndex)])
			scaffFile.close()
			scaffNum += 1
	#        for structID in indCluster:
	#            print structID
	#            if structID:
	#                scaffFile.write(indPDB[int(structID)])
	#    scaffNum += 1 


def genScaffoldTCL(parameters):
	
	scaffParameters = open(parameters.scaffParams,"r")
	scaffLines = scaffParameters.readlines()

	alignmentRes = scaffLines[0]
	alignmentRes = alignmentRes.rstrip()
	alignmentRes = alignmentRes.replace("alignmentRes ","")

	scaffDIR = parameters.resultsDIR + "/" + parameters.variant + "/scaffold/"
	genScaff = open(parameters.clusterTCL, "w+")
	for tFile in os.listdir(scaffDIR):
		if tFile.endswith(".pdb"):
			if "cluster" in tFile:
				pdbClustFile = scaffDIR + tFile
				pdbScaffFile = pdbClustFile
				pdbScaffFile = pdbScaffFile.replace(".pdb",".scaffold.pdb")
				genScaff.write("mol new %s waitfor all\n" % pdbClustFile)
#                genScaff.write("set domain [atomselect top %s]\n" % alignmentRes)
				genScaff.write("set domain [atomselect top all]\n")
				genScaff.write("set avePos [measure avpos $domain]\n")
				genScaff.write("$domain set {x y z} $avePos\n")
				genScaff.write("$domain writepdb %s\n" % pdbScaffFile)
				
	genScaff.write("quit\n")

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
				

def parseADconfig(parameters):
	paramFile = "%s/variantSimulations/%s/config/%s.autodock" % \
						(parameters.runDIR,parameters.protein,parameters.bindingID)
	autodockParams = open(paramFile, "r")
	paramData = autodockParams.readlines()
	parameters.ADsearchSpace = [paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip(),
								paramData.pop(0).rstrip()]

#Function not used?
# def parseFlexConfig(parameters):
#     #TODO input from config folder
#     parameters.flexRes = paramData.pop(0)
#     parameters.flexRes = parameters.flexRes.rstrip()
#     parameters.flexRes = parameters.flexRes.split(" ")


def getFlexRes(pdbFile,flexRes):
	flexConfig = open(flexRes, "r").readlines()
	flexRes = flexConfig.pop(0)
	flexRes = flexRes.rstrip()
	flexRes = flexRes.split(" ")
	
	
	currScaff = open(pdbFile,"r").readlines()
	currScaff.pop(0)
	flexRes = [int(x) for x in flexRes]
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
	if parameters.scaffID:
		scaffParameters = open(parameters.scaffParams,"r")
	else:
		print("no scaff params specified")
		sys.exit()
		
	scaffLines = scaffParameters.readlines()

	clusterRes = scaffLines[1]
	clusterRes = clusterRes.rstrip()
	clusterRes = clusterRes.replace("clusterRes ","")

	if not os.path.isdir("%s/variantSimulations/%s/bin/" \
						 % (parameters.runDIR, parameters.protein)):
		os.makedirs("%s/variantSimulations/%s/bin/" \
					 % (parameters.runDIR, parameters.protein))

	
	alignmentConfig = "%s/variantSimulations/%s/bin/%s.%s.alignStuct.TCL" \
							% (parameters.runDIR, parameters.protein,
							   parameters.protein, parameters.variant)
	alignmentTCL = open(alignmentConfig, "w+")
	alignmentTCL.write("mol new %s\n" % parameters.templatePDB)
	alignmentTCL.write("set ref [atomselect top %s]\n" % clusterRes)
	alignmentTCL.write("mol new %s\n" % currScaff)
	alignmentTCL.write("set target [atomselect top %s]\n" % clusterRes)
	alignmentTCL.write("set allP [atomselect top all]\n")
	alignmentTCL.write("set M [measure fit $target $ref]\n")
	alignmentTCL.write("$allP move $M\n")
	alignmentTCL.write("$allP writepdb %s\n" % currScaff)
	alignmentTCL.write("quit\n")
	alignmentTCL.close()
	alignmentCommand = "%s -e %s" % (parameters.VMDpath, alignmentConfig)
	os.system(alignmentCommand)

def genVinaConfig(parameters):
	vinaConfig = open(parameters.vinaConfig,"w+")
	if not os.path.isfile(parameters.flexConfig):
		vinaConfig.write("receptor = %s\n" % parameters.scaff1out)
	else:
		vinaConfig.write("receptor = %s\n" % parameters.scaffRigid)
		vinaConfig.write("flex = %s\n" % parameters.scaffFlex)
	vinaConfig.write("ligand = %s\n" % parameters.currDrugPath)
	vinaConfig.write("out = %s/%s.pdbqt\n" % (parameters.vinaOutDir, parameters.vinaBase))
	vinaConfig.write("log = %s/%s.log\n" % (parameters.vinaOutDir, parameters.vinaBase))
	vinaConfig.write("exhaustiveness = %s\n" % parameters.vinaExh)
	for searchParam in parameters.ADsearchSpace:
		vinaConfig.write("%s\n" % searchParam)

#moving to runDrugSearch                
#	vinaCommand = "%s --config %s" % (parameters.VINApath, parameters.vinaConfig)
#	os.system(vinaCommand)





#### start of program move to main()
def runVarMDsim(parameters):
	print("Performing varMDsim")
	if parameters.newStruct:
		if os.path.isdir("%s/variantSimulations/%s" % (parameters.runDIR,parameters.protein)) and not parameters.clean:
			print("ERROR - %s is already created... resubmit with new protein name" %parameters.protein)
			sys.exit()
		if parameters.clean:
			shutil.rmtree('%s/variantSimulations/%s' %(parameters.runDIR,parameters.protein))
			os.makedirs("%s/variantSimulations/%s" % (parameters.runDIR,parameters.protein))
		os.makedirs("%s/variantSimulations/%s/structures" % (parameters.runDIR,parameters.protein))
		os.system("cp %s %s/variantSimulations/%s/structures/%s.template.pdb" \
				  % (parameters.newStruct,parameters.runDIR,
					 parameters.protein, parameters.protein)) 

	if os.path.isfile(parameters.simPDB):
		if os.path.isfile(parameters.simPSF):
			print("using %s %s as initial structure" \
				% (parameters.protein, parameters.variant))

			#if no variant structure exists, use clean template to create
		else:
			print("%s does not exist" % parameters.simPSF)
			print("generating new pdb and psf for simulation")
			genNAMDstructFiles(parameters)
	else:
		print("%s does not exist" % parameters.simPDB)
		print("generating new pdb and psf for simulation")
		genNAMDstructFiles(parameters)

	if parameters.simLength:
		print("Performing Variant %.3f ns Simulation" % parameters.simLength)
		
		runNAMD(parameters)
                runNAMDcommand = "%s +p%i %s > %s.log" % \
			         (parameters.NAMDpath, parameters.simProc,
			          parameters.NAMDconfig, parameters.NAMDout)
	        os.system(runNAMDcommand)
                if not os.path.isfile("%s.dcd" % parameters.NAMDout):
                        print("NAMD run failed")
                        sys.exit()
		
		if parameters.singleRun:
			genSingleRunTCL(parameters)
			if parameters.cgcRun:
				CGCpdb = parameters.NAMDout + ".pdb"
				CGClog = parameters.NAMDout + ".log"
				cwd = os.getcwd()
				os.system("mv %s %s" % (CGCpdb, cwd))
				os.system("mv %s %s" % (CGClog, cwd))

	else:
		print("simulation length not specified")        
		sys.exit()

def runVarScaffold(parameters):
	print("Performing varScaffold")
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

	if os.path.isfile(parameters.scaffParams):
		#check if variant specified, otherwise perform clustering for all variants
		if parameters.variant:
			print("generating Scaffold for %s ONLY" % (parameters.variant))
			variantList = (parameters.variant,)
		else:
			print("generating all Variant Scaffolds")
			variantList = os.listdir(parameters.resultsDIR)

## TODO include option to analyze trajectory clusters to determine rmsd threshold
		for varSimResult in variantList:
			parameters.variant = varSimResult
			parameters.simPSF = "%s/variantSimulations/%s/structures/%s.%s.psf" \
								% (parameters.runDIR,parameters.protein,
								   parameters.protein,parameters.variant)
			parameters.scaffoldTCL =  "%s/variantSimulations/%s/bin/%s.%s.%s.genScaffold.tcl" % \
									  (parameters.runDIR, parameters.protein,
									   parameters.protein, parameters.variant, parameters.scaffID)
			parameters.trajAnalysisTCL =  "%s/variantSimulations/%s/bin/%s.%s.%s.analyzeTraj.tcl" % \
										  (parameters.runDIR, parameters.protein,
										   parameters.protein, parameters.variant, parameters.scaffID)
			parameters.clusterTCL =  "%s/variantSimulations/%s/bin/%s.%s.%s.genRepScaffold.tcl" % \
										  (parameters.runDIR, parameters.protein,
										   parameters.protein, parameters.variant, parameters.scaffID)
			parameters.scaffBASE = "%s/variantSimulations/%s/results/%s/scaffold/%s.%s.%s" % \
								   (parameters.runDIR, parameters.protein,parameters.variant,
									parameters.protein, parameters.variant, parameters.scaffID)


			print("generating Scaffold for %s" % (parameters.variant))
			#will not generate new TCL if previous scaffID log exists
			scaffLOG = parameters.scaffBASE + ".log"
			print(scaffLOG)
			if not os.path.isfile(scaffLOG):
				if not parameters.clustPDBtraj:
					runClusterTCL(parameters)
				else:
					parameters.scaffoldTCL =  "%s/variantSimulations/%s/bin/%s.%s.%s.genPDBtrajScaffold.tcl" % \
											  (parameters.runDIR, parameters.protein,
											   parameters.protein, parameters.variant, parameters.scaffID)
					runPDBclustTCL(parameters)
	                        vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	                        os.system(vmdClustCommand)  

					

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
				print("%s exists. Remove to regenerate scaffold stats" % scaffLOG)

	else:
		print("Scaffold Parameters %s does not exist" % (parameters.scaffParams))
		print("Exiting")


def runDrugSearch(parameters):

	print("Performing drugSearch")
	if parameters.newBindingConfig:
		if not os.path.isfile(parameters.drugBindConfig):            
			print("using new config %s" % parameters.newBindingConfig)
			configDIR = "%s/variantSimulations/%s/config" % \
						(parameters.runDIR,parameters.protein)
			if not os.path.isdir(configDIR):
				os.makedirs(configDIR)
				
			os.system("cp %s %s/%s.autodock" % (parameters.newBindingConfig, configDIR,
												parameters.bindingID))
		else:
			print("%s already exists - remove or choose new bindingID" % parameters.drugBindConfig)
			sys.exit()

	if os.path.isfile(parameters.drugBindConfig):
		parseADconfig(parameters)
	else:
		print("%s does not exist" % parameters.drugBindConfig)
		sys.exit()

	parameters.vinaOutDir = "%s/variantSimulations/%s/results/%s/drugBinding/" % \
							(parameters.runDIR,
							 parameters.protein,
							 parameters.variant)
	if not os.path.isdir(parameters.vinaOutDir):
		os.makedirs(parameters.vinaOutDir)

		
	parameters.flexConfig = "%s/variantSimulations/%s/config/%s.flex" % \
							(parameters.runDIR,parameters.protein,parameters.bindingID)

	if parameters.flexBinding:
		if not os.path.isdir("%s/variantSimulations/%s/config" % \
							 (parameters.runDIR,parameters.protein)):
			os.makedirs("%s/variantSimulations/%s/config" % \
						(parameters.runDIR,parameters.protein))
			
		
		if not os.path.isfile(parameters.flexConfig):
			os.system("cp %s %s" \
					  % (parameters.flexBinding, parameters.flexConfig))
		else:
			print("%s already exists - remove or choose new bindingID" % parameters.flexConfig)
			sys.exit()

	if parameters.inputScaff:
		print("using input scaffolds")
		variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/scaffold/"
		#if not os.path.exists(variantDIR):
		#    os.makedirs(variantDIR)
		os.makedirs(variantDIR)
		inputCount = 1
		for pdbFile in parameters.inputScaff:
#            pdbNewLoc = variantDIR + os.path.basename(pdbFile)
#            pdbNewLoc = os.path.splitext(pdbNewLoc)[0] + ".scaffold.pdb"
			inputID = "input%s" % str(inputCount)
			pdbNewLoc = "%s/%s.%s.%s.scaffold.pdb" % (variantDIR, parameters.protein,
													  parameters.variant, inputID)
			inputCount += 1
#            print pdbNewLoc
			os.system("cp %s %s" % (pdbFile,pdbNewLoc))


	if parameters.drugLibrary:
		if parameters.singleDrug:
			#todo - add single drug binding
			print("binding single drug %s" % parameters.singleDrug)
                        if not os.path.isdir(parameters.drugLibPath):
			        print("drug library does not exist! Creating it!")
			        os.makedirs(parameters.drugLibPath)
                        os.system("cp %s %s" % (parameters.singleDrug,
                                                parameters.drugLibPath))
			print("using small molecules in %s" % parameters.drugLibPath)

		else:
		        if not os.path.isdir(parameters.drugLibPath):
			        print("drug library does not exist ")
			        sys.exit()
		        else:
			        print("using small molecules in %s" % parameters.drugLibPath)

			
							

			
	if parameters.bindingID:
		print("using search space defined in %s" % parameters.drugBindConfig)


		if not os.path.isfile(parameters.templatePDB):
			if parameters.bindingTemplate:
				if not os.path.isdir("%s/variantSimulations/%s/structures/" %\
								 (parameters.runDIR,parameters.protein)):
					os.makedirs("%s/variantSimulations/%s/structures" % \
								(parameters.runDIR,parameters.protein))
				os.system("cp %s %s/variantSimulations/%s/structures/%s.template.pdb" \
						  % (parameters.bindingTemplate, parameters.runDIR,
							 parameters.protein, parameters.protein))
			else:
				print("specify binding template")
				sys.exit()

	   
				
		if parameters.bindSingleVar:
			bindingVar = [parameters.variant,]
		else:
			bindingVar = os.listdir(parameters.resultsDIR)
			
		for var in bindingVar:
			#regenerate pdbqt for all scaffold.pdb files
			print("binding to variant %s" % var)
			#todo check if pdbqt exist already
			scaffDIR = "%s/%s/scaffold/" % (parameters.resultsDIR, var)
			for scaffPDB in os.listdir(scaffDIR):
				if scaffPDB.endswith("scaffold.pdb"):
					print(scaffPDB)
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
					os.system(prepBaseScaff)

					if os.path.isfile(parameters.flexConfig):
						print("Found Flex Config: performing flexable residue binding")
						
						scaffFlexRes = getFlexRes(currScaffPath,parameters.flexConfig)
						print(scaffFlexRes)
						parameters.scaffFlex = scaffBase + ".flex.pdbqt"
						parameters.scaffRigid = scaffBase + ".rigid.pdbqt"
						prepFlexScaff = "%s %s/prepare_flexreceptor4.py -r %s -s %s -g %s -x %s" \
									 % (parameters.PYTHONSHpath, parameters.ADTpath,
										parameters.scaff1out, scaffFlexRes,
										parameters.scaffRigid, parameters.scaffFlex)
						print(prepFlexScaff)
						os.system(prepFlexScaff)

					for dFile in os.listdir(parameters.drugLibPath):
						if dFile.endswith(".pdbqt"):
							parameters.currDrugPath = parameters.drugLibPath + dFile
							parameters.drugBase = os.path.splitext(dFile)[0]
							parameters.vinaBase = os.path.basename("%s.%s.%s" % \
																   (scaffBase,parameters.bindingID,
																	parameters.drugBase))
							parameters.vinaConfig = "%s/variantSimulations/%s/config/%s.vina" % \
													(parameters.runDIR,
													 parameters.protein,
													 parameters.vinaBase)

							genVinaConfig(parameters)
                                                        vinaCommand = "%s --config %s" % (parameters.VINApath, parameters.vinaConfig)
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
		print("specify bindingID")
		sys.exit()
		

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

def main():

	parameters = argParse()
	parameters.checkRequiredArgs()
	parameters.setDefault()

	if parameters.mode == "varMDsim":
		runVarMDsim(parameters)
		print("MD Simulation finished!")
	elif parameters.mode == "varScaffold":
		runVarScaffold(parameters)       
		print("Scaffold clustering finished!")
	elif parameters.mode == "drugSearch":
		runDrugSearch(parameters)
		print("Drug binding simulation finished!")
	else:
		print("invalid mode selected")
	
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


