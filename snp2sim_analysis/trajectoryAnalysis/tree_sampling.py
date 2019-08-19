import os
import sys
import argparse
import subprocess
import shutil
import yaml

class samplingTree():
	def __init__(self, parameters, run, scaff, parent):
		self.num = run
		self.scaff = scaff
		self.child = []
		self.parent = parent
		self.parameters = parameters
		if self.parent:
			self.depth = self.parent.depth + 1
		else:
			self.depth = 0
		if self.depth < self.parameters.depth:
			self.extend = True
		else:
			self.extend = False
	def addChild(self, node):
		self.child.append(node)
	def mark(self):
		self.extend = False
	def __str__(self):
		if not self.child:
			return "[%d]" % self.num
		else:
			return "[%d: %s]" % (self.num, ", ".join([str(x) for x in self.child]))

class argParse():
	def __init__(self):
		parser = argparse.ArgumentParser(
		#version="0.1",
		description="snp2sim - Molecular Simulation of Somatic Variation",
		epilog="written by Matthew McCoy, mdm299@georgetown.edu, and Vikram Shivakumar"
		)

		parser.add_argument("--config",
							help="use parameters from config.yaml. OVERWRITES command line parameters",
							action="store",
							)
		parser.add_argument("--depth",
							help="maximum depth to search for scaffolds",
							action="store",
							type=int,
							)
		parser.add_argument("--newStruct",
							help="initial structure",
							action="store",
							)
		parser.add_argument("--initialTraj",
							help="initial trajectory (optional)",
							action="store",
							)
		parameters, unknownparams = parser.parse_known_args()

		if (parameters.config):
				self.__dict__.update(yaml.load(open(parameters.config)))
		else:
			print("No config file given.")
			sys.exit(1)

		self.__dict__.update(parameters.__dict__)

	def setDefaults(self):
		self.programDir = os.path.dirname(os.path.realpath(__file__))
		if "." in self.protein:
			self.protein.replace(".", "_")
		if isinstance(self.varAA,list) and isinstance(self.varResID,list):
			if len(self.varAA) == len(self.varResID):
				self.variant = [str(self.varResID[x]) + self.varAA[x] for x in range(len(self.varAA))]
				self.variant = "_".join(self.variant)
			else:
				print("number of ResIDs does not match number of variant AAs.")
				print("Using WT structure for simulation")
				self.variant = "wt"
		elif self.varResID and self.varAA:
			self.variant = str(self.varResID) + self.varAA
		else:
			print("varResID and varAA not specified")
			print("Using WT structure for simulation")
			self.variant = "wt"
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

def runInstance(parameters, node):
	trajDir = parameters.trajDIR
	scaffDir = parameters.scaffoldDIR
	if node.num == 0 and os.path.isdir(trajDir):
		shutil.rmtree(trajDir)
	if node.num == 0 and os.path.isdir(scaffDir):
		shutil.rmtree(scaffDir)
	if os.path.isdir(trajDir):
		os.rename(trajDir, trajDir + "_" + str(node.num - 1))
	if os.path.isdir(scaffDir):
		os.rename(scaffDir, scaffDir + "_" + str(node.num - 1))
	stashStruct(parameters, node)
	if not os.path.isdir(trajDir + "_" + str(node.num)):
                trajCommand = "python ../../snp2sim.py --config %s --mode varMDsim --newStruct %s --simID %d" %(parameters.config, node.scaff, node.num)
		trajCommand = "python %s/../../snp2sim.py --config %s --mode varMDsim --newStruct %s --simID %d" %(parameters.programDir, parameters.config, node.scaff, node.num)
		try:
			run_out = subprocess.check_output(trajCommand, shell = True)
			for line in run_out.decode().split("\n"):
				print(line)
		except subprocess.CalledProcessError as e:
			print("Error in varMDsim run %d: " %parameters.num)
			for line in e.output.decode().split("\n"):
				print(line)
			sys.exit(1)
	else:
		os.rename(trajDir + "_" + str(node.num), trajDir)
	if not os.path.isdir(scaffDir + "_" + str(node.num)):
                scaffCommand = "python ../../snp2sim.py --config %s --mode varScaffold --scaffID %d" %(parameters.config, node.num)
		scaffCommand = "python %s/../../snp2sim.py --config %s --mode varScaffold --scaffID %d" %(parameters.programDir, parameters.config, node.num)
		try:
			run_out = subprocess.check_output(scaffCommand, shell = True)
			for line in run_out.decode().split("\n"):
				print(line)
		except subprocess.CalledProcessError as e:
			print("Error in varScaffold run %d: \n" %parameters.num)
			for line in e.output.decode().split("\n"):
				print(line)
			sys.exit(1)
	else:
		os.rename(scaffDir + "_" + str(node.num), scaffDir)
	os.rename(trajDir, trajDir + "_" + str(node.num))
	os.rename(scaffDir, scaffDir + "_" + str(node.num))

def initialScaff(parameters):
	trajDir = parameters.trajDIR
	scaffDir = parameters.scaffoldDIR
	if os.path.isdir(trajDir):
		shutil.rmtree(trajDir)
	if os.path.isdir(scaffDir):
		shutil.rmtree(scaffDir)
	os.makedirs(trajDir)
	shutil.copy(parameters.initialTraj, trajDir)
	parameters.initialTraj = trajDir + parameters.initialTraj
	
	scaffCommand = "python %s/../../snp2sim.py --config %s --mode varScaffold --scaffID 0" %(parameters.programDir, parameters.config)
	try:
		run_out = subprocess.check_output(scaffCommand, shell = True)
		print(run_out.decode('ascii'))
	except subprocess.CalledProcessError as e:
		print("Error in varScaffold run %d: \n" %parameters.num)
		print(e.output.decode('utf-8'))
		sys.exit(1)

	shutil.move(trajDir[:-1], trajDir[:-1] + "_0")
	shutil.move(scaffDir[:-1], scaffDir[:-1] + "_0")
	curscaff = samplingTree(parameters, 0, parameters.initialTraj, None)
	curscaff.mark()
	scaffDir = os.path.join(parameters.resultsDIR, "scaffold_0")
	for file in sorted(os.listdir(scaffDir)):
		if file.endswith("scaffold.pdb"):
			parameters.num += 1
			scaff = scaffDir + "/" + file
			curscaff.addChild(samplingTree(parameters, parameters.num, scaff, curscaff))
	return curscaff
def stashStruct(parameters, node):
	structdir = parameters.structDIR
	shutil.move(structdir, structdir[:-1] + "_%d" % (node.num - 1))
	bindir = parameters.binDIR
	shutil.move(bindir, bindir[:-1] + "_%d" % (node.num - 1))
	configdir = parameters.configDIR
	shutil.move(configdir, configdir[:-1] + "_%d" % (node.num - 1))
	#if not os.path.isdir(os.path.join(structdir, "structs_%d" % (node.num - 1))):
		#os.mkdir(os.path.join(structdir, "structs_%d" % (node.num - 1)))
	#for file in sorted(os.listdir(structdir)):
		#if os.path.isfile(file):
			#os.system("mv %s %s" %(os.path.join(structdir, file), os.path.join(structdir, "structs_%d" % (node.num - 1), file)))
def checkStopCondition(parameters):
	#true if stop growing leaves
	if hasattr(parameters, "rmsdThresh"):
		oldthresh = parameters.rmsdThresh
	else:
		oldthresh = 0
	updateRMSDThresh(parameters)
	if parameters.rmsdThresh == oldthresh:
		return True
	else:
		return False
def updateRMSDThresh(parameters):
	scaffDir = parameters.scaffoldDIR
	for x in range(parameters.num):
		curRun = scaffDir + "_" + str(x)
		for file in sorted(os.listdir(curRun)):
			if file.endswith("scaffold.pdb"):
				scaffList.append(os.path.join(curRun, file))

	calcPairwiseRMSD(scaffList)
	with open(parameters.matrix, "r") as f:
		parameters.rmsd = [list(map(int, x.split(','))) for x in f.readlines()]
		parameters.rmsdThresh = max([max(x) for x in parameters.rmsd])

def calcPairwiseRMSD(parameters, pdbs):
	parameters.matrix = parameters.configDIR + "%s.%s.scaffpairwiseRMSD.csv" % (parameters.protein, parameters.variant)
	parameters.varPSF = parameters.structDIR + "%s.%s.UNSOLVATED.psf" \
							% (parameters.protein, parameters.variant)
	alignmentRes = parameters.alignmentResidues
	clusterRes = parameters.clusterResidues
	clustTCL = open(matrix ,"w")
	clustTCL.write("package require csv\n")
	for scaff in pdbs:
		clustTCL.write("mol addfile %s waitfor all\n" % scaff)

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

	clustTCL.write("set output [open %s w]\n" % parameters.matrix)
	clustTCL.write("set back [atomselect top \""+clusterRes+"\" frame 0]\n")
	clustTCL.write("set refclustRes [atomselect top \""+clusterRes+"\" frame 0]\n")

	#calculates RMSD and writes to a file
	clustTCL.write("for {set i 0} {$i < $nf} {incr i} {\n")
	clustTCL.write("set row {}\n")
	clustTCL.write("for {set j 0} {$j < $nf} {incr j} {\n")
	clustTCL.write("$refclustRes frame $i \n")
	clustTCL.write("$back frame $j \n")
	clustTCL.write("set M [measure rmsd $back $refclustRes] \n")
	clustTCL.write("lappend row $M}\n")
	clustTCL.write("puts $output [::csv::join $row]}\n")

	clustTCL.write("close $output\n")
	clustTCL.write("quit")
	clustTCL.close()

	makeTableCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
	os.system(makeTableCommand) 


def samplingNodes(parameters, tree):
	if not tree.child:
		if tree.extend:
			return tree
		else:
			return None
	bigList = []
	cur = [samplingNodes(parameters, x) for x in tree.child]
	for l in cur:
		if l:
			bigList.append(l)
	if tree.extend:
		bigList.append(tree)
		bigList.sort(key=lambda x: x.num)
	else:
		bigList.sort(key=lambda x: x.num)
	return bigList

def main():
	parameters = argParse()
	parameters.setDefaults()

	parameters.num = 0
	if parameters.newStruct:
		top = samplingTree(parameters, parameters.num, parameters.newStruct, None)
	elif parameters.initialTraj:
		top = initialScaff(parameters)
	else:
		print("newStruct or initialTraj required")
		sys.exit(1)

	curList = samplingNodes(parameters, top)
	while curList:
		curscaff = curList.pop(0)
		runInstance(parameters, curscaff)
		curscaff.mark()
		if not checkStopCondition():
			scaffDir = parameters.resultsDIR + "scaffold_%d" % curscaff.num
			for file in sorted(os.listdir(scaffDir)):
				if file.endswith("scaffold.pdb"):
					parameters.num += 1
					scaff = scaffDir + "/" + file
					curscaff.addChild(samplingTree(parameters, parameters.num, scaff, curscaff))
		curList = samplingNodes(parameters, top)

	print("Finished with sampling!")


main()







