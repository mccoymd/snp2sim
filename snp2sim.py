#!/usr/bin/env python

import os
import sys
import argparse
import copy
import random

def _parseCommandLine():

    parser = argparse.ArgumentParser(
        version="0.1",
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
    
    parser.add_argument("--simLength",
                        help="varTraj simulation length in ns",
                        action="store",
                        type=int,
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
                        help="amino acid to change",
                        action="store",
                        type=str,
                        )

    cmdlineparameters, unknownparams = parser.parse_known_args()
    parameters = copy.copy(cmdlineparameters)

    return parameters


def genNAMDstructFiles(parameters):
    if not os.path.exists("./variantSimulations/%s/bin/" % parameters.protein):
        os.makedirs("./variantSimulations/%s/bin/" % parameters.protein)
    if os.path.isfile(parameters.templatePDB):
        print "generating solvated/ionized PDB and PSF from %s" \
            % parameters.templatePDB

        if os.path.isfile(parameters.varPDB) and os.path.isfile(parameters.varPSF):
            print "unsolvated PDB and PSF exist"
        else:
            genStructTCL(parameters)
            genStructCommand = "%s -e %s" % (parameters.vmdPath, parameters.varStructTCL)
            os.system(genStructCommand)

        genSolvTCL(parameters)
        genSolvCommand = "%s -e %s" % (parameters.vmdPath, parameters.solvTCL)
        os.system(genSolvCommand)
    else:
        print "no template exists"
        print "use --new to specify clean PDB (only cannonical aa) "
        sys.exit()

    return

def genStructTCL(parameters):
    print "generating struct files"
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
                   "H":"HIS","R":"ARG","N":"ASN","D":"ASP","T":"THR"}
        structFile.write("mutator -psf %s -pdb %s -o %s -ressegname PROT -resid %s -mut %s\n" \
                         % (parameters.wtPSF, parameters.wtPDB, \
                            parameters.varPrefix, parameters.varResID, longAA.get(parameters.varAA)))
    
    structFile.write("quit\n")

    return

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
    return

def genNAMDconfig(parameters):
    if not os.path.exists("./variantSimulations/%s/config/" % parameters.protein):
        os.makedirs("./variantSimulations/%s/config/" % parameters.protein)
    configFile = open(parameters.NAMDconfig,"w+")
    configFile.write("structure          %s\n" % parameters.simPSF)
    configFile.write("coordinates        %s\n" % parameters.simPDB)
    configFile.write("set outputname     %s\n" % parameters.NAMDout)
    configFile.write("paraTypeCharm on\n")
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
    configFile.write("langevinHydrogen    off    ;# don't couple langevin bath to hydrogens\n")
    
    
    
    return

parameters = _parseCommandLine()

### Hardcoded Parameters

if parameters.protein:
    print parameters.protein
    parameters.vmdPath = "vmd"
    parameters.namdPath = "namd"
    if parameters.mode == "varMDsim":
        if not parameters.simID:
            parameters.simID = random.randint(1,1000)
            
        parameters.simTopology = ("./simParameters/top_all36_prot.rtf",
                                  "./simParameters/toppar_water_ions_namd.str")
        parameters.simParameters = ("./simParameters/par_all36_prot.prm",
                                    "./simParameters/toppar_water_ions_namd.str")
        if parameters.varResID and parameters.varAA:
            parameters.variant = parameters.varResID + parameters.varAA
        else:
            print "varResID and varAA not specified"
            print"Using WT structure for simulation"
            parameters.variant = "wt"
            
        parameters.simPDB = "./variantSimulations/%s/structures/%s.%s.pdb" \
                            % (parameters.protein,parameters.protein,parameters.variant)
        parameters.simPSF = "./variantSimulations/%s/structures/%s.%s.psf" \
                            % (parameters.protein,parameters.protein,parameters.variant)
        parameters.templatePDB = "./variantSimulations/%s/structures/%s.template.pdb" \
                                 % (parameters.protein,parameters.protein)
        parameters.wtPDB = "./variantSimulations/%s/structures/%s.wt.UNSOLVATED.pdb" \
                             % (parameters.protein, parameters.protein)
        parameters.wtPSF = "./variantSimulations/%s/structures/%s.wt.UNSOLVATED.psf" \
                             % (parameters.protein, parameters.protein)
        parameters.varPrefix = "./variantSimulations/%s/structures/%s.%s.UNSOLVATED" \
                               % (parameters.protein, parameters.protein, parameters.variant)
        parameters.varPDB = "./variantSimulations/%s/structures/%s.%s.UNSOLVATED.pdb" \
                             % (parameters.protein, parameters.protein, parameters.variant)
        parameters.varPSF = "./variantSimulations/%s/structures/%s.%s.UNSOLVATED.psf" \
                             % (parameters.protein, parameters.protein, parameters.variant)
        parameters.wtStructTCL = "./variantSimulations/%s/bin/%s.wt.genStructFiles.tcl" \
                                 % (parameters.protein, parameters.protein)
        parameters.varStructTCL = "./variantSimulations/%s/bin/%s.%s.genStructFiles.tcl" \
                                  % (parameters.protein, parameters.protein, parameters.variant)
        parameters.solvTCL = "./variantSimulations/%s/bin/%s.%s.genSolvStruct.tcl" \
                             % (parameters.protein, parameters.protein, parameters.variant)
        parameters.solvBoundary = "./variantSimulations/%s/bin/%s.%s.solvBoundary.txt" \
                                  % (parameters.protein, parameters.protein, parameters.variant)
        parameters.structPrefix = "./variantSimulations/%s/structures/%s.%s" \
                                  % (parameters.protein, parameters.protein, parameters.variant)
        parameters.NAMDconfig = "./variantSimulations/%s/config/%s.%s.%s.NAMD" \
                                  % (parameters.protein, parameters.protein,
                                     parameters.variant, parameters.simID)
        parameters.NAMDout = "./variantSimulations/%s/results/%s/trajectory/%s.%s.%s" \
                             % (parameters.protein, parameters.variant,
                                parameters.protein, parameters.variant, parameters.simID)
                                
        
else:
    print "no protein specified"
    sys.exit()

#todo: revamp to define workflow from config

if parameters.mode == "varMDsim":
    print "Performing varMDsim"
    if parameters.newStruct:
        if os.path.isdir("./variantSimulations/%s" % parameters.protein):
            print "ERROR - %s is already created... resubmit with new protien name"
            sys.exit()
        os.makedirs("./variantSimulations/%s/structures" % parameters.protein)
        os.system("cp %s ./variantSimulations/%s/structures/%s.template.pdb" \
                  % (parameters.newStruct, parameters.protein, parameters.protein)) 

    if os.path.isfile(parameters.simPDB):
        if os.path.isfile(parameters.simPSF):
            print "using %s %s as initial structure" \
                % (parameters.protein, parameters.variant)

            #if no variant structure exists, use clean template to create
        else:
            print "%s does not exist" % parameters.simPSF
            print "generating new pdb and psf for simulation"
            genNAMDstructFiles(parameters)
    else:
        print "%s does not exist" % parameters.simPDB
        print "generating new pdb and psf for simulation"
        genNAMDstructFiles(parameters)

    if parameters.simLength:
        print "Performing Variant %ins Simulation" % parameters.simLength
        genNAMDconfig(parameters)
        #todo run NAMD simulation
    else:
        print "simulation length not specified"        
        sys.exit()

#######

elif parameters.mode == "varScaffold":
    #check for mode parameters
    print "Performing varScaffold"

    #check if variant specified, otherwise perform clustering for all variants
    #if new clustParamID, create config
    #check if clustParamID config exists
    #create TCL
    #create template structures
    


    
elif parameters.mode == "drugSearch":
    #check for mode parameters
    #check for mode parameters
    print "Performing drugSearch"
    
elif parameters.mode == "varAnalysis":
    #check for mode parameters
    #check for mode parameters
    print "Performing varAnalysis"
    
else:
    print "invalid mode selected"
    
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


