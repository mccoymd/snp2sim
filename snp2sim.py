#!/usr/bin/env python

import os
import sys
import argparse
import copy

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

    cmdlineparameters, unknownparams = parser.parse_known_args()
    parameters = copy.copy(cmdlineparameters)

    return parameters

def genNAMDstructFiles(parameters):
    parameters.templatePDB = "./variantSimulations/%s/structures/%s.template.pdb" \
                             % (parameters.protein,parameters.protein)

    if not os.path.exists("./variantSimulations/%s/bin/" % parameters.protein):
        os.makedirs("./variantSimulations/%s/bin/" % parameters.protein)

    if os.path.isfile(parameters.templatePDB):
        print "generating solvated/ionized PDB and PSF from %s" \
            % parameters.templatePDB

            
        if os.path.isfile("./variantSimulations/%s/bin/%s.%s.genStructFiles.tcl" \
                          % (parameters.protein, parameters.protein, parameters.variant)):
            print "%s.%s.genStructFiles.tcl exists" \
                % (parameters.protein, parameters.variant)
        else:
            parameters = genStructTCL(parameters)

        #TODO generate solvated/ionized PDB and PSF using vmd




        parameters = solvStructTCL(parameters)
        sys.exit()
    else:
        print "no template exists"
        print "remove ./variantSimulations/%s directory" % parameters.protein
        print "use --new to specify clean PDB (only cannonical aa) "
        sys.exit()

    return parameters

def genStructTCL(parameters):
    parameters.basePDB = "./variantSimulations/%s/structures/%s.wt.UNSOLVATED.pdb" \
                         % (parameters.protein, parameters.protein)
    parameters.basePSF = "./variantSimulations/%s/structures/%s.wt.UNSOLVATED.psf" \
                         % (parameters.protein, parameters.protein)

    if not os.path.isfile("./variantSimulations/%s/bin/%s.wt.genStructFiles.tcl" \
                      % (parameters.protein, parameters.protein)):
        print "generating wt struct file %s.wt.genStructFiles.tcl" \
                      % (parameters.protein)


        structFileName = "./variantSimulations/%s/bin/%s.wt.genStructFiles.tcl" \
                         % (parameters.protein, parameters.protein)
        structFile = open(structFileName,"w+")
        
        structFile.write("package require psfgen\n")
        for tFile in parameters.simTopology:
            structFile.write("topology %s\n" % tFile)
        structFile.write("pdbalias residue HIS HSE\n")
        structFile.write("segment PROT {pdb %s}\n" % parameters.templatePDB)
        structFile.write("coordpdb %s PROT\n" % parameters.templatePDB)
        structFile.write("guesscoord\n")
        structFile.write("writepdb %s\n" % parameters.basePDB)
        structFile.write("writepsf %s\n" % parameters.basePSF)
        structFile.write("quit\n")
                        
    varStructFileName = "./variantSimulations/%s/bin/%s.%s.genStructFiles.tcl" \
                        % (parameters.protein, parameters.protein,parameters.variant)
                         
    if not os.path.isfile(varStructFileName):
        outPrefix = "./variantSimulations/%s/structures/%s.%s.UNSOLVATED" \
                    % (parameters.protein, parameters.protein, parameters.variant)
        varStructFile = open(varStructFileName,"w+")
        varStructFile.write("mutator -psf %s -pdb %s -o %s -ressegname PROT -resid %s -mut %s\n" \
                            % (parameters.basePSF, parameters.basePDB, \
                               outPrefix, parameters.varResID, parameters.varAA))
        varStructFile.write("quit\n");

    
    return(parameters)

def solvStructTCL(parameters):

    return(parameters)

def genNAMDconfig(parameters):
    #todo generate NAMD config file
    return(parameters)

parameters = _parseCommandLine()

### Hardcoded Parameters
parameters.simTopology = ("./simParameters/top_all36_prot.rtf",
                          "./simParameters/toppar_water_ions_namd.str")
parameters.simParameters = ("./simParameters/par_all36_prot.prm",
                            "./simParameters/toppar_water_ions_namd.str")



#todo: revamp to define workflow from config
if parameters.protein:
    print parameters.protein
else:
    print "no protein specified"
    sys.exit()

if parameters.mode == "varMDsim":
    print "Performing varMDsim"
    if parameters.newStruct:
        if os.path.isdir("./variantSimulations/%s" % parameters.protein):
            print "ERROR - %s is already created"
            sys.exit()
        os.makedirs("./variantSimulations/%s/structures" % parameters.protein)
        os.system("cp %s ./variantSimulations/%s/structures/%s.template.pdb" \
                  % (parameters.newStruct, parameters.protein, parameters.protein)) 

    if parameters.variant:
        #todo - check that variant format is correct "wt" or "x###x"

        ## parameters of resid and aa hardcoded in...
        
        #check if variant stucture files have already been generated
        parameters.variantPDB = "./variantSimulations/%s/structures/%s.%s.pdb" \
                                % (parameters.protein,parameters.protein,parameters.variant)
        parameters.variantPSF = "./variantSimulations/%s/structures/%s.%s.psf" \
                                % (parameters.protein,parameters.protein,parameters.variant)
        if os.path.isfile(parameters.variantPDB):
            if os.path.isfile(parameters.variantPSF):
                print "using %s %s as initial structure" \
                    % (parameters.protein, parameters.variant)

            #if no variant structure exists, use clean template to create
            else:
                print "%s does not exist" % parameters.variantPSF
                parameters = genNAMDstructFiles(parameters)
        else:
            print "%s does not exist" % parameters.variantPDB
            parameters = genNAMDstructFiles(parameters)
    else:
        print "no variant specified... using wildtype"
        parameters.variant = "wt"
        parameters.variantPDB = "./variantSimulations/%s/structures/%s.%s.pdb" \
                                % (parameters.protein,parameters.protein,parameters.variant)
        parameters.variantPSF = "./variantSimulations/%s/structures/%s.%s.psf" \
                                % (parameters.protein,parameters.protein,parameters.variant)
        if os.path.isfile(parameters.variantPDB):
            if os.path.isfile(parameters.variantPSF):
                print "using %s %s as initial structure" \
                    % (parameters.protein, parameters.variant)

            #if no variant structure exists, use clean template to create
            else:
                print "%s does not exist" % parameters.variantPSF
                parameters = genNAMDstructFiles(parameters)
        else:
            print "%s does not exist" % parameters.variantPDB
            parameters = genNAMDstructFiles(parameters)


    
    
    if parameters.simLength:
        print "Performing Variant %ins Simulation" % parameters.simLength
        parameters = genNAMDconfig(parameters)
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


