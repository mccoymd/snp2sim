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
    if os.path.isfile(parameters.templatePDB):
        print "generating solvated/ionized PDB and PSF from %s" \
            % parameters.templatePDB
        #TODO generate tcl to create variant structure pdb and psf
        #TODO generate solvated/ionized PDB and PSF using vmd
    else:
        print "no template exists"
        print "remove ./variantSimulations/%s directory" % parameters.protein
        print "use --new to specify clean PDB (only cannonical aa) "
        sys.exit()

    return parameters

def genNAMDconfig(parameters):
    #todo generate NAMD config file
    return(parameters)

parameters = _parseCommandLine()

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


