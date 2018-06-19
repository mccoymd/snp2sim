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
                        help="varTraj"
                        "mdScaffold"
                        "drugSearch"
                        "varAnalysis",
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
        print "no template exists "
        sys.exit()

    return parameters

def genNAMDconfig(parameters):
    #todo generate NAMD config file
    return(parameters)

parameters = _parseCommandLine()

#todo: revamp to define workflow from config
if parameters.mode == "varTraj":
    print "Performing varTraj"
    if parameters.protein:
        print "simulating %s" % parameters.protein
        if parameters.variant:
            #todo - check that variant format is correct "wt" or "x###x"
            #check if variant stucture files have already been generated
            parameters.variantPDB
            = "./variantSimulations/%s/structures/%s.%s.pdb" \
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
            #todo add wildtype simulation
            print "WARNING: auto wildtype simulation in development"
            sys.exit()
    else:
        print "no protein specified"
        sys.exit()
    

    
    if parameters.simLength:
        print "Performing Variant %ins Simulation" % parameters.simLength
        parameters = genNAMDconfig(parameters)
        #todo run NAMD simulation
    else:
        print "simulation length not specified"
        sys.exit()

    #check for required files

elif parameters.mode == "mdScaffold":
    #check for mode parameters
    #check for mode parameters
    print "Performing mdScaffold"
    
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


