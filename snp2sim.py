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
                        help="name of protein system"
                        action="store",
                        type=str,
                        )
    parser.add_argument("--variant",
                        help="name of protein system"
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


parameters = _parseCommandLine()

##future work - define workflow from config

if parameters.mode == "varTraj":
    #check for mode parameters
    if parameters.simLength:
        print "Performing Variant %ins Simulation" % parameters.simLength
    else:
        print "simulation length not specified"
        exit

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


