#!/usr/bin/env python

import os
import sys
import argparse
import copy
import random
import multiprocessing

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
    parser.add_argument("--simID",
                        help="amino acid to change",
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
    parser.add_argument("--clustPDBtraj",
                        help="cluster results from multiple varScafold singleRun",
                        action="store_true",
                        )
    parser.add_argument("--loadPDBtraj", nargs='+',
                        help="pdb trajectory files to import, list one after another",
                        action="store",
                        )

    cmdlineparameters, unknownparams = parser.parse_known_args()
    parameters = copy.copy(cmdlineparameters)

    return parameters


def genNAMDstructFiles(parameters):
    if not os.path.exists("%s/variantSimulations/%s/bin/" % (parameters.runDIR,parameters.protein)):
        os.makedirs("%s/variantSimulations/%s/bin/" % (parameters.runDIR,parameters.protein))
    if os.path.isfile(parameters.templatePDB):
        print "generating solvated/ionized PDB and PSF from %s" \
            % parameters.templatePDB

        if os.path.isfile(parameters.varPDB) and os.path.isfile(parameters.varPSF):
            print "unsolvated PDB and PSF exist"
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
        print "boundary file does not exist"
        print "remove solvated/ionized PDB and PSF and resubmit"
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
                     

    
    return

def genSingleRunTCL(parameters):
    print "generating %s" % parameters.singleRunTCL
    pdbTCL = open(parameters.singleRunTCL,"w+")
    pdbTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
    variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
    print "using DCD files in %s" % variantDIR
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
    
    return

def genClusterTCL(parameters):
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
    print scaffLOG
    if os.path.isfile(scaffLOG):
        print "%s exists. Remove to recalculate Scaffolds" % scaffLOG
        return
        
    print "generating %s" % parameters.scaffoldTCL
    clustTCL = open(parameters.scaffoldTCL,"w+")
    clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
    variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
    print "using DCD files in %s" % variantDIR
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
    
    return

def genPDBclustTCL(parameters):
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
    print scaffLOG
    if os.path.isfile(scaffLOG):
        print "%s exists. Remove to recalculate Scaffolds" % scaffLOG
        return

    print "generating %s" % parameters.scaffoldTCL

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
#    clustTCL.write("mol new %s waitfor all\n" % parameters.simPSF)
#    variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
#    print "using DCD files in %s" % variantDIR
#    for tFile in os.listdir(variantDIR):
#        if tFile.endswith(".dcd"):
#            dcdFile = variantDIR + tFile
#            clustTCL.write("mol addfile %s waitfor all\n" % dcdFile)

#    clustTCL.write("package require pbctools\n")
#    clustTCL.write("pbc wrap -center com -centersel \"protein\" -compound residue -all\n")
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
    
    return

def sortPDBclusters(parameters):
    allPDBoutput = parameters.scaffBASE + ".all.pdb"
    scaffLOG = parameters.scaffBASE + ".log"    

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

    clustLogfile = open(scaffLOG, "r")
    clustLogLines = clustLogfile.readlines()
    clusterMembership = clustLogLines[1]
    clusterMembership = clusterMembership.rstrip()
    clusterMembership = clusterMembership.replace("} {","_")
    clusterMembership = clusterMembership.replace("{","")
    clusterMembership = clusterMembership.replace("}","")
    clusterMembership = clusterMembership.split("_")
    del clusterMembership[-1]
    
    scaffNum = 1;
    for indCluster in clusterMembership:
        indCluster = indCluster.split(" ")
        if len(indCluster) > 0.1*len(indPDB):
            scaffFileName = parameters.scaffBASE + ".cluster" + str(scaffNum) + ".pdb"
            scaffFile = open(scaffFileName,"w+")
            scaffFile.write(pdbHeader)
            for structID in indCluster: 
                scaffFile.write(indPDB[int(structID)])
            scaffNum += 1 

    return


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
                genScaff.write("set domain [atomselect top %s]\n" % alignmentRes)
                genScaff.write("set avePos [measure avpos $domain]\n")
                genScaff.write("$domain set {x y z} $avePos\n")
                genScaff.write("$domain writepdb %s\n" % pdbScaffFile)
                
    genScaff.write("quit\n")
    return

#### start of program move to main()

parameters = _parseCommandLine()

### Hardcoded Parameters

#parameters.runDIR = os.getcwd()
parameters.runDIR = os.path.abspath(__file__)
parameters.runDIR = os.path.dirname(parameters.runDIR)
print(parameters.runDIR)

if parameters.protein:
    print parameters.protein
    if not parameters.VMDpath:
        parameters.VMDpath = "vmd"
    if not parameters.NAMDpath:
        parameters.NAMDpath = "namd2"
        
    if parameters.mode == "varMDsim":
        if not parameters.simID:
            parameters.simID = str(random.randint(1,1000))

        if parameters.varResID and parameters.varAA:
            parameters.variant = parameters.varResID + parameters.varAA
        else:
            print "varResID and varAA not specified"
            print"Using WT structure for simulation"
            parameters.variant = "wt"

    if parameters.mode == "varScaffold":
        if not parameters.scaffID:
            print "no scaff ID - exiting"
            sys.exit()


        
    parameters.simTopology = ("%s/simParameters/top_all36_prot.rtf" % parameters.runDIR,
                              "%s/simParameters/toppar_water_ions_namd.str" % parameters.runDIR)
    parameters.simParameters = ("%s/simParameters/par_all36_prot.prm" % parameters.runDIR,
                                "%s/simParameters/toppar_water_ions_namd.str" %parameters.runDIR)

    if not parameters.simProc:
        parameters.simProc = multiprocessing.cpu_count()    
    
    parameters.simPDB = "%s/variantSimulations/%s/structures/%s.%s.pdb" \
                        % (parameters.runDIR,parameters.protein,
                           parameters.protein,parameters.variant)
    parameters.simPSF = "%s/variantSimulations/%s/structures/%s.%s.psf" \
                        % (parameters.runDIR,parameters.protein,
                           parameters.protein,parameters.variant)
    parameters.templatePDB = "%s/variantSimulations/%s/structures/%s.template.pdb" \
                             % (parameters.runDIR,parameters.protein,
                                parameters.protein)
    parameters.wtPDB = "%s/variantSimulations/%s/structures/%s.wt.UNSOLVATED.pdb" \
                       % (parameters.runDIR,parameters.protein,
                          parameters.protein)
    parameters.wtPSF = "%s/variantSimulations/%s/structures/%s.wt.UNSOLVATED.psf" \
                       % (parameters.runDIR,parameters.protein,
                          parameters.protein)
    parameters.varPrefix = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED" \
                           % (parameters.runDIR,parameters.protein,
                              parameters.protein, parameters.variant)
    parameters.varPDB = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED.pdb" \
                        % (parameters.runDIR, parameters.protein,
                           parameters.protein, parameters.variant)
    parameters.varPSF = "%s/variantSimulations/%s/structures/%s.%s.UNSOLVATED.psf" \
                        % (parameters.runDIR,parameters.protein,
                           parameters.protein, parameters.variant)
    parameters.wtStructTCL = "%s/variantSimulations/%s/bin/%s.wt.genStructFiles.tcl" \
                             % (parameters.runDIR,parameters.protein, parameters.protein)
    parameters.varStructTCL = "%s/variantSimulations/%s/bin/%s.%s.genStructFiles.tcl" \
                              % (parameters.runDIR,parameters.protein,
                                 parameters.protein, parameters.variant)
    parameters.singleRunTCL = "%s/variantSimulations/%s/bin/%s.%s.genSingleRunPDB.tcl" \
                              % (parameters.runDIR,parameters.protein,
                                 parameters.protein, parameters.variant)
    parameters.solvTCL = "%s/variantSimulations/%s/bin/%s.%s.genSolvStruct.tcl" \
                         % (parameters.runDIR,parameters.protein,
                            parameters.protein, parameters.variant)
    parameters.solvBoundary = "%s/variantSimulations/%s/config/%s.%s.solvBoundary.txt" \
                              % (parameters.runDIR,parameters.protein,
                                 parameters.protein, parameters.variant)
    parameters.structPrefix = "%s/variantSimulations/%s/structures/%s.%s" \
                              % (parameters.runDIR,parameters.protein,
                                 parameters.protein, parameters.variant)
    parameters.NAMDconfig = "%s/variantSimulations/%s/config/%s.%s.%s.NAMD" \
                            % (parameters.runDIR, parameters.protein, parameters.protein,
                               parameters.variant, parameters.simID)
    parameters.NAMDout = "%s/variantSimulations/%s/results/%s/trajectory/%s.%s.%s" \
                         % (parameters.runDIR, parameters.protein, parameters.variant,
                            parameters.protein, parameters.variant, parameters.simID)
    parameters.scaffParams = "%s/variantSimulations/%s/config/%s.scaff" \
                         % (parameters.runDIR, parameters.protein, parameters.scaffID)
    parameters.resultsDIR =  "%s/variantSimulations/%s/results" % \
                             (parameters.runDIR, parameters.protein)
    parameters.trajDIR = "%s/variantSimulations/%s/results/%s/trajectory/" \
                         % (parameters.runDIR, parameters.protein, parameters.variant)
                            

    if parameters.loadPDBtraj:
        variantDIR = parameters.resultsDIR + "/" + parameters.variant + "/trajectory/"
        if not os.path.exists(variantDIR):
            os.makedirs(variantDIR)

        for pdbFile in parameters.loadPDBtraj:
            pdbNewLoc = parameters.trajDIR + os.path.basename(pdbFile)
            print pdbNewLoc
            os.system("mv %s %s" % (pdbFile,pdbNewLoc))
    

else:
    print "no protein specified"
    sys.exit()

#todo: revamp to define workflow from config

if parameters.mode == "varMDsim":
    print "Performing varMDsim"
    if parameters.newStruct:
        if os.path.isdir("%s/variantSimulations/%s" % (parameters.runDIR,parameters.protein)):
            print "ERROR - %s is already created... resubmit with new protien name"
            sys.exit()
        os.makedirs("%s/variantSimulations/%s/structures" % (parameters.runDIR,parameters.protein))
        os.system("cp %s %s/variantSimulations/%s/structures/%s.template.pdb" \
                  % (parameters.newStruct,parameters.runDIR,
                     parameters.protein, parameters.protein)) 

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
        print "Performing Variant %.3f ns Simulation" % parameters.simLength
        genNAMDconfig(parameters)

        runNAMDcommand = "%s +p%i %s > %s.log" % \
                         (parameters.NAMDpath, parameters.simProc,
                          parameters.NAMDconfig, parameters.NAMDout)
        if not os.path.isdir("%s/variantSimulations/%s/results/%s/trajectory/" \
                             % (parameters.runDIR, parameters.protein, parameters.variant)):
            os.makedirs("%s/variantSimulations/%s/results/%s/trajectory/" \
                             % (parameters.runDIR, parameters.protein, parameters.variant))


        print "running NAMD with %i processors" % parameters.simProc
        os.system(runNAMDcommand)
        
        if parameters.singleRun:
            genSingleRunTCL(parameters)
            genPDBcommand = "%s -e %s" % (parameters.VMDpath, parameters.singleRunTCL)
            os.system(genPDBcommand)
            
        
    else:
        print "simulation length not specified"        
        sys.exit()

#######

elif parameters.mode == "varScaffold":
    print "Performing varScaffold"

    
    
    if parameters.newScaff:
        if not os.path.isdir("%s/variantSimulations/%s/config" % \
                             (parameters.runDIR,parameters.protein)):
            os.makedirs("%s/variantSimulations/%s/config" % \
                        (parameters.runDIR,parameters.protein))
        if not os.path.isfile(parameters.scaffParams):
            os.system("cp %s %s/variantSimulations/%s/config/%s.scaff" \
                      % (parameters.newScaff,parameters.runDIR,
                         parameters.protein, parameters.scaffID))

        else:
            print "Scaff config for %s exists." % (parameters.scaffParams)
            print "Select new scaffID or remove existing scaff config"
            sys.exit()

    if os.path.isfile(parameters.scaffParams):
        #check if variant specified, otherwise perform clustering for all variants
        if parameters.variant:
            print "generating Scaffold for %s ONLY" % (parameters.variant)
            variantList = (parameters.variant,)
        else:
            print "generating all Variant Scaffolds"
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


            print "generating Scaffold for %s" % (parameters.variant)
            #will not generate new TCL if previous scaffID log exists
            scaffLOG = parameters.scaffBASE + ".log"
            print scaffLOG
            if not os.path.isfile(scaffLOG):
                if not parameters.clustPDBtraj:
                    genClusterTCL(parameters)
                    vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
                    os.system(vmdClustCommand)
                else:
                    parameters.scaffoldTCL =  "%s/variantSimulations/%s/bin/%s.%s.%s.genPDBtrajScaffold.tcl" % \
                                              (parameters.runDIR, parameters.protein,
                                               parameters.protein, parameters.variant, parameters.scaffID)
                    genPDBclustTCL(parameters)
                    vmdClustCommand = "%s -e %s" % (parameters.VMDpath, parameters.scaffoldTCL)
                    os.system(vmdClustCommand)                    

                sortPDBclusters(parameters)
                genScaffoldTCL(parameters)
                vmdScaffCommand = "%s -e %s" % (parameters.VMDpath, parameters.clusterTCL)
                os.system(vmdScaffCommand)
#                else:
#                    print "only calculating PDB trajectory from DCD"
#                    trajID = str(random.randint(1,1000))
#                    oldTrajName = parameters.scaffBASE + ".all.pdb"
#                    newTrajName = parameters.scaffBASE + ".run" + trajID + ".pdb"
#                    os.system("mv %s %s" % (oldTrajName, newTrajName))
#                    sys.exit()

                
            else:
                print "%s exists. Remove to regenerate scaffold stats" % scaffLOG

    else:
        print "Scaffold Parameters %s does not exist" % (parameters.scaffParams)
        print "Exiting"


    
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


