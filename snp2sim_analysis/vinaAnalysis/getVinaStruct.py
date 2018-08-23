#!/usr/bin/env python

import os
import sys

vinaDir = sys.argv[1]
for vinaResultFile in os.listdir(vinaDir):
    if vinaResultFile.endswith(".pdbqt"):
        vinaResultFile = vinaDir + "/" + vinaResultFile
        vinaPDBFile = vinaResultFile.replace(".pdbqt",".pdb")
        vinaResults = open(vinaResultFile, "r")
        vinaResultsContent = vinaResults.readlines()
        vinaPDB = open(vinaPDBFile, "w")
        pose = 1
        for line in vinaResultsContent:
            if line.startswith("MODEL"):
                vinaPDB.write("MODEL %s\n" % pose)
                pose = pose + 1
            elif line.startswith("HETATM"):
                vinaPDB.write(line)
            elif line.startswith("ATOM"):
                vinaPDB.write(line)
            elif line.startswith("ENDMDL"):
                vinaPDB.write(line)

        vinaResults.close()
        vinaPDB.close()
