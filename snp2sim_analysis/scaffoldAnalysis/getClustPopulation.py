#!/usr/bin/env python

import os
import sys

clustLogDir = sys.argv[1]

print "Variant\tClusterID\tPercent"
for clustLogFile in os.listdir(clustLogDir):
    if clustLogFile.endswith(".log"):
        clustStats = clustLogFile.split(".")
        clustLogFile = clustLogDir + "/" + clustLogFile
        clustLog = open(clustLogFile,"r")
        clustLogLines = clustLog.readlines()
        clusterMembership = clustLogLines[1]
        clusterMembership = clusterMembership.rstrip()
        clusterMembership = clusterMembership.replace("} {","_")
        clusterMembership = clusterMembership.replace("{","")
        clusterMembership = clusterMembership.replace("}","")
        clusterMembership = clusterMembership.split("_")
        totalStruct = 0
        for cluster in clusterMembership:
            clusterStruct = cluster.split(" ")
            totalStruct = totalStruct + len(clusterStruct)
        clusterNum = 1
        for cluster in clusterMembership:
            clusterStruct = cluster.split(" ")
            clustPercent = len(clusterStruct) / float(totalStruct)
            print "%s\t%i\t%.2f" % (clustStats[1], clusterNum, clustPercent)
            clusterNum = clusterNum + 1
