# snp2sim
Molecular Simulation of Somatic Variation

Usage:

Command Line Options:
* mode
varTraj - build variant specific trajectory
mdScaffold - cluster trajectorys into representative scaffold
drugSearch - bind drug library to variant scaffolds
varAnalysis - analyze drugSearch results for multiple variants

* protein
name used to refer to group of variants

*variant
variant to simulate

Mode specific options:
varTraj
  --simulation time

mdScaffold
  --alignmentResID
  --clusterResID
  --clusterThresh

drugSearch
  --searchSpace
  --library

varAnalysis

--------
additional tools (in development)
  buildSearchSpace - defines protein orientation and search parameters
  buildDrugLibrary - creates new library from set of small molecules pdb files

