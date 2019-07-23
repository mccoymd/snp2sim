# SNP2SIM
Molecular Simulation of Protein Structure Variants
Questions: Matthew McCoy - mdm299 "at" georgetown.edu

## Installation Instructions
The SNP2SIM workflow is a Python script which controls the 
execution of molecular simulations to generate
variant specific structural scaffolding for small molecule 
docking. The script runs on a linux operating system and 
controls the execution of third party software detailed below.
  
### Nanoscale Molecular Dynamics (NAMD)
NAMD is used to execute the molecular dynamics simulations,
and can be downloaded from 
https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD
The binaries should be added to the users PATH environment variable as "namd2".

Alternatively, the user can specify the location of the executable 
using the `--NAMDpath "path to NAMD binary"` command line option.

### Visual Molecular Dynamics (VMD)
VMD is used to maniputulate and analyze protein structure files 
and can be downloaded from 
https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
The binaries should be added to the users PATH environment variable as "vmd".

Alternatively, the user can specify the location of the executable 
using the `--VMDpath "path to VMD binary"` command line option.

### AutoDock Vina and AutoDockTools
AutoDock Vina is used for the small molecule binding simulations
and can be downloaded from 
http://vina.scripps.edu/download.html
The AutoDock Vina executable should be added to the users PATH environment variable 
as "vina".

Alternatively, the user can specify the location of the executable 
using the `--VINApath "path to AutoDock Vina binary"` command line option.

From the included AutoDock Tools, the python wrapper script "pythonsh" 
should be added to to the users PATH environment 

AutoDockTools scripts prepare_receptor4.py and prepare_flexreceptor4.py 
should be installed to directory path 
`/opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/`

Alternatively, the user can specify the location of the pythonsh executable  
can be specified using the `--PYTHONSHpath "path to pythonsh"` command line option
and the lovation of the Ulititys24 folder included with the AutoDock Tools distribution
can be specified with the `--ADTpath "path to Utilities24 directory` command line option.

### R dependencies
The varScaffold module and varAnalysis modules use R scripts to generate visualizations and cluster MD trajectories.
Along with the base installation of R, the following packages are necessary dependencies for the modules:

varScaffold dependencies: data.table, ggplot2, plotly, htmlwidgets, fpc, markovchain

varAnalysis dependencies: shiny, ggplot2, highcharter, viridis, plotly, webshot

## Example Usage
Usage:python snp2sim.py ""options""

The workflow is configured to store intermediate files and
results in a predefined directory structure. If the required
trajectory/scaffold files are not present in the snp2sim directory,
they must be specified through the command line.
  
The files used to run a case study using PD-L1 are provided
in the "example" directory.

### Generating Structural Trajectories using varMDsim

#### Input Files
The input file is a single chain of a PDB structure.

#### Command line options

WT simulation
note: If the name used as the --protein already exists in the output directory, 
the --newStruct option is not necessary and the workflow will utilize the
configureation files that have been previously generated. 

`python snp2sim.py --mode varMDsim --protein PDL1 --varResID 115 --varAA T  --newStruct example/PDL1.Vtype.pdb --simLength 0.1`


#### Config file options

General options:

```
# General options
#################
 # varMDsim, varScaffold, drugSearch, or varAnalysis
mode: varMDsim

 #Name of protein system
protein: PDL1

 #Variant as "wt" or "x###x" - (OPTIONAL for varAnalysis)
variant: 

#varResID + varAA override variant field. If one is filled, the other must be filled, or variant field used
 #variant residue ID from PDB template
varResID: 115

 #amino acid to change
varAA: T

 #if Run is on CGC platform, move pdb and log to snp2sim root for processing
cgcRun:

runDIR: /Path/To/SNP2SIM/Working/Dir/

#Overwrite previous run with same protein name.
clean: 


```

varMDsim specific options:

```
# Mode specific options

# varMDsim options
##################

 #path to cleaned PDB file (protein structure w/ cannonical aa)
newStruct: example/PDL1.Vtype.pdb

 #varTraj simulation length in ns
simLength: 0.1

 #amino acid to change
simID:

 #number of processors to run simulation
simProc:

 #output summary PDB trajectory only
singleRun:

 #Only generate initial structures, without MD simulation
genStructures:
```


#### Output files

PROTEIN refers to the value of the --protein command line option.
The module will output tcl scripts for the generation of the solvated 
protein structure in the "variantSimulations/PROTEIN/bin/" directory. 

NAMD configs are output to the "variantSimulations/PROTEIN/config/" directory.

NAMD input structure files (solvated psf and pdb files)
are output to the "variantSimulations/PROTEIN/structures/" directory.

VARIANT refers to the concatenated values of --varResID and --varAA

NAMD trajectory results are found in the 
"variantSimulations/PROTEIN/VARIANT/trajectory" directory.

### Generating Variant Scaffolds using varScaffold

#### Input Files
NAMD output, if it exists in the file structure defined by the 
--protein --varResID and --varAA command line options.

The input can also be a list of PDB formatted trajectory files using 
the --clustPDBtraj command line option.

#### Command line options

`python snp2sim.py --mode varScaffold --protein PDL1 --varResID 115 --varAA T --scaffID bindingRes `


#### Config file options

General options:

```
# General options
#################
 # varMDsim, varScaffold, drugSearch, or varAnalysis
mode: varScaffold

 #Name of protein system
protein: PDL1

 #Variant as "wt" or "x###x" - (OPTIONAL for varAnalysis)
variant: 

#varResID + varAA override variant field. If one is filled, the other must be filled, or variant field used
 #variant residue ID from PDB template
varResID: 115

 #amino acid to change
varAA: T

 #if Run is on CGC platform, move pdb and log to snp2sim root for processing
cgcRun:

runDIR: /Path/To/SNP2SIM/Working/Dir/

#Overwrite previous run with same protein name.
clean: 


```

varScaffold specific options:

```
# varScaffold
#############

 #unique identifier for run
scaffID: pdl_scaff_1

 #Use PDB trajectories (possibly from singleCGC runs)
clustPDBtraj:

 #pdb trajectory files to import, list one after another
loadPDBtraj:

#Cluster Parameters - follow VMD atomselection format
 #Residues to superimpose trajectory 
alignmentResidues: backbone and resid 19 to 131
 #Residues to consider when clustering trajectory
clusterResidues: backbone and resid 19 20 54 56 66 68 115 116 117 121 122 123 124 125

 #use a table of structural features to perform clustering, otherwise use a pairwise distance matrix
featureTableMethod: true
 #if mds method used, path to pairwise RMSD table, otherwise automatically generated
PairwiseRMSD: 

 #create a script that colors the cluster residues in VMD by cluster
colorTrajectory: true

```


#### Output files

The output files include a representative PDB for each cluster, 
as well as a log file with the cluster assignments for all 
trajectory structures are stored in the 
"variantSimulations/PROTEIN/VARIANT/scaffold" directory.

### Generating Small Molecule Docking Results using drugSearch

#### Input Files
The drugSearch module takes a list of PDB formatted structures, 
as well as a reference structure used to define the search space 
location and dimensions (provided in the config). The user can also specify 
residues to use 

#### Command line options

`python snp2sim.py --mode drugSearch --protein PDL1 --varResID 115 --varAA T 
--bindingTemplate example/PDL1.Vtype --newBindingConfig example/pdl1SearchSpace.txt  
--drugLibrary pdl1-SMI 
--inputScaff example/exampleResults/PDL1.115T.bindingRes.cl1.scaffold.pdb example/exampleResults/PDL1.115T.bindingRes.cl2.scaffold.pdb`

#### Config file options

General options:

```
# General options
#################
 # varMDsim, varScaffold, drugSearch, or varAnalysis
mode: drugSearch

 #Name of protein system
protein: PDL1

 #Variant as "wt" or "x###x" - (OPTIONAL for varAnalysis)
variant: 

#varResID + varAA override variant field. If one is filled, the other must be filled, or variant field used
 #variant residue ID from PDB template
varResID: 115

 #amino acid to change
varAA: T

 #if Run is on CGC platform, move pdb and log to snp2sim root for processing
cgcRun:

runDIR: /Path/To/SNP2SIM/Working/Dir/

#Overwrite previous run with same protein name.
clean: 


```

drugBinding specific options:

```
# drugSearch
############

 #path to PDB file used to create search space to align scaffold
bindingTemplate:

 #exhaustiveness parameter of autodock vina
vinaExh: 50

 #path to new binding config file
newBindingConfig:

 #Automatically determine the search box based on the 
autoSearchSpace: true

 #list of residues that comprise the search space if autoSearchSpace is true, e.i 22 21 19 103 113 125
searchResidues: 54 117 115 56 123 113 66 68 58 121 63 122 76 73

 #list of residue numbers in flexible binding pocket, e.i 22 21 19 103 113 125
flexBinding: 54 117 115 56 123 113 66 68 58 121 63 122 76 73

 #name of snp2sim drug library
drugLibrary: pdl1-SMI 

 #path to single drug PDBQT
singleDrug:

 #pdb scaffold files to import, list one after another
inputScaff:
 - example/exampleResults/PDL1.115T.bindingRes.cl1.scaffold.pdb 
 - example/exampleResults/PDL1.115T.bindingRes.cl2.scaffold.pdb

 #only bind single variant scaffolds
bindSingleVar:

 #path to dir with ligand PDBs, which will be converted to drug library
ligandPDB:

 #number of times to run the docking simulation, to get an uncertainty measurement of the binding energy
numTrials: 10
```

#### Output files
The log files containing the binding affities for the top
ligand poses, as well as the PDBQT files that contain the 
coordinates of the ligand and flexable residues are depositied in the 
"variantSimulations/PROTEIN/results/VARIANT/drugBinding" directory.

Additionally, the PDBQT files used for the AutoDock Vina
simulations can be found in the 
"variantSimulations/PROTEIN/results/VARIANT/scaffold" directory.
  

### Analyzing the results of the workflow

#### Input Files

Results of drugBinding are searched in the drugBinding results directory for each variant, along with scaffold results.

#### Config file options

General options:

```
# General options
#################
 # varMDsim, varScaffold, drugSearch, or varAnalysis
mode: varAnalysis

 #Name of protein system
protein: PDL1

runDIR: /Path/To/SNP2SIM/Working/Dir/

#Overwrite previous run with same protein name.
clean: 
```

varAnalysis specific options:

```
# varAnalysis
############

 #Variants to be used in variants (must have a results/variant/drugBinding directory with .pdbqt results)
 #Format ###AA (### = residue number and AA = variant amino acid) or wt
 #must include wt
analysisVariants: 
  - 53P
  - 68L
  - 86W
  - 94M
  - 95R
  - 97V
  - 115T
  - wt

 #Analyze results from specific drug libraries
analysisDrugLibrary: pdl1-LMW

 #Analyze results from a specific ligand
analysisDrug:

 #Enable if using scaffolding results from old versions of SNP2SIM
legacyScaff:

 #Enable if using drugBinding results from custom PDB proteins, and not varScaffold results
customScaff: 

```
#### Output files

The output files can be found in "variantSimulations/PROTEIN/results/analysis/list_of_variants_and_drugs". 
The results include a summary file with the drug binding results in tsv format. 
Select plots of the results are included in the Figures folder of the analysis directory, along with a script `interactive_visualize.sh` which starts an RShiny app with the drugBinding results loaded to viusalize interactively.


## Command Line and Config Options:
### General:
--mode "string"
  select which snp2sim module to run
  appropriate values: varMDsim, varScaffold, drugSearch

--protein "string"
  user specified name for simulation
  used in naming output directories and file

### Optional:
--varResID "integer"
  Position in input PDB file to mutate to varAA

--varAA "character"
  Single letter code for one of 20 canonical amino acids to mutate varResID

--NAMDpath "string"
  path to NAMD executatble

--VMDpath "string"
  path to VMD executable

--VINApath "string"
  path to AutoDock Vina execuatable

--PYTHONSHpath "string"
  path to AutoDockTools script PYTHONSH

--ATDpath "string"
  path th AutoDockTools Utilitys directory

--cgcRun
  If flag is set, output files will be copied to
  working directory (for output on CGC platform app)

### Config template

```
# General options
#################
 # varMDsim, varScaffold, drugSearch, or varAnalysis
mode: 

 #Name of protein system
protein:

 #Variant as "wt" or "x###x" - (OPTIONAL for varAnalysis)
variant:

#varResID + varAA override variant field. If one is filled, the other must be filled, or variant field used
 #variant residue ID from PDB template
varResID: 

 #amino acid to change
varAA: 

 #if Run is on CGC platform, move pdb and log to snp2sim root for processing
cgcRun:

runDIR: 

#Overwrite previous run with same protein name.
clean: 


#Paths to executables. Ignore if using Docker.

 #path to VMD executable
VMDpath: 

 #path to AutoDock Vina executable
VINApath:

 #path to NAMD executable
NAMDpath:

 #path to autodock tools python executable
PYTHONSHpath:

 #path to autodock utilities directory
ADTpath:
```
### Options specific to snp2sim modules
### varMDsim:
  Usage Notes:
  If a variant is now specified using --varResID and --varAA,
  the simulation will run on the unmutated structure

  Module Options
    --newStruct "file path"
      PDB file containing only the protein structure

    --simLength "real number"
      length of the NAMD simulation in nanoseconds

    --simID "string" (optional - default to random numeric ID)
      Unique identier for the simulation run for generating
      multiple independent trajectories for input to varScaffold

    --simProc "integer" (optional - default is to run on all processors)
      number of processors to run NAMD simulation

    --singleRun (optional)
      flag that will generate PDB trajectories from the results.
      Used to run SNP2SIM
#### Config template

```
# varMDsim options
##################

 #path to cleaned PDB file (protein structure w/ cannonical aa)
newStruct:

 #varTraj simulation length in ns
simLength:

 #amino acid to change
simID:

 #number of processors to run simulation
simProc:

 #output summary PDB trajectory only
singleRun:

 #Only generate initial structures, without MD simulation
genStructures:
```

### varScaffold:
  Usage Notes:
  The config file specifies the alignment and clustering parameters

  Module Options
  
    --scaffID "string"
      Specifies naming of scaffolding output logfile containing
      cluster assignements for alignment and clustering parameters

    --clustPDBtraj (optional)
      if PDB clusters generated from varMDsim --singleRun are used as input
      this flag must be set and trajectory files specified with --loadPDBtraj

    --loadPDBtraj "file path(s)" (optional)
      paths to PDB trajectory files, separated by single space

#### Config template

```
# varScaffold
#############

 #unique identifier for run
scaffID:

 #Use PDB trajectories (possibly from singleCGC runs)
clustPDBtraj:

 #pdb trajectory files to import, list one after another
loadPDBtraj:

#Cluster Parameters - follow VMD atomselection format
 #Residues to superimpose trajectory 
alignmentResidues:
 #Residues to consider when clustering trajectory
clusterResidues:

 #use a table of structural features to perform clustering, otherwise use a pairwise distance matrix
featureTableMethod: 
 #if mds method used, path to pairwise RMSD table, otherwise automatically generated
PairwiseRMSD: 

 #create a script that colors the cluster residues in VMD by cluster
colorTrajectory: 
```

### drugSearch:
  Usage Notes:
  This module can use a manually defined search space provided to --newBindingConfig
  defined in reference to a template structure (original input PDB)
  Values can be determined using methods specified in the AutoDock Tutorial

  Binding Config:

  "center_x = ""x coordinate""

  "center_y = ""y coordinate""

  "center_z = ""z coordinate""

  "size_x = ""dimension of x""

  "size_y = ""dimension of y""

  "size_z = ""dimension of z""

  This search space can also be automatically determined given a list of residues to include in the search space.
  

  Additionally, the flexable residues can be supplied to --flexBinding as
  a sting of integers corresponding to the numeric residue ID of
  binding residues

  Multiple drugs can be bound by copying ligand PDBQT files into
  ./snp2sim/drugLibraries/

  Module Options:
    --bindingTemplate "file path"
      path to PDB file used to specify AutoDock Vina search space

    --newBindingConfig "file path"
      path to config file defining search space for AutoDock Vina

    --flexBinding "file path"
      path to file containing a line of intergers corresponding to
      the residue ID of the binding residues.

    --drugLibrary "string"
      Name of drug library in snp2sim/drugLibraries or identifier for
      ligands input through the command line

    --singleDrug "file path"
      Path to PDBQT formated ligand file

    --inputScaff "file path(s)" (optional)
      path to PDB files of variant specific scaffolds

    --vinaExh "integer" (optional - default 50)
      used to specifiy exhaustiveness of AutoDock Vina search
  
#### Config template

```
 #path to PDB file used to create search space to align scaffold
bindingTemplate:

 #exhaustiveness parameter of autodock vina
vinaExh:

 #path to new binding config file
newBindingConfig:

 #Automatically determine the search box based on the 
autoSearchSpace: 

 #list of residues that comprise the search space if autoSearchSpace is true, e.i 22 21 19 103 113 125
searchResidues:

 #list of residue numbers in flexible binding pocket, e.i 22 21 19 103 113 125
flexBinding:

 #name of snp2sim drug library
drugLibrary:

 #path to single drug PDBQT
singleDrug:

 #pdb scaffold files to import, list one after another
inputScaff:

 #only bind single variant scaffolds
bindSingleVar:

 #path to dir with ligand PDBs, which will be converted to drug library
ligandPDB:

 #number of times to run the docking simulation, to get an uncertainty measurement of the binding energy
numTrials:
```

### varAnalysis:
  Usage Notes: The input to the analysis module is a list of variants to analyze, as well as optionally specific drugs or drug libraries to include.
  
#### Config template

```
# varAnalysis
############

 #Variants to be used in variants (must have a results/variant/drugBinding directory with .pdbqt results)
 #Format ###AA (### = residue number and AA = variant amino acid) or wt
 #must include wt
analysisVariants: 
  - list
  - of
  - variants

 #Analyze results from specific drug libraries
analysisDrugLibrary: 

 #Analyze results from a specific ligand
analysisDrug:

 #Enable if using scaffolding results from old versions of SNP2SIM
legacyScaff:

 #Enable if using drugBinding results from custom PDB proteins, and not varScaffold results
customScaff: 

```

