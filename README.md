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

## Example Usage
Usage:python snp2sim.py ""options""

The workflow is configured to store intermediate files and
results in a predefined directory structure. If the required
trajectory/scaffold files are not present in the snp2sim directory,
they must be speficied through the command line.
  
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
If NAMD output exists in the file structure defined by the 
--protein --varResID and --varAA command line options.

The input can also be a list of PDB formatted trajectory files using 
the --clustPDBtraj command line option.

Additionally, a clustering config is required. See the template file
in the example directory for the format.

#### Command line options

`python snp2sim.py --mode varScaffold --protein PDL1 --varResID 115 --varAA T  --newScaff example/pdl1ScaffConfig.0.7.txt --scaffID bindingRes `

#### Output files

The output files include a representative PDB for each cluster, 
as well as a log file with the cluster assignments for all 
trajectory structures are stored in the 
"variantSimulations/PROTEIN/VARIANT/scaffold" directory.

### Generating Small Molecule Docking Results using drugSearch
The drugSearch module requires the most user preparation,
specifically using AutoDockTools to define the search space
on a reference structure (typically the initial structure 
used to generate the variant scaffolds.)

#### Input Files
The drugSearch module takes a list of PDB formatted structures, 
as well as a reference structure used to define the search space 
location and dimensions (provided in the config). The user can also specify 
residues to use 

#### Command line options

`python snp2sim.py --mode drugSearch --protein PDL1 --varResID 115 --varAA T 
--newScaff example/pdl1ScaffConfig.0.7.txt --scaffID bindingRes
--bindingTemplate example/PDL1.Vtype --newBindingConfig example/pdl1SearchSpace.txt 
--flexBinding example/pdl1FlexRes.txt --drugLibrary pdl1-SMI 
--inputScaff example/exampleResults/PDL1.115T.bindingRes.cl1.scaffold.pdb example/exampleResults/PDL1.115T.bindingRes.cl2.scaffold.pdb

#### Output files
The log files containing the binding affities for the top
ligand poses, as well as the PDBQT files that contain the 
coordinates of the ligand and flexable residues are depositied in the 
"variantSimulations/PROTEIN/results/VARIANT/drugBinding" directory.

Additionally, the PDBQT files used for the AutoDock Vina
simulations can be found in the 
"variantSimulations/PROTEIN/results/VARIANT/scaffold" directory.
  
## Command Line Options:
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

### Options specific to snp2sim modules
#### varMDsim:
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


#### varScaffold:
  Usage Notes
  The imput config file specifies the alignment and clustering parameters
  It is recommend that multiple itterations with different RMSD thresholds
  are run to determine the value which optimizes the number of variant scaffolds.
  
  Text file formatted as:
  "alignementRes """VMD atomselect command"""
  "clusterRes """VMD atomselect command"""
  "rmsdThresh ""real number""

  Module Options
    --newScaff "file path"
      path to config file with alignment and clustering parameters (see above)
      
    --scaffID "string"
      Specifies naming of scaffolding output logfile containing
      cluster assignements for alignment and clustering parameters

    --clustThresh "real number" (optional - default 0.09)
      The minimum porportion of the trajectory population required
      for scaffold identification

    --clustPDBtraj (optional)
      if PDB clusters generated from varMDsim --singleRun are used as input
      this flag must be set and trajectory files specified with --loadPDBtraj

    --loadPDBtraj "file path(s)" (optional)
      paths to PDB trajectory files, separated by single space


#### drugSearch:
  Usage Notes:
  This module requires the input of the search space parameters provided to --newBindingConfig
  defined in reference to a template structure (original input PDB)
  Values can be determines using methods specified in the AutoDock Tutorial
  Binding Config:
  "center_x = ""x coordinate""
  "center_y = ""y coordinate""
  "center_z = ""z coordinate""
  "size_x = ""dimension of x""
  "size_y = ""dimension of y""
  "size_z = ""dimension of z""

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
  



