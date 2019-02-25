# SNP2SIM
Molecular Simulation of Protein Structure Variants
Questions: Matthew McCoy - mdm299 <at> georgetown.edu

## Installation Instructions
The SNP2SIM workflow is a Python script which generates variant specific structural scaffolding for small molecule docking simulations.

Dependencies:
  Python
  Installed to PATH
    NAMD executable "namd2"
    VMD executatble "vmd"
    AutoDock Vina executable "vina"
    pythonsh (from AutoDockTools) executable "pythonsh"
  Alternatively, AutoDockTools scripts prepare_receptor4.py, prepare_flexreceptor4.py
  installed at /opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/

## Example Usage
Usage:python snp2sim.py <<options>>

  The workflow is configured to store intermediate files and
  results in a predefined directory structure. If the required
  trajectory/scaffold files are not present in the snp2sim directory,
  they must be speficied through the command line.

## Command Line Options:
### General:
--mode <string>
  select which snp2sim module to run
  appropriate values: varMDsim, varScaffold, drugSearch

--protein <string>
  user specified name for simulation
  used in naming output directories and file

### Optional:
--varResID <integer>
  Position in input PDB file to mutate to varAA

--varAA <character>
  Single letter code for one of 20 canonical amino acids to mutate varResID

--NAMDpath <string>
  path to NAMD executatble

--VMDpath <string>
  path to VMD executable

--VINApath <string>
  path to AutoDock Vina execuatable

--PYTHONSHpath <string>
  path to AutoDockTools script PYTHONSH

--ATDpath <string>
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
    --newStruct <file path>
      PDB file containing only the protein structure

    --simLength <real number>
      length of the NAMD simulation in nanoseconds

    --simID <string> (optional - default to random numeric ID)
      Unique identier for the simulation run for generating
      multiple independent trajectories for input to varScaffold

    --simProc <integer> (optional - default is to run on all processors)
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
  >alignementRes "<<VMD atomselect command>>"
  >clusterRes "<<VMD atomselect command>>"
  >rmsdThresh <<real number>>

  Module Options
    --newScaff <file path>
      path to config file with alignment and clustering parameters (see above)
      
    --scaffID <string>
      Specifies naming of scaffolding output logfile containing
      cluster assignements for alignment and clustering parameters

    --clustThresh <real number> (optional - default 0.09)
      The minimum porportion of the trajectory population required
      for scaffold identification

    --clustPDBtraj (optional)
      if PDB clusters generated from varMDsim --singleRun are used as input
      this flag must be set and trajectory files specified with --loadPDBtraj

    --loadPDBtraj <file path(s)> (optional)
      paths to PDB trajectory files, separated by single space


#### drugSearch:
  Usage Notes:
  This module requires the input of the search space parameters provided to --newBindingConfig
  defined in reference to a template structure (original input PDB)
  Values can be determines using methods specified in the AutoDock Tutorial
  Binding Config:
  >center_x = <<x coordinate>>
  >center_y = <<y coordinate>>
  >center_z = <<z coordinate>>
  >size_x = <<dimension of x>>
  >size_y = <<dimension of y>>
  >size_z = <<dimension of z>>

  Additionally, the flexable residues can be supplied to --flexBinding as
  a sting of integers corresponding to the numeric residue ID of
  binding residues

  Multiple drugs can be bound by copying ligand PDBQT files into
  ./snp2sim/drugLibraries/

  Module Options:
    --bindingTemplate <file path>
      path to PDB file used to specify AutoDock Vina search space

    --newBindingConfig <file path>
      path to config file defining search space for AutoDock Vina

    --flexBinding <file path>
      path to file containing a line of intergers corresponding to
      the residue ID of the binding residues.

    --drugLibrary <string>
      Name of drug library in snp2sim/drugLibraries or identifier for
      ligands input through the command line

    --singleDrug <file path>
      Path to PDBQT formated ligand file

    --inputScaff <file path(s)> (optional)
      path to PDB files of variant specific scaffolds

    --vinaExh <integer> (optional - default 50)
      used to specifiy exhaustiveness of AutoDock Vina search
  


### Analysis (in development)

