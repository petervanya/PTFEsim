# PTFEsim

PFTE simulation in LAMMPS
09/11/15

Collection of scripts to generate a data file of Nafion 
(and other PTFE membranes) with or without electrodes
for further DPD simulation

DPD parametrisation for Nafion, beads A to P
* A: CF2 CF2 CF2 CF2 CF2 CF2
* B: O CF2 C(CF3)F O CF2
* C: CF3 SO3H
* W: 6 H2O
* E: electrodes from carbon black, 54 carbon atoms per bead (optional)
* P: platinum beads (optional)

Default cell size: a few (5-20) DPD sizes (one = 8.14 AA)

## Workflow
1. Create `input.yaml` file with parameters using `default_ptfe_input.sh`
2. Run `gen_ptfe.py` with appropriate options to create LAMMPS data file
3. Run LAMMPS
Data file contains entries in SI units.

## Dependencies
* Numpy
* `sudo pip install docopt` for nice command line reading
