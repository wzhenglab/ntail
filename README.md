## Coarse-grained model for 
## "Transient Non-local Interactions Dominate the Dynamics of Measles Virus NTAIL"

Simulation_with_HOOMD.py: script for coarse-grained simulations, requiring HOOMD-Blue v2.9.3 and azplugins

usage: python Simulation_with_HOOMD.py $1 $2 $3 $4 $5
  - $1: sequence file name
  - $2: salt concentration in unit of M
  - $3: pairwise interaction strength epsilon 
  - $4: temperature in unit of K
  - $5: harmonic dihedral potential strength

e.g: python Simulation_with_HOOMD.py ntail 0.15 0.16 298 0.2

