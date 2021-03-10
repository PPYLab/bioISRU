# bioISRU
Repository of code used in process modeling of biotechnology enabled ISRU for Martian biofuel

Scripts contained in this repository were used for process modeling as well as flux balance analysis

Process Modeling:
process_model.py - defines functions for modeling and simulation of process unit operations
proces_specifications.py - uses process_model to determine process specifications, namely water, power and mass usage
bioISRU.py - provides sample script for implementation of process model

Flux Balance Analysis:
code to be implemented with COBRA 3.0 Toolbox (https://opencobra.github.io/cobratoolbox/stable/)
code utilized E. coli genome scale model iML1515 (http://bigg.ucsd.edu/models/iML1515)

add12PDO.m - adds pathway for 1,2-propanediol production to genome scale model
add13BDO.m - adds pathway for 1,3-butanediol production to genome scale model
add23BDO.m - adds pathway for 2,3-butanediol production to genome scale model
fuel_theoretical_yields.m - calculates maximum theoretical yields for 1,2-PDO, 1,3-BDO, and 2,3-BDO
