Repository of code used in process modeling of biotechnology enabled ISRU for Martian biofuel


Source code contained in this repository were used for process modeling


Process Modeling: 
process_model.py - defines functions for modeling and simulation of process unit operations 
process_specifications.py - uses process_model to determine process specifications, namely water, power and mass usage 
bioISRU.py - provides sample script for implementation of process model


Sample Dataset:
Sample_data_Peralta-Yahya.csv - contains sample inputs for determining the mass, power, water and land requirements for varying inputs of biomass productivity, digester yield, fermenter yield, and material density
Details about inputs
Productivity: Cyanobacteria biomass productivity (g/m^2/day)
Digester: Weight fractions of toral cyanobacteria into glucose 
Fermenter: g of 2,3-BDO /g glucose
PBR: Density of the biofilm material (PBR/ cotton)
Reactor: Density of the reactor (steel/ HDPE)


System Requirements:
Source code built and tested in Python 3.6 on Windows 10
Requires installation of Python 3.6 (or later)
Requires download and installation of Phasepy Python package version 0.0.49 or later (https://pypi.org/project/phasepy/)
Source code relies on Pandas, Numpy and Scipy packages


Installation Guide:
Once Python is installed, source code can be run without further installation


Demo:
Running the bioISRU.py script will give a demo on how the program runs
Output is a csv file (demo_output.csv) containing water, mass, power and land usage for input specifications
To run the demo, make sure that the Sample_data_Peralta-Yahya.csv file is in the same working directory as the source code
The included csv file (expected_demo_output.csv) can be used to compare user's output with the expected source code results
Demo should complete in 1-2 hours


Instructions for Use:
To determine process specifications for varying inputs, the Sample_data_Peralta-Yayha.csv file can be edited to test for desired combinations of biomass productivity, digester yield, fermenter yield, and material density
Deeper analysis of the process can be obtained using the process_specifications.py script to look at water, power and mass contributions for individual unit operations




Flux Balance Analysis: 
code to be implemented with COBRA 3.0 Toolbox (https://opencobra.github.io/cobratoolbox/stable/) in Matlab
Source code written and tested in Matlab 2017b
code utilized E. coli genome scale model iML1515 (http://bigg.ucsd.edu/models/iML1515)


add12PDO.m - adds pathway for 1,2-propanediol production to genome scale model 
add13BDO.m - adds pathway for 1,3-butanediol production to genome scale model 
add23BDO.m - adds pathway for 2,3-butanediol production to genome scale model 
fuel_theoretical_yields.m - calculates maximum theoretical yields for 1,2-PDO, 1,3-BDO, and 2,3-BDO