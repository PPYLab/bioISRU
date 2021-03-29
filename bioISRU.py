# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 16:43:10 2021

@author: nkruyer3
"""
#Example script for going from desired process specifications to determining water, power and mass requirements of the
#process - for biofilm growth mode

from process_model import find_land
from process_specifications import specs
from process_model import BDO_production_film
import pandas as pd

#Returns water, power and mass requirement under modeled biomass productivity, digester yeild, fermentation and material densities
def wpm(biomass_productivity, digester_yield, fermentation_yield, PBR_density, reactor_density):
    land = find_land(5000, biomass_productivity, digester_yield, fermentation_yield)[0] #can change initial guess as needed
    flow = BDO_production_film(land, biomass_productivity, digester_yield, fermentation_yield)[1]
    water, power, mass = specs(land, flow, biomass_productivity, digester_yield, fermentation_yield, PBR_density, reactor_density)
    
    return water, power, mass, land

#Load sample data set into Pandas dataframe
sample = pd.read_csv('Sample_data_Peralta-Yahya.csv', nrows = 9)

#calculate specifications (s) for sample data loaded into 'sample' dataframe
s = [wpm(sample['Productivity'][i], sample['Digester'][i], sample['Fermenter'][i], sample['PBR'][i], sample['Reactor'][i]) for i in range(len(sample))]

water = [s[i][0] for i in range(len(s))]
power = [s[i][1] for i in range(len(s))]
mass =  [s[i][2] for i in range(len(s))]
land =  [s[i][3] for i in range(len(s))]

sample['Water'] = water
sample['Power'] = power
sample['Mass'] = mass
sample['Land'] = land

sample.to_csv('expected_demo_output.csv')