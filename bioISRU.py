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

#Returns water, power and mass requirement under modeled biomass productivity, digester yeild, fermentation and material densities
def wpm(biomass_productivity, digester_yield, fermentation_yield, PBR_density, reactor_density):
    land = find_land(5000, biomass_productivity, digester_yield, fermentation_yield) #can change initial guess as needed
    flow = BDO_production_film(land, biomass_productivity, digester_yield, fermentation_yield)[1]
    water, power, mass = specs(land, flow, biomass_productivity, digester_yield, fermentation_yield, PBR_density, reactor_density)
    
    return water, power, mass