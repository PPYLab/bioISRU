# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:44:39 2021

@author: nkruyer3
"""
import numpy as np
from process_model import BDO_production_film
from process_model import membrane_sep as MEM

#Calculates water usage for biofilm process
def water_film(land, flow):
    PBR_water = land*1.84*2.18
    enz_water = flow*60*48
    ferm_water = flow*72
    
    return PBR_water, enz_water, ferm_water

#Calculates power usage for biofilm process
def power_film(land, flow, enz, ferm, prod, enz_yield, ferm_yield):
    PBR_power = 3.883e-5*land*1.84 #Source - Ozkan et al, 2012
    enz_power = enz/1000*1.5
    ferm_power = ferm/1000*1.5
    LLE_power = flow*60/1000*1.5
    mem_power = mem_heat(land,prod, enz_yield, ferm_yield)/60
    
    return sum([PBR_power, enz_power, ferm_power, LLE_power, mem_power])

#Calculates power usage in the membrane separator
def mem_heat(land, prod, enz_yield, ferm_yield):
    res = BDO_production_film(land, prod, enz_yield, ferm_yield)
    flows = res[2].x
    mem_out = MEM(0.5, flows[6:9])
    h2o_flow = mem_out.y[0,:]
    bdo_flow = mem_out.y[1,:]
    but_flow = mem_out.y[2,:]

    H_h2o = 2.4 #kJ/g
    H_bdo = 0.644 #kJ/g
    H_but = 0.739 #kJ/g

    h2o_heat = [-H_h2o*(h2o_flow[i]-h2o_flow[i-1]) for i in range(1,len(h2o_flow))]
    bdo_heat = [-H_bdo*(bdo_flow[i]-bdo_flow[i-1]) for i in range(1,len(bdo_flow))]
    but_heat = [-H_but*(but_flow[i]-but_flow[i-1]) for i in range(1,len(but_flow))]

    heat = sum(h2o_heat) + sum(bdo_heat) + sum(but_heat) #kJ/min
    
    return heat

#calculates payload mass for the biofilm growth mode
#input = cyano cultivation land, biomass productivity, volume of enzyme reactor, volume of fermenter...
#volume of lle unit, density of biofilm growth substrate, density of reactor material
#returns mass of each unit operation
def mass_film(land, prod, enz, ferm, lle, PBR_dens, reac_dens):
    PBR_mass = (land*1.84)*0.0003*100*100*100*PBR_dens/1000/1000 
    frame = frame_mass(land)
    DAP_mass = prod*500*land*0.002/22.2*132.06/1000000*1.2
    NH3_mass = ((prod*500*land*0.029/22.2*1.2)-(DAP_mass*1000000/132.06*2))*17.031/1000000
    TE_mass = 1.84*land*2.18*2.12837/1000000 #mass of trace elements, based on Cogne, 2002
    enz_mass = tank_mass(enz, reac_dens) + 3*enz*50/1000000 #adds in total mass of enzymes needed
    ferm_mass = tank_mass(ferm, reac_dens)
    LLE_mass = tank_mass(lle, reac_dens) + 0.001*60*24*500*810/1000000 #add mass of butanol needed
    mem_mass = 0.000623 #assume membrane mass does not fluctuate much, <<1% of total mass so small changes won't majorly impact total mass
    sto_tank_mass = 0.545 #2,3-BDO and O2 required = constant --> storage tank same'
    h2o_treat_mass = 0.008
  
    return PBR_mass, frame, DAP_mass, NH3_mass, TE_mass, enz_mass, ferm_mass, LLE_mass, mem_mass, sto_tank_mass, h2o_treat_mass

#Function to calcualte mass of reactor tank - based on common specificiations
def tank_mass(vol, density):
    v = (vol*1.2)/1000
    r = np.cbrt(v/(5*np.pi))
    h = r*5
    SA = 2*np.pi*r*h + 2*np.pi*r**2
    steel_v = SA*0.006
    mass = steel_v*density
    
    return mass

#Calculates mass of cyanobacterial support frame
def frame_mass(land):
    mass = 1.154*land/10000*1.38*0.38
    return mass

#Determines total water, power and mass usage for complete process using biofilm growth
#Inputs = cyano cultivation land, flow rate to enzyme digester, biomass productivity, digester yield, fermenter yield, ...
#biofilm growth substrate density, reactor material density
#Returns water in L, power in kW and mass in tons
def specs(land, flow, prod, enz_yield, ferm_yield, PBR_dens, reac_dens):
    PBR, enz, ferm = water_film(land,flow)
    water = sum([PBR, enz, ferm])
    power = power_film(land, flow, enz, ferm, prod, enz_yield, ferm_yield)
    mass = sum(mass_film(land, prod, enz, ferm, flow*60, PBR_dens, reac_dens))
    
    return water, power, mass