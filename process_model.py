# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:24:40 2020

@author: nkruyer3
"""
#Unit operation models for cyanobacteria to 2,3-BDO bio-ISRU process

import numpy as np
import phasepy as pp #Phasepy package from Chapparo & Mejia, 2020
from scipy.integrate import solve_ivp
from scipy.optimize import root

#Biomass productivity model equation for suspended growth - from Karemore et al, 2020
def prod_sus(F,a,Ek,y,E0,k,D,R0,Cc,t1,t2):
     p1 = a*Ek*y
     p2 = np.exp(-k*D)
     p3 = (Ek + E0)/(Ek + E0*p2)
     p4 = t1/D
     p5 = R0*Cc*y*t2
     
     return F*((p1*np.log(p3)*p4) - p5)
 
#Biomass productivity model equation for biofilm growth - adapted from Karemore et al, 2020
def prod_biofilm(F,a,Ek,y,E0,k,D,R0,Cc,t1,t2):
    p1 = a*Ek*y
    p2 = np.exp(-k*D)
    p3 = (Ek + E0)/(Ek + E0*p2)
    p4 = t1
    p5 = R0*Cc*y*t2
    
    return F*((p1*np.log(p3)*p4) - p5) 

#Calculate temperature dependence of Ek and R0 productivity model parameters - from Karemore et al, 2020
def Ek_temp(T):
    K = T+273.15
    return 5228028746771.481*np.exp(-60000/(8.3145*K))
def R0_temp(T):
    K = T+273.15
    return 4497.016238*np.exp(-27000/(8.3145*K))/60

#Calculates productivity for cyanobacteria grown in suspension
#Returns total grams of biomass produced per day (g_day) as well as productivity in g/m^2 of land area/day (pr_gmd)
#input = cyano cultivation required land, F ratio, width and height of PBR
#Can calculate productivity without land_area
def cyano_cultivation_sus(land_area, F, PBR_w, PBR_h):
    T = 25; a = 0.061; y = 22.2e-6; 
    E0_M = 142.8; k = 175; 
    Cc = 18000; t1 = 43200; t2 = 84600;
       
    pr = prod_sus(F,a,Ek_temp(T),y,E0_M,k,PBR_w,R0_temp(T),Cc,t1,t2) #gives productivity in g/m^3/day
    v = 48*PBR_w*PBR_h*24000/1250 #calculates L/m^2 of land area, based on 24 PBRs/1250 m^2
    pr_gmd = (pr/1000)*v #converts calculated productivity into g/m^2/day
    g_day = pr_gmd*land_area
    
    return g_day, pr_gmd


#Calculates productivity for cyanobacteria grown in a biofilm
#Returns total grams of biomass produced per day (g_day) as well as productivity in g/m^2 of land area/day (pr_gmd)
#input = cyano cultivation required land, F ratio, biofilm thickness
#Can calculate productivity without land_area
def cyano_cultivation_film(land_area, F, film_depth):
    T = 25; a = 0.061; y = 22.2e-6; 
    E0_M = 142.8; k = 175*50; 
    Cc = 135; t1 = 43200; t2 = 84600;
    
    pr_gmd = prod_biofilm(F,a,Ek_temp(T),y,E0_M,k,film_depth,R0_temp(T),Cc,t1,t2) #productivty in g/m^2/day
    g_day = pr_gmd*land_area
    
    return g_day, pr_gmd

#Models biomass concentrator unit operation as a function of feed_flow, feed concentation, cross flow membrane flux
#and desired outlet concentration
#Returns flowrate of concentrated cyanobacteria, pure water and required cross flow membrane area
def biomass_concentrator(feed_flow, feed_conc, flux, out_conc):
    dil_fact = out_conc/feed_conc
    cy_flow_out = feed_flow/dil_fact
    H2O_flow_out = feed_flow - cy_flow_out
    mem_area = H2O_flow_out*60/flux
    
    return cy_flow_out, H2O_flow_out, mem_area

#Models enzyme digester based on inlet cyanobacteria concentration and glucose yield
#Returns concentration of glucose out of the digester (g/L)
def enzyme_digester(cy_conc, gluc_yield):
    gluc = cy_conc*gluc_yield
    
    return gluc #g/L

#Models E. coli fermenters based on inlet glucose concentration and 2,3-BDO yeidl
#Returns concentration of 2,3-BDO out of digester (g/L)
def fermenter(gluc_conc, BDO_yield):
    BDO = gluc_conc*BDO_yield
    
    return BDO #g/L

#Converts flows and concentrations to mole fraction for input into LLE_unit
#Inputs are flowrates of H2O, 2,3-BDO and butanol
#Returns mole fraction of each component
def mol_frac(H2O,BDO,BUT):
    
    #total moles of each component
    mol_H2O = H2O*(1/18) #assuming H2O input as g/min - converts to mol/min
    mol_BDO = BDO*(1/90.121) #assuming BDO input as g/min - converts to mol/min
    mol_BUT = BUT*(1/74.121) #assuming BUT input as g/min - converts to mol/min
    
    tot_mol = mol_H2O + mol_BDO + mol_BUT
    
    #mole fraction of each component
    y_H2O = mol_H2O / tot_mol
    y_BDO = mol_BDO / tot_mol
    y_BUT = mol_BUT / tot_mol
    
    return y_H2O, y_BDO, y_BUT

#Converts mole fractions into mass flow rates
#Inputs of flowrate (L/min) and mole fraction vector [y_H2O, y_2,3-BDO, y_butanol]
#Output of mass flowrate (g/min)
def flow_conc(flow, mol_frac):
    #turns from mol_frac list into concentration (g/L)
    h2o = mol_frac[0]
    bdo = mol_frac[1]
    but = mol_frac[2]
    
    tot_g = h2o*18 + bdo*90.121 + but*74.121
    
    x_h2o = (h2o*18)/tot_g
    x_bdo = (bdo*90.121)/tot_g
    x_but = (but*74.121)/tot_g
    
    dens = x_h2o*1000 + x_bdo*1000 + x_but*810
    mass_flow_tot = dens*flow
    
    mass_flow = [mass_flow_tot*x_h2o, mass_flow_tot*x_bdo, mass_flow_tot*x_but]
    
    return mass_flow

#Simluates the LLE unit operation using Phasepy python package

def LLE_unit(Z):
   
    #Defines physical properties of H2O, 2,3-BDO and butanol in Phasepy package
    H2O = pp.component(name = 'H2O', Tc = 647.13, Pc = 220.55, Zc = 0.229, Vc = 55.948, w = 0.344861, Ant = [11.64785144, 3797.41566067,-46.7783044], GC = {'H2O':1})
    BDO = pp.component(name = 'BDO', Tc = 639.26, Pc = 50.875, Vc = 290, w = 1.12, Ant = [6.07439, 2616.746, -24.565], GC = {'CH3':2, 'CH':2, 'OH(P)':2})
    BUT = pp.component(name = 'BUT', Tc = 563.1, Pc = 44.2, Vc = 275, w = 0.593 , Ant = [4.54607, 1351.555,-93.34], GC = {'CH3':1, 'CH2':3, 'OH(P)':1})
    
    #binary interaction parameters regressed in AspenTech
    kij = np.array([[0,0.07254,0.1261],[0.07254,0,0.0233324],[0.1261,0.0233324,0]])

    mix = pp.mixture(H2O, BDO)
    mix.add_component(BUT)
    mix.unifac()
    mix.kij_cubic(kij)

    eos = pp.rkseos(mix, mixrule = 'qmr')

    x0 = np.array([0.95,0.03,0.02])
    y0 = np.array([0.1, 0.1, 0.8])

    a = pp.equilibrium.lle(x0,y0,Z,298,1,eos)
    
    h2o_rich = a[0]
    but_rich = a[1]
    split = a[2]
    
    return h2o_rich, but_rich, split

#Defines flux equations for component flux through membrane in membrane sepration unit operation
def flux(A,F):
    H2O = F[0]
    BDO = F[1]
    BUT = F[2]
    
    l = 0.0005 #cm, membrane thickness
    p1 = 0.81 #density of 1-butanol
    
    X = BUT/(BUT + BDO + H2O) #mass fraction of BUT
    
    #permeability equations based on Shao et al (2009)
    Pb = 5*np.exp(0.8*X)*1e-8*60 #cm^2/min 
    Pbdo = 0.4*np.exp(X)*1e-8*60 #cm^2/min
    Pw = 20*(1 + X)*1e-8*60 #cm^2/min
    
    #correction to prevent negative flowrates - assumes concentrations < 1e-5 can be approximated as 0
    if H2O > 1e-5:
        dH2OdA = -(Pw/l)*(H2O/((BUT/p1) + BDO + H2O)) #g/cm^2/s
    else:
        dH2OdA = 0
    if BDO > 1e-5:
        dBDOdA = -(Pbdo/l)*(BDO/((BUT/p1) + BDO + H2O)) #g/cm^2/s
    else:
        dBDOdA = 0
    if BUT > 1e-5:
        dBUTdA = -(Pb/l)*(BUT/((BUT/p1) + BDO + H2O)) #g/cm^2/s
    else:
        dBUTdA = 0
    
    
    return dH2OdA, dBDOdA, dBUTdA

#Models membrane separation unit operation for 2,3-BDO concentration
#Solves differential equations defines in flux function
#returns mole flows of H2O, 2,3-BDO and butanol
def membrane_sep(A, C):
    #A = membrane area in m^2
    #C = concentration matrix [H2O, BDO, BUT] in g/min
    area = [0,A*100*100]
    F = solve_ivp(flux,area,C)
        
    return F
    #return y_h2o, y_bdo, y_but
    
#Uses previously defined functions for LLE and membrane separator, implementing
#system of equations simulating recyling of the butanol rich stream from the 
#membrane separator to the inlet of the LLE unit
def sep_equations(vars, *inputs):
    #Stream 1 = BDO input from fermenter, Stream 2 = fresh butanol input, Stream 3 = mixed recycle + fermenter + fresh but
    #Stream 4 = H2O rich stream from LLE, Stream 5 = But rich stream from LLE
    #Stream 6 = Recycle stream (but rich) from membrane, Stream 7 = product stream (bdo rich) from membrane
    
    x1w, x1b, x2t, x4w, x4b, x4t, x5w, x5b, x5t, x6w, x6b, x6t, x7w, x7b, x7t = vars 
    
    h2o_in, bdo_in, but_in, A = inputs
    
    eq1 = x1w - h2o_in #initial water flow - g/min
    eq2 = x1b - bdo_in #initial BDO flow - g/min
    eq3 = x2t - but_in #initial butanol flow - g/min
        
    x3w = x1w + x6w
    x3b = x1b + x6b
    x3t = x2t + x6t
    
    mass_in = x3w + x3b + x3t
    mol_in = (x3w/18) + (x3b/90.121) + (x3t/74.121)
    y3w, y3b, y3t = mol_frac(x3w, x3b, x3t)
    mol_frac_in = np.array([y3w, y3b, y3t])
    h2o_rich, but_rich, lle_split = LLE_unit(mol_frac_in)
    
    but_rich_mol = mol_in*lle_split
    but_rich_mass = but_rich_mol*(but_rich[0]*18 + but_rich[1]*90.121 + but_rich[2]*74.121)
    h2o_rich_mol = mol_in - but_rich_mol
    h2o_rich_mass = mass_in - but_rich_mass
    
    eq7 = x4w - h2o_rich[0]*h2o_rich_mol*18
    eq8 = x4b - h2o_rich[1]*h2o_rich_mol*90.121
    eq9 = x4t - h2o_rich[2]*h2o_rich_mol*74.121
    
    eq10 = x5w - but_rich[0]*but_rich_mol*18
    eq11 = x5b - but_rich[1]*but_rich_mol*90.121
    eq12 = x5t - but_rich[2]*but_rich_mol*74.121
    
    eq13 = x6w + x7w - x5w
    eq14 = x6b + x7b - x5b
    eq15 = x6t + x7t - x5t
    
    mem_out = membrane_sep(A, [x5w, x5b, x5t])
    mem_h2o = mem_out.y[0]
    mem_bdo = mem_out.y[1]
    mem_but = mem_out.y[2] 
    
    eq16 = x7w - mem_h2o[-1]
    eq17 = x7b - mem_bdo[-1]
    eq18 = x7t - mem_but[-1]
    
    return [eq1, eq2, eq3, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18]

#Solves system of equations defined in sep_equations function
#Returns stream flows and compositions for LLE unit, membrane separator and associated streams
def separations(h2o, bdo, but, area):
    x0 = [4740, 17.86, 8.1, 4740, 1, 1, 0.5, 15, 100, 10, 10, 10, 0.1, 20, 2]
    sep_out = root(sep_equations, x0, args = (h2o, bdo, but, area))
    
    return sep_out

#Models entire process for biofilm growth mode
#Input land requiremnt of cyano cultivation
#Output = total 2,3-BDO production, flowrate of concentrated cyanobacteria into the enzyme digester and stream results for separations
def BDO_production_sus(land):
    g_day, pr_gmd = cyano_cultivation_sus(land, 1.84, 0.045, 1)
    
    but_flow = 0.001
    mem_area = 0.5
    
    cy_conc = 1 #g/L
    concentrated_cy = 20 #g/L
    flow_rate = (pr_gmd*land)/(24*60)/(concentrated_cy/cy_conc) #L/mins
    glucose = enzyme_digester(concentrated_cy, 0.3)
    BDO = fermenter(glucose, 0.432)*flow_rate #g/min output
    
    h2o_flow = flow_rate*1000
    but = but_flow*810
    sep_out_arr = separations(h2o_flow, BDO, but, mem_area)
    
    sep_out = sep_out_arr.x
    
    purity = sep_out[13]/(sep_out[12] + sep_out[13] + sep_out[14])
    recovery = sep_out[13]/BDO
    production = sep_out[13]*60*24*500/1000000
    
    #flow of concentrated cyano into digester
    flow = sep_out[0]/1000 #(L/min)
    
    return production, flow, sep_out_arr

#Models entire process for biofilm growth mode
#Input land requiremnt of cyano cultivation, biomass, productivity, digester yield and fermenter yield
#Output = total 2,3-BDO production, flowrate of concentrated cyanobacteria into the enzyme digester and stream results for separations
def BDO_production_film(land, pr_gmd, enz_yield, ferm_yield):
    concentrated_cy = 20
    flow_rate = (pr_gmd*land)/(24*60)/(concentrated_cy) #L/min - idea is to harvest algae @ 20 g/L and not need a concentrator step
    glucose = enzyme_digester(concentrated_cy, enz_yield)
    BDO = fermenter(glucose, ferm_yield)*flow_rate #g/min output
    
    but_flow = 0.001
    mem_area = 0.5
    
    h2o_flow = flow_rate*1000
    but = but_flow*810
    sep_out_arr = separations(h2o_flow, BDO, but, mem_area) #input in g/time
    
    sep_out = sep_out_arr.x #flows are in g/min
    
    #also define purity and recovery % - not currently  used
    purity = sep_out[13]/(sep_out[12] + sep_out[13] + sep_out[14])
    recovery = sep_out[13]/BDO
    production = sep_out[13]*60*24*500/1000000
    
    #flow of concentrated cyano into digester
    flow = sep_out[0]/1000 #(L/min)
    
    return production, flow, sep_out_arr

#Determines cyanobacterial cultivation land requirement based on biological process specifications (cyano productivity, digester yield and fermentation yield)
#Requires starting guess for required land
def find_land(start_area, pr_gmd, enz_yield, ferm_yield):
    land = start_area
    prod = 5
    while prod < 15:
        land += 5000
        res = BDO_production_film(land, pr_gmd, enz_yield, ferm_yield)
        prod = res[0]
    land = land - 5000
    while prod < 15:
        land += 1000
        res = BDO_production_film(land, pr_gmd, enz_yield, ferm_yield)
        prod = res[0]
    land = land - 1100
    prod = 5
    while prod < 15:
        land += 100
        res = BDO_production_film(land, pr_gmd, enz_yield, ferm_yield)
        prod = res[0]
    land = land - 110
    prod = 5
    while prod < 15:
        land += 10
        res = BDO_production_film(land, pr_gmd, enz_yield, ferm_yield)
        prod = res[0]
    land = land - 11
    prod = 5
    while prod < 15:
        land += 1
        res = BDO_production_film(land, pr_gmd, enz_yield, ferm_yield)
        prod = res[0]
    print(land)
        
    return land, prod, res[1]