# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:37:08 2019

@author: Antonia Praetorius
"""
'''
Author: Xiaoyu Zhang (xiaoyu.zhang@aces.su.se)
Date: 2023-12-12 10:21:29
LastEditTime: 2024-06-28 15:24:50
Description: This file was modified for the global constants used in modeling.
'''


#this files stores all constants needed for running the model (i.e. for 
#particulate & environmental compartment objects & calcualting rate processes)

from Inputs.GenerateGenericRiverImputFile import * 


k_B_J_K = 1.38*10**-23 #Boltzmann constant k_B (in J/K)
g_m_s2 = 9.81 #gravitational acceleration on earth (in m/s2)

density_w_21C_kg_m3 = 998 #density of water at 21 degree C (in kg?m2)

mu_w_21C_mPas = 0.9764 #dynamic viscosity of water at 21 degree C (in mPa*s)
mu_w_21C_kg_ms = mu_w_21C_mPas/1000 #dynamic viscosity of water at 21 degree C (in kg/(m*s))


v_w_21C_m2_s = mu_w_21C_kg_ms/density_w_21C_kg_m3 #kinematic viscosity of water at 21 degree C (in m2/s)

tc = 0.05 # critical shear stress (tc): empirical value (~ 0.03 to 0.06), here set as 0.05


### Sediment parameters
#v_sed_trans (kg/s)
#sed_porosity (dimensionless)
#sed_density (g/cm3)

if imputRiver == 'Generic':
    v_sed_trans = 3.0
    sed_porosity = 0.85
    sed_density = 2.5
