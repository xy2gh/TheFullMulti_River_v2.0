# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:58:47 2020

@author: A. Praetorius and PradoDomercq
"""
'''
Author: Xiaoyu Zhang (xiaoyu.zhang@aces.su.se)
Date: 2023-12-12 10:21:29
LastEditTime: 2024-09-02 15:51:24
Description: This file was modified for updating the new parameterization and
             calculating the rate constants of modeling MPs in the aquatic environment.
'''


#extension of RS_generator module containing functions to calculate all rate constants 
#Modification of advection and addition of sed transport for rivers

import math
import pandas as pd

#import file storing required constants
from helpers.GlobalConstants import *
from Inputs.GenerateGenericRiverImputFile import *

    
    
def degradation(t_half_d):
    
    #degradation estimations
    """ relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    #degradation half-life of MPs used as input is in days
        
    #degradation rate constant 
    k_deg = math.log(2)/(t_half_d*24*60*60) 

    # k_deg = 0
    # NOTE: turned off for other simulations
    
    return k_deg




    
def fragmentation(process_df,idx,particle,sizeBinIdx,aggState):

    t_frag_d=process_df.t_frag_d.loc[idx]

    if aggState== "A":#Free particel
        MP_radius_m=particle.radius_m
        MP_volume_m3=particle.volume_m3
        MP_diameter_um= particle.diameter_um
    elif (aggState == "C") or (aggState == "B"): #Biofouled or Heteroaggregated
        MP_radius_m=particle.parentMP.radius_m
        MP_volume_m3=particle.parentMP.volume_m3
        MP_diameter_um= particle.parentMP.diameter_um
    else: #Biofouled and heteroaggregated
        MP_radius_m=particle.parentMP.parentMP.radius_m
        MP_volume_m3=particle.parentMP.parentMP.volume_m3
        MP_diameter_um= particle.parentMP.parentMP.diameter_um


    if t_frag_d == "NAN":
        k_frag = 0
        fragments_formed=0
    else:
        if sizeBinIdx == "a" and aggState == "A":
            k_frag = 0
            fragments_formed=0
        else:
            volume_fragment = 4/3*math.pi*(MP_radius_m/10)**3 #!!!only works for bins 10 times smaller!!!
            fragments_formed = MP_volume_m3/volume_fragment
            k_frag = (1/(float(t_frag_d)*24*60*60))*MP_diameter_um/1000
             
    # k_frag = 0
    # fragments_formed=0
    # NOTE: turned off for other simulations
    
    return (k_frag, fragments_formed)

    
    


def settling(particle, compartment, comp_dict, settlingMethod):
    
    MP_radius_m=particle.radius_m
    MP_diameter_m = 2*MP_radius_m
    MP_volume_m3=particle.volume_m3
    
    MP_density_kg_m3=particle.density_kg_m3
    deltarho_g_m3 = (MP_density_kg_m3-density_w_21C_kg_m3)/1000
    
    #settling calculations
    #Settling occurs in all compartments but in sediment (compartment 4)
    if comp_dict[compartment].name == "sediment":
        k_set = 0
        
    elif deltarho_g_m3 < 0:
        k_set = 0
        
    else:
        
        if settlingMethod == "Stokes":
            vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
            k_set = vSet_m_s/comp_dict[compartment].depth_m
        
        else:
            print('Error: cannot pair with the right settling method (Stokes Law)')

            
    return k_set
 


        

def rising(particle, compartment, comp_dict, settlingMethod):
    
    MP_radius_m=particle.radius_m
    MP_diameter_m = 2*MP_radius_m
    MP_volume_m3=particle.volume_m3
    
    MP_density_kg_m3=particle.density_kg_m3
    deltarho_g_m3 = deltarho_g_m3 = (MP_density_kg_m3-density_w_21C_kg_m3)/1000
    
    # Reynolds number - Re
    Re = vflow_m_s*MP_diameter_m/v_w_21C_m2_s
    
    
    #rising calculations
    #Rising only occus in the flowing water and stagnant water compartments (2 and 3)
    if (comp_dict[compartment].name == "surface") or (comp_dict[compartment].name == "sediment"): 
        k_rise = 0
        
    elif deltarho_g_m3 > 0:
        k_rise = 0
        
    else:
        
        if settlingMethod == "Stokes":
            vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
            k_rise = -vSet_m_s/comp_dict[compartment].depth_m
        
        else: 
            print('Error: cannot pair with the right settling method (Stokes Law)')

        
    return k_rise





def heteroagg(process_df,idx,particle,SPM1,G,T_K,compartment,aggState,settlingMethod):
    alpha=process_df.alpha.loc[idx]
    
    # MP info
    MP_radius_m=particle.radius_m
    MP_diameter_m = 2*MP_radius_m
    MP_volume_m3=particle.volume_m3
    MP_density_kg_m3=particle.density_kg_m3
    MP_deltarho_g_m3 = (MP_density_kg_m3-density_w_21C_kg_m3)/1000
    # Reynolds number - Re
    MP_Re = vflow_m_s*MP_diameter_m/v_w_21C_m2_s
    
    # SPM info    
    SPM_radius_m=SPM1.radius_m
    SPM_density_kg_m3=SPM1.density_kg_m3
    SPM_concNum_part_m3=SPM1.concNum_part_m3


    # Calculating the settling velocity
    if settlingMethod == "Stokes":
        MP_vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
        
    else: 
        print('Error: cannot pair with the right settling method (Stokes Law)')
    
    # Default SPM settling velocity based on Stokes Law
    SPM_vSet_m_s = 2/9*(SPM_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(SPM_radius_m)**2    
    
    #Heteroaggregation only occurs for pristine (A) and biofouled MPs (C) and on the water compartments (1, 2 and 3)
    #heteroaggregation for B and D is limited by alpha values given as NA    
    if (compartment == "4") or (aggState == "B") or (aggState =="D"):
        k_hetAgg = 0
    else:
    
        #heteroaggregation rate constants
        """heteroaggregation requires to particles to collide and interact
        favorably for the collision to result in attachment
        the heteroaggregation rate constants is therefore composed of two
        parts, 1) a collision rate constant and 2) and attachment 
        efficiency (alpha) (representing the probability of attachment).
        For heteroaggregation a common simplifaction is the assumption that
        SPM concentration is not signficantly affected by the heteroaggre-
        gation process. Therefore, a pseudo first-order heteroaggregation 
        rate constant is obtained by multiplying collision rate with alpha 
        and with the SPM number concentration"""
        
        #first the different collision mechanisms are calculated
        k_peri = (2*k_B_J_K*T_K)/(3*mu_w_21C_kg_ms)*(MP_radius_m + SPM_radius_m)**2/(MP_radius_m * SPM_radius_m)
        #perikinetic contributions to collision rate constant (Brownian motion)
        
        k_ortho = 4/3*G*(MP_radius_m + SPM_radius_m)**3
        #orthokinetic contributions to collision rate constant (caused by fluid motion)
        
        k_diffSettling = math.pi*(MP_radius_m + SPM_radius_m)**2 * abs(MP_vSet_m_s-SPM_vSet_m_s)
        #differential settling contributions to collision rate constant
    
        k_coll = k_peri + k_ortho + k_diffSettling
        #the collision rate constant
        
        k_hetAgg = alpha*k_coll*SPM_concNum_part_m3
        #the pseudo first-order heteroaggregation rate constant

    # k_hetAgg = 0
    # NOTE: turned off for other simulations
    
    
    return k_hetAgg





def breakup(process_df,idx,particle,SPM1,G,T_K,compartment,aggState,settlingMethod):

    # MP info
    MP_radius_m=particle.radius_m
    MP_diameter_m = 2*MP_radius_m
    MP_volume_m3=particle.volume_m3
    MP_density_kg_m3=particle.density_kg_m3
    MP_deltarho_g_m3 = (MP_density_kg_m3-density_w_21C_kg_m3)/1000
    # Reynolds number - Re
    MP_Re = vflow_m_s*MP_diameter_m/v_w_21C_m2_s
    
    # SPM info    
    SPM_radius_m=SPM1.radius_m
    SPM_density_kg_m3=SPM1.density_kg_m3
    SPM_concNum_part_m3=SPM1.concNum_part_m3
    
    
    # Calculating the settling velocity
    if settlingMethod == "Stokes":
        MP_vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
    
    else: 
        print('Error: cannot pair with the right settling method (Stokes Law)')
    
    # Default SPM settling velocity based on Stokes Law
    SPM_vSet_m_s = 2/9*(SPM_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(SPM_radius_m)**2    
    
    
    # Breackup doesnt occur in the sediment compartment and only for MP aggregates (B and D), 
    #however Kbreackup is calculated based on Kheter of A and C
    
    if (compartment == "4") or (aggState == "A") or (aggState =="C"):
        k_aggBreakup = 0
    else:
    
        #first the different collision mechanisms are calculated
        
        k_peri = (2*k_B_J_K*T_K)/(3*mu_w_21C_kg_ms)*(MP_radius_m + SPM_radius_m)**2/(MP_radius_m * SPM_radius_m)
        #perikinetic contributions to collision rate constant (Brownian motion)
        
        k_ortho = 4/3*G*(MP_radius_m + SPM_radius_m)**3
        #orthokinetic contributions to collision rate constant (caused by fluid motion)
        
        k_diffSettling = math.pi*(MP_radius_m + SPM_radius_m)**2 * abs(MP_vSet_m_s-SPM_vSet_m_s)
        #differential settling contributions to collision rate constant
    
        k_coll = k_peri + k_ortho + k_diffSettling
        #the collision rate constant
        
        k_hetAgg = process_df.alpha.loc[idx-1]*k_coll*SPM_concNum_part_m3
        #the pseudo first-order heteroaggregation rate constant
    
        k_aggBreakup = (1/10) * k_hetAgg
    
     
    # k_aggBreakup = 0
    # NOTE: turned off for other simulations

        
    return k_aggBreakup





def advection(compartments_prop,comp_dict, compartment,riverSection,river_flows):
     #advective transport
    
    # Based on Praetorius et al. 2012: Kflow = v_riv_flow*(Aw1/Vw1)
    #Being v_riv_flow the river flow velocity in ms-1, Aw1 is the crossectional 
    #area of the flowing water and Vw1 the volume of the box of moving water.
    #dimensions of the river we estimated resudence time of 26 days in flowing
    #water and 28 min in the surface watercompartment
    
    #RIVER SECTION DEPENDENT WITH VARYING DISCHARGE
    #calculate Cross sectional area of the flowing river
    depths=compartments_prop[compartments_prop.nameRS == "RS"+riverSection].depth_m
    RS_width=compartments_prop[compartments_prop.nameRS == "RS"+riverSection].width_m
    CrossAreaRS_m2=float(sum(depths[0:3])*(RS_width[0:1]))
    
    flow_df=river_flows[river_flows.Region_I == int(riverSection)+1]
    discharge_m3_s= pd.Series(flow_df["q(m3/h)"]/60/60)
    
    if comp_dict[compartment].name == "flowingWater" or comp_dict[compartment].name == "surface":
        k_adv_series = discharge_m3_s*(comp_dict[compartment].CrossArea_m2/CrossAreaRS_m2)/ comp_dict[compartment].volume_m3
        k_adv = tuple(k_adv_series)
        
    else:
        k_adv = 0

    
    return k_adv





def mixing(flowingWater, compartment, comp_dict, vflow_m_s): 
    
    diffusity_m2_s = 1*10**-3 #ref(1)
    Avertical_m2 = lengthTotal_km *1000/numRS * widthRS
    
    k_mix_up = diffusity_m2_s/Avertical_m2 * (vflow_m_s/10**-3)**2
    k_mix_down = diffusity_m2_s/Avertical_m2 * (vflow_m_s/10**-3)**2
    
        
    if comp_dict[compartment].name == "flowingWater":
            k_mix=(k_mix_up,k_mix_down)

    elif comp_dict[compartment].name == "surface":
            k_mix = k_mix_up*(flowingWater.volume_m3)/(comp_dict[compartment].volume_m3)
        
    elif comp_dict[compartment].name == "stagnantWater":
            k_mix = k_mix_down*(flowingWater.volume_m3)/(comp_dict[compartment].volume_m3)
            
    elif comp_dict[compartment].name == "sediment":
            k_mix = 0

    
    return k_mix

''' 
    k_mix is used for the mixing process between flowing & surface water and flowing & stagnant water, 
    conducting a turbulent vertical coefficient in rivers (1).
        
    References:
    (1): <Handbook on Mixing in Rivers> Edited by J.C. Rutherford
    (Water and Soil Miscellaneous Publication No. 26. 1981. 60pp.ISSN 0110-4705)
'''





def biofilm(compartment, process_df, comp_dict, idx, aggState):
    #Biofilm formation taken into account for pristin and heteroaggregated MPs (A and B)
    #only takes place in the water compartments ( 1, 2 and 3)
    
    if (aggState == "A") or (aggState=="B"):
    
        if comp_dict[compartment].name == "surface":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
                
        if comp_dict[compartment].name == "flowingWater":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
            
        if comp_dict[compartment].name == "stagnantWater":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
            
        if comp_dict[compartment].name == "sediment":
                k_biof = 0
    else:
        k_biof = 0 
            
        #assume it takes x days for biofilm coverage to grow
        #need to update!!    # k_biof = 1/30/24/60/60 #assume it takes 30 days for biofilm coverage to grow
    # #need to update!!
    

    # k_biof = 0 
    # NOTE: turned off for other simulations


    return k_biof





def defouling(compartment, process_df, comp_dict, idx, aggState):
    #Defouling = degradation of Biofilm. for biofouled and heteroaggregated and biofouled particles (C and D)
    #only takes place in the water compartments ( 1, 2 and 3)
    
    if (aggState == "C") or (aggState=="D"):
        if comp_dict[compartment].name == "sediment":
            k_defoul = 0
    
        else:
            if type(process_df.t_biof_degrad_d.loc[idx]) == str:
                k_defoul = 0
            else:
                k_defoul = 1/process_df.t_biof_degrad_d.loc[idx]/24/60/60

    else:
        k_defoul = 0 
        #assume it takes x days for biofilm coverage to be degraded #k_defoul = 0  # for mass balance test
    
    
    # k_defoul = 0
    # NOTE: turned off for other simulations
    
    return k_defoul





#for the sediment compartment rate constants for resuspension and
# burial in deep sediment are calculated & degradation rate assigned

def resusp(particle, compartment, comp_dict, vflow_m_s):
    #MP_radius_m=particle.radius_m
    #vFlow_m_s=river_flows[river_flows.Region_I == int(riverSection)]
    
    MP_radius_m=particle.radius_m

    if comp_dict[compartment].name == "sediment":
        k_resusp = (5.0 * 10**-6) * (vflow_m_s/ 10**-3)**2 / comp_dict[compartment].depth_m
             
    else:
        k_resusp = 0

    
    return k_resusp





def burial(particle, compartment, comp_dict, settlingMethod, vflow_m_s):
    #MP_radius_m=particle.radius_m
    #vFlow_m_s=river_flows[river_flows.Region_I == int(riverSection)]
    
    MP_radius_m=particle.radius_m
    MP_diameter_m = 2*MP_radius_m
    MP_volume_m3=particle.volume_m3
    
    MP_density_kg_m3=particle.density_kg_m3
    deltarho_g_m3 = (MP_density_kg_m3-density_w_21C_kg_m3)/1000
    
    
    #settling calculations
    if comp_dict[compartment].name == "sediment":
                
        # calculating the resuspension
        k_resusp = (9 * 10**-5) * (vflow_m_s/ 10**-3)**2 / comp_dict[compartment].depth_m           
        
    
        if deltarho_g_m3 < 0:
            k_burial = 0
            
        else:
            if settlingMethod == "Stokes":
                vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
                
                k_burial = vSet_m_s/compDepth[2] - k_resusp
                
                if k_burial < 0:
                    k_burial = 0
                
    else:
        k_burial = 0

    
    return k_burial





def sedTransport(compartment, comp_dict):
    #MP_radius_m=particle.radius_m
    #vFlow_m_s=river_flows[river_flows.Region_I == int(riverSection)]
    
    if comp_dict[compartment].name == "sediment":
        
        #!!! Rolling and Sliding transport
        m_sed_kg=(1-sed_porosity)*sed_density*(10**3)*comp_dict[compartment].volume_m3
        k_sed_trans=v_sed_trans/m_sed_kg
                    
    else:
        k_sed_trans=0
        
    return k_sed_trans





def diffusion(particle, compartment, comp_dict): 
    
    MP_radius_m=particle.radius_m
    D_parameter = k_B_J_K*float(T_K)/(6*math.pi*mu_w_21C_kg_ms*MP_radius_m)

    if comp_dict[compartment].name == "stagnantWater": 
        #Stokes–Einstein–Sutherland equation
        k_diffusion = D_parameter*comp_dict[compartment].depth_m/comp_dict[compartment].volume_m3
    
    if comp_dict[compartment].name == "sediment":
        k_diffusion = D_parameter*comp_dict[compartment].depth_m/comp_dict[compartment].volume_m3 * compDepth[2]/compDepth[3]
    
    else:
        k_diffusion = 0

        
    return k_diffusion