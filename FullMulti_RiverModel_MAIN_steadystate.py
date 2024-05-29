# -*- coding: utf-8 -*-

"""
Author: Xiaoyu Zhang (xiaoyu.zhang@aces.su.se)
Date: 2024-02-22 16:18:40
LastEditTime: 2024-03-20 16:54:50
Description: 
    This file is modified from Prado Domercq's raw script for simulating the time varying situation.
    Running this script can gain the steady-state concentration of microplastic particles over distance. 
    
"""

# Import libraries and modules
from Functions.readImputParam import readProcessparam, microplasticData, readCompartmentData
from helpers.GlobalConstants import *
from helpers.helpers import *
from datetime import datetime, timedelta
from pathlib import Path
from Functions.fillInteractions_df_fun_v2_0 import*
from Functions.reshape_RC_df_fun import*
from Functions.RC_estimation_function_v2_XZ import*
from Functions.dilutionVol_calculator_func import*
from Functions.objectGenerationRiver_func import*
from matplotlib import ticker
import matplotlib.ticker as ticker
from cycler import*
from celluloid import Camera
import os
import pandas as pd
import itertools
import numpy as np
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from scipy.integrate import odeint
import warnings
warnings.filterwarnings('ignore')


# Import imput files
data_folder = Path("Inputs/")

process_df = readProcessparam(data_folder / "process_paramRiver.txt")
MP_prop = microplasticData(data_folder / "microplasticsSizeClass.txt")
compartments_prop = readCompartmentData(
    data_folder / "compartmentsGenericRiverSec_prop.txt")

river_flows = pd.read_csv(data_folder / "flow_connectivity.csv")

# Add river section depth field
RSdepth = []
for row in range(len(compartments_prop)):
    RSdepth.append(round(sum(compartments_prop.depth_m[0:4]), 2))
compartments_prop["depthRS_m"] = RSdepth


# MODEL SET UP

# RIVER COMPARTMENTS
compartments = ["Surface Water", "Flowing Water", "Stagnant Water", "Sediment"]
riverComp = ["1", "2", "3", "4"]
# MICROPLASTICS FORMS
MPforms = ["A", "B", "C", "D"]
MPformslabels = ["Free", "Heteroaggregated",
                 "Biofilm-covered", "Biofilm-heteroaggregated"]
# SIZE BINS
sizeBin = ["a", "b", "c", "d", "e"]
# Detection limit for MPs via Fourier Transform Infrared Spectroscopy is 20um
sizeBinLabel = ["0.1um", "1um", "10um", "100um", "1000um"]
# MPS RIVER PROCESSES (FATE AND TRANSPORT) LIST
processList = ["degradation", "fragmentation", "heteroagg", "breakup", "settling", "rising",
               "advection", "mixing", "biofilm", "resusp", "burial", "sedTransport", "defouling"]
# RIVER SECTIONS
numberRS = len(compartments_prop)/len(riverComp)
listRS = [*range(0, int(numberRS), 1)]
riverSect = [str(item) for item in listRS]
riverLengths = [str(it) for it in compartments_prop["length_m"]]
riverSectLength = riverLengths[0::4]
RS_cumLength_m = []
for d in range(len(riverSectLength)):
    if d == 0:
        RS_cumLength_m.append(float(riverSectLength[d]))
    else:
        RS_cumLength_m.append(
            float(riverSectLength[d])+float(RS_cumLength_m[d-1]))


# DEFINE RUN PARAMETERS
# - Solver (Dynamic or SteadyState).
# - mode (Standard). Monthly under developement
# - mode2 (Timelimit or raw): "Timelimit" mode sets up a time limit of 30min on the processes that exceeds that speed, while "raw" mode leaves the rate constant as calcualted. The raw version can straing the solver due to time.
# - record (True or False) : if "True" the results, RC and intercations dataframes will be recorded in the Results folder.
SOLVER = "SteadyState"
mode = "Standard"
mode2 = "raw"
record = "True"

# DEFINE SCENARIO
# - imputRiver: define wich river to study. (Generic)
# - composition: define MP composition to study. The composition is defined in the microplastics imput file microplasticsSizeClass.txt. Examples given in such file: PE, PA, PVC, and PS.
# - imputMP: define imput location and wich MP form is emmited by indicating the river section number, river compartment number, MP aggregation state and size bin: ex. 02Ae (RS=0:1, comp= 2:flowing water, MPtype:A:FreeMP, sizeBin:e:1000um)
# - imputFlow: define number of particles per minute entering the system
# - imputPulse: define number of particles entering the systems in a pulse (if any)
imputRiver = "Generic"

composition = "PE"
imputMP = "02Ae"

imputFlow = 60
imputPulse = 0


# Set simulation time:
# - t0: starting time (seconds)
# - daysSimulation: total length of simulation (days)
# - stepSize (seconds)
t0 = 0
daysSimulation = 365
tmax = 24*60*daysSimulation*60
sec_day = 24*60*60
stepSize = 60*60*24  # time step of 1day
timesteps = int(sec_day*daysSimulation/stepSize)

start_date = "2022-01-01"

t_span = np.linspace(0, tmax, int(timesteps)+1, dtype=int)
t_span_date = pd.to_datetime(t_span, unit="s", origin=pd.Timestamp(start_date))


# Define model run results file name
# Set up current date label

date_time_str = '2020-01-01 00:00'
DayStart = datetime.strptime(date_time_str, '%Y-%m-%d %H:%M')
LastDay = DayStart + timedelta(minutes=365*24*60)
date = DayStart

daterun = date.today()
daterun_label = daterun.strftime("%Y_%m_%d")

# Model run label
if imputFlow != 0:
    runtitle = daterun_label + "_" + composition + \
        "_inflow_" + imputMP+"_"+str(imputFlow)+"n_s_" + "in_River" + imputRiver + "_" + mode2
else:
    runtitle = daterun_label+"_" + composition + "_pulse_" + \
        imputMP+"_"+str(imputPulse)+"n_s_" + "in_River" + imputRiver + "_" + mode2
runtitle

# Generate list of species (combination of river section-compartment-MPform-sizeFraction)
# Generate COMBINATIONS
combinations = list(itertools.product(riverSect, riverComp, MPforms, sizeBin))
# Generate raw list of combinations and lists of concentrations (C) and inflows (I)
CombList = []
Ilist = []
Clist = []

def convertTuple(tup):
    str = ''.join(tup)
    return str

for e in combinations:
    Clist.append("C_" + convertTuple(e))
    Ilist.append("I_" + convertTuple(e))
    CombList.append(convertTuple(e))


# MODEL RUN

"""Estimate RATE CONSTANTS & HALF LIFE TIMES"""
# Estimate rate constants for all "species" = combination of RS*riverComp*MPforms*sizeBins (len (Clist))

RC_df = RC_estimation_function_v2_0(processList, CombList, Clist, MP_prop, compartments_prop, process_df, numberRS, composition, mode2, mode, date, riverComp, MPforms, sizeBin, river_flows)

# Reshape table of RC
RC_df_tidy = reshape_RC_df(RC_df, CombList)
RC_df_final = RC_df_tidy.pivot_table(
    index=["Compartment", "MP_form", "SizeFrac"], columns='Process', values='k_s-1', aggfunc=lambda x: x.tolist())

# Half live time dataframe
Half_life_df_h = RC_df_tidy.pivot_table(index=["Compartment", "MP_form", "SizeFrac"],
                                        columns=['Process'],
                                        values='t1/2_h')

interactions_df = fillInteractions_fun_v2_0(RC_df, Clist, river_flows)

# Heatmaps of half-life times
T_hm_array = [Half_life_df_h.loc["Surface Water"], Half_life_df_h.loc["Flowing Water"],
              Half_life_df_h.loc["Stagnant Water"], Half_life_df_h.loc["Sediment"]]
T_hm_array_ts = [x.transpose() for x in T_hm_array]
T_hm_array_ts = [y.replace(0, 0.0001) for y in T_hm_array_ts]
figuresHM = []
HM_titles = []
for y in range(len(T_hm_array_ts)):
    T_hm_array_ts[y].index = pd.CategoricalIndex(T_hm_array_ts[y].index, categories=["degradation", "fragmentation", "mixing",
                                                 "biofilm", "heteroagg", "breakup", "advection", "settling", "rising", "resusp", "burial", "sedTransport", "defouling"])
    T_hm_array_ts[y].sort_index(level=0, inplace=True)
    figHM = plt.figure(figsize=(10, 5))
    log_norm = LogNorm(vmin=T_hm_array_ts[y].min(
    ).min(), vmax=T_hm_array_ts[y].max().max())
    ax = sns.heatmap(T_hm_array_ts[y], mask=(T_hm_array_ts[y] == 0.0001), annot=False, annot_kws={
                     "size": 6}, fmt=".2g", cmap='GnBu', linewidth=0.01, linecolor='white', square=True, robust=True, cbar_kws={'label': 't_half (hours)'}, norm=log_norm)
    ax.invert_xaxis()
    ax.set_title(compartments[y])
    # Drawing the frame
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
    ax.set_facecolor('lightgrey')
    plt.show()
    figuresHM.append(figHM)
    figTitle = "halfLife_hm_" + compartments[y] + "_" + composition #!!! 03.14 XZ revised the title
    HM_titles.append(figTitle)



"""SOLVE SYSTEM OF ODES"""

# Initial number of particles in the system
PartNum_t0 = pd.DataFrame(index=Clist, columns=['number of particles'])
for p in range(len(PartNum_t0)):
    PartNum_t0.iloc[p][0] = 0
PartNum_t0.loc["C_"+imputMP] += imputPulse / 60

# Inflow of particles as particles per second
Ilist = []
for C in Clist:
    Ilist.append("I" + C[1:])
inflow_vector = pd.DataFrame(index=Ilist, columns=["number of particles"])
inflow_vector.loc[:, :] = 0
inflow_vector.loc["I_"+imputMP] = imputFlow / 60  # transformed to particles per sec


# intitial conditions
N0 = PartNum_t0['number of particles'].to_numpy(dtype="float")
I = inflow_vector['number of particles'].to_numpy(dtype="float")
E = - ( N0 + I )
# time points
# time = np.linspace(0, tmax, int(timesteps)+1, dtype=int)  # in seconds
time = daysSimulation


# Solve ODEs
if SOLVER == 'SteadyState':
    k = interactions_df.to_numpy()
    Nfinal=np.linalg.solve(k, E)
    Nfinal_T = np.transpose(Nfinal)
    NFinal_num = pd.DataFrame(data=Nfinal_T.reshape(1,-1), index=['number of particles'], columns=Clist)
    

elif SOLVER == "Dynamic":
    print("Solver - Dynamic could be applied in the other .py file")



# Vector of volumes corresponding to the compartments of the river
dilution_vol_m3 = volumesVector(Clist, compartments_prop)
ConcFinal_num_m3 = pd.DataFrame(data=0, index=['number of particles'], columns=Clist)
dilution_vol_array = np.array(dilution_vol_m3)
for ind in range(len(NFinal_num)):
    ConcFinal_num_m3.iloc[ind] = NFinal_num.iloc[ind]/dilution_vol_array[ind] #[0]    
    

# Substitute values smaller than 10-15 to 0
ConcFinal_num_m3 = ConcFinal_num_m3.apply(
    lambda x: [y if y >= 1e-15 else 0 for y in x])

# Concentration in mass per volume
volume = RC_df.loc["volume_m3"].to_numpy()
density = RC_df.loc["density_kg_m3"].to_numpy()
ConcFinal_mg_m3 = ConcFinal_num_m3*volume*density*10**6



"""RESULTS"""

# CREATE RESULTS AND FIGURES & TABLES FOLDERS

# Set current working directory
cwd = os.getcwd()
os.chdir(cwd+"/Results")

# create folder for the day if it doesnt already exists
path = cwd+"/Results"
os.path.isdir(path)
old_path = (daterun_label)
new_path = os.path.isdir(old_path)
if not new_path:
    os.makedirs(old_path)
    print("Created Folder : ", old_path)
else:
    print(old_path, "folder already exists.")

results_path= cwd+"/Results/"+old_path

os.chdir(results_path)

Fig_folder= "/Figures_SS"
os.path.isdir(results_path+Fig_folder)

new_path = os.path.isdir(results_path+Fig_folder)
if not new_path:
    os.mkdir("Figures_SS")
    print("Created Folder : ", Fig_folder)
else:
    print(Fig_folder, "folder already exists.")
    
Tab_folder= "/Tables_SS"
os.path.isdir(results_path+Tab_folder)

new_path = os.path.isdir(results_path+Tab_folder)
if not new_path:
    os.mkdir("Tables_SS")
    print("Created Folder : ", Tab_folder)
else:
    print(Tab_folder, "folder already exists.")

results_figures_path= results_path+Fig_folder
results_tables_path= results_path+Tab_folder


# Function to extract concentration values by size fraction
def extract_SizeBins(comp, MPform): 
    Aa = []
    Ab = []
    Ac = []
    Ad = []
    Ae = []
    for i in range(len(listRS)):
        Aa.append(ConcFinal_num_m3.iloc[0, Clist.index(
            "C_"+str(listRS[i])+comp+MPform+"a")])
        Ab.append(ConcFinal_num_m3.iloc[0, Clist.index(
            "C_"+str(listRS[i])+comp+MPform+"b")])
        Ac.append(ConcFinal_num_m3.iloc[0, Clist.index(
            "C_"+str(listRS[i])+comp+MPform+"c")])
        Ad.append(ConcFinal_num_m3.iloc[0, Clist.index(
            "C_"+str(listRS[i])+comp+MPform+"d")])
        Ae.append(ConcFinal_num_m3.iloc[0, Clist.index(
            "C_"+str(listRS[i])+comp+MPform+"e")])
    return [Aa, Ab, Ac, Ad, Ae]


# Function to extract lists from a list by criteria
def listofindex(criteria, Clist):
    lista = [[] for x in range(len(criteria))]
    for i in range(len(lista)):
        lista[i] = [n for n in Clist if criteria[i] in n[-3:]]
    return lista


# Extract list of indexes needed for plotting
list_of_indexesMpType = listofindex(MPforms, Clist)
list_of_indexesCompartments = listofindex(riverComp, Clist)
list_ofindexesSizeBins = listofindex(sizeBin, Clist)


# Distribution of MPs per aggregation state and compartment in the steady state
MpTypeNum_ss = pd.DataFrame(index=['number of particles'], columns=[m+" (Total number)" for m in MPformslabels]+["Total"])

# Relative abundance of MPs aggregation states in the whole system in the steady state
RelativeAbun_MPtype_ss = pd.DataFrame(0, columns=[m+" (%)" for m in MPformslabels], index=MpTypeNum_ss.index)

compNum_ss = pd.DataFrame(index=['number of particles'], columns=[m+" (Total number)" for m in compartments])

# Fractionation of MPs per compartmet for the whole system in the steady state
RelativeAbun_Comp = pd.DataFrame(0, columns=[m+" (%)" for m in compartments], index=MpTypeNum_ss.index)



# Convert concentration to particle number
PartNum_ss = ConcFinal_num_m3.loc['number of particles']*dilution_vol_m3
MpTypeNum_ss.loc['number of particles', 'Total'] = sum(PartNum_ss)
PartNum_ss = PartNum_ss.to_frame()
for mp in range(1, 1+len(MPforms)):
    MpTypeNum_ss.loc['number of particles', MPformslabels[mp-1] + " (Total number)"] = sum(
        PartNum_ss.loc[list_of_indexesMpType[mp-1], :]['number of particles'].to_list())
    if MpTypeNum_ss.loc['number of particles', 'Total'] == 0:
        RelativeAbun_MPtype_ss.loc['number of particles', MPformslabels[mp-1] + " (%)"] = 0
    else:
        RelativeAbun_MPtype_ss.loc['number of particles', MPformslabels[mp-1] + " (%)"] = round(
            (MpTypeNum_ss.loc['number of particles', MPformslabels[mp-1] + " (Total number)"] / MpTypeNum_ss.loc['number of particles', 'Total']) * 100, 2)
for com in range(1, 1+len(compartments)):
    compNum_ss.loc['number of particles', compartments[com-1] + " (Total number)"] = sum(
        PartNum_ss.loc[list_of_indexesCompartments[com-1], :]['number of particles'].to_list())
    if MpTypeNum_ss.loc['number of particles', 'Total'] == 0:
        RelativeAbun_Comp.loc['number of particles', compartments[com-1] + " (%)"] = 0
    else:
        RelativeAbun_Comp.loc['number of particles', compartments[com-1] + " (%)"] = round(
            (compNum_ss.loc['number of particles', compartments[com-1] + " (Total number)"] / MpTypeNum_ss.loc['number of particles', 'Total']) * 100, 2)



# Estimate total number of particles per size fraction in each compartment for the whole river
SizeFracNum_comp = pd.DataFrame(index=sizeBinLabel, columns=compartments)
for siz in range(len(sizeBin)):
    for co in range(len(SizeFracNum_comp.columns)):
        list_sb = [n for n in list_of_indexesCompartments[co]
                   if sizeBin[siz] in n[-3:]]
        SizeFracNum_comp.iloc[siz, co] = sum(
            NFinal_num[list_sb].loc['number of particles'].to_list())
        
        

Total_num_water = sum(sum(SizeFracNum_comp[compartments[0:3]].values))
Total_num_sediment = sum(SizeFracNum_comp[compartments[-1]].values)
RelativeAbun_SizeFrac = pd.DataFrame(index=sizeBinLabel, columns=[
                                     "Relative abundance in water (%)", "Relative abundance in sediment (%)"])
for siz in range(len(sizeBin)):
    if sum(SizeFracNum_comp.iloc[siz][compartments[0:3]].values) == 0:
        RelativeAbun_SizeFrac.iloc[siz, 0] = 0
    else:
        RelativeAbun_SizeFrac.iloc[siz, 0] = round(
            (sum(SizeFracNum_comp.iloc[siz][compartments[0:3]].values)/Total_num_water)*100, 2)
    if SizeFracNum_comp.iloc[siz][compartments[-1]] == 0:
        RelativeAbun_SizeFrac.iloc[siz, 1] = 0
    else:
        RelativeAbun_SizeFrac.iloc[siz, 1] = round(
            (SizeFracNum_comp.iloc[siz][compartments[-1]]/Total_num_sediment)*100, 2)

# Estimate total number of particles per size fraction in each compartment for the downstream section (last 200 km )
SizeFracNum_comp_downs = pd.DataFrame(index=sizeBinLabel, columns=compartments)
for siz in range(len(sizeBin)):
    for co in range(len(SizeFracNum_comp.columns)):
        list_sb = [n for n in list_of_indexesCompartments[co]
                   if sizeBin[siz] in n[-3:]]
        list_sb_downstream = list_sb[20*4-4*4:]
        SizeFracNum_comp_downs.iloc[siz, co] = sum(
            NFinal_num[list_sb_downstream].loc['number of particles'].to_list())
        

SizeFracNum_comp_downs
Total_num_water_downs = sum(
    sum(SizeFracNum_comp_downs[compartments[0:3]].values))
Total_num_sediment_downs = sum(SizeFracNum_comp_downs[compartments[-1]].values)
RelativeAbun_SizeFrac_downs = pd.DataFrame(index=sizeBinLabel, columns=[
                                           "Relative abundance in water (%)", "Relative abundance in sediment (%)"])
for siz in range(len(sizeBin)):
    if sum(SizeFracNum_comp_downs.iloc[siz][compartments[0:3]].values) == 0:
        RelativeAbun_SizeFrac_downs.iloc[siz, 0] = 0
    else:
        RelativeAbun_SizeFrac_downs.iloc[siz, 0] = round(
            (sum(SizeFracNum_comp_downs.iloc[siz][compartments[0:3]].values)/Total_num_water_downs)*100, 2)
    if SizeFracNum_comp_downs.iloc[siz][compartments[-1]] == 0:
        RelativeAbun_SizeFrac_downs.iloc[siz, 1] = 0
    else:
        RelativeAbun_SizeFrac_downs.iloc[siz, 1] = round(
            (SizeFracNum_comp_downs.iloc[siz][compartments[-1]]/Total_num_sediment_downs)*100, 2)



"""PLOTS"""

""" Multyplots graphs: Concentration vs distance over time """

# Set x values (distance in km)
# Distance values
x = [d/1000 for d in RS_cumLength_m]
compartmentsLabel = ["Surface\n Water",
                     "Flowing\n Water", 
                     "Stagnant\n Water", 
                     "Sediment"]

# Choose style and colour palette
palette = plt.get_cmap('Set2')
plt.style.use('seaborn-white')
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='ivory', alpha=0.5)
props2 = dict(boxstyle='round', facecolor='white', alpha=0.5)

# Move to figures folder
os.chdir(results_figures_path)
os.chdir(cwd)

### Static plot (t=tfinal) ###

fig2, axs = plt.subplots(len(compartments), len(MPforms), figsize=(
    15, 10), sharex='col', sharey="row", squeeze="True")

labels = ['0.1 um', '1 um', '10 um', '100 um', '1000 um']
if imputFlow == 0:
    fig2.suptitle(composition + " plastic particles fate along the River " + imputRiver + " (pulse= " +
                  str(imputPulse)+" particles)", fontsize=20)
else:
    fig2.suptitle(composition + " plastic particles after " +
                  str(int(daysSimulation)) + " days", fontsize=18,  y=0.95)



qnum = ['number of particles']
for j in range(len(compartments)):
    if j == 3:
        for k in range(len(MPforms)):
            # Plot
            y = extract_SizeBins(riverComp[j], MPforms[k])
            axs[j, k].plot(x, [e / (10**6*sed_density) for e in y[0]],
                           linewidth=2.5, color=palette(0), label='0.1 um')
            axs[j, k].plot(x, [e / (10**6*sed_density) for e in y[1]],
                           linewidth=2.5, color=palette(1), label='1 um')
            axs[j, k].plot(x, [e / (10**6*sed_density) for e in y[2]],
                           linewidth=2.5, color=palette(2), label='10 um')
            axs[j, k].plot(x, [e / (10**6*sed_density) for e in y[3]],
                           linewidth=2.5, color=palette(3), label='100 um')
            axs[j, k].plot(x, [e / (10**6*sed_density) for e in y[4]],
                           linewidth=2.5, color=palette(4), label='1000 um')

            if k == 3:
                axs[j, k].text(1.2, 0.5,  compartmentsLabel[j]+"\n"+str(RelativeAbun_Comp.loc['number of particles', compartments[j] + " (%)"]) +" %",
                               fontsize=15, rotation=0, va='center', ha='center', bbox=props2, transform=axs[j, k].transAxes)

            axs[j, k].set_yscale('log')
            axs[j, k].set_ylim(10**-9, 10**3)
            if k == 0:
                axs[j, k].set_ylabel("Conc (Num/g)", fontsize=15)
            axs[j, k].yaxis.set_major_locator(
                ticker.LogLocator(base=10.0, numticks=4))
            axs[j, k].set_xlim(x[0], x[-1])
            axs[j, k].tick_params(
                axis='x', labelsize=12, direction='inout', length=6, width=1, grid_alpha=0.5)
            axs[j, k].tick_params(
                axis='y', labelsize=10, direction='inout', length=6, width=1, grid_alpha=0.5)
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
    else:
        for k in range(len(MPforms)):
            # Plot
            y = extract_SizeBins(riverComp[j], MPforms[k])
            axs[j, k].plot(x, y[0], linewidth=2.5,
                           color=palette(0), label='0.1 um')
            axs[j, k].plot(x, y[1], linewidth=2.5,
                           color=palette(1), label='1 um')
            axs[j, k].plot(x, y[2], linewidth=2.5,
                           color=palette(2), label='10 um')
            axs[j, k].plot(x, y[3], linewidth=2.5,
                           color=palette(3), label='100 um')
            axs[j, k].plot(x, y[4], linewidth=2.5,
                           color=palette(4), label='1000 um')
            if j == 0: #!!! 03.15 XZ modified the figure text label
                #axs[j, k].text(0.5, 1.1, MPformslabels[k], fontsize=15,
                #               transform=axs[j, k].transAxes, ha='center')
                axs[j, k].text(0.5, 1.1, MPformslabels[k] + "\n " + str(RelativeAbun_MPtype_ss.loc['number of particles', MPformslabels[k] + " (%)"]
                                                                        ) + " %", fontsize=15, bbox=props2, transform=axs[j, k].transAxes, ha='center')            
            if k == 3:
                axs[j, k].text(1.2, 0.5, compartmentsLabel[j]+"\n"+str(RelativeAbun_Comp.loc['number of particles', compartments[j] + " (%)"]) +" %",
                               fontsize=15, rotation=0, va='center', ha='center', bbox=props2, transform=axs[j, k].transAxes)
            if k == 0:
                axs[j, k].set_ylabel("Conc (Num/$m^3$)", fontsize=15)
            axs[j, k].set_yscale('log')

            if j == 0:
                axs[j, k].set_ylim(10**-6, 10**6)
                axs[j, k].yaxis.set_major_locator(
                    ticker.LogLocator(base=10.0, numticks=4))
            elif j == 1:
                axs[j, k].set_ylim(10**-6, 10**6)
                axs[j, k].yaxis.set_major_locator(
                    ticker.LogLocator(base=10.0, numticks=4))
            elif j == 2:
                axs[j, k].set_ylim(10**-6, 10**6)
                axs[j, k].yaxis.set_major_locator(
                    ticker.LogLocator(base=10.0, numticks=4))
            axs[j, k].set_xlim(x[0], x[-1])

            axs[j, k].tick_params(
                axis='x', labelsize=10, direction='inout', length=6, width=1, grid_alpha=0.5)
            axs[j, k].tick_params(
                axis='y', labelsize=10, direction='inout', length=6, width=1, grid_alpha=0.5)
            from matplotlib import ticker
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            axs[j, k].minorticks_on()


# Axis titles
#plt.text(0.02, 0.5, "Concentration of particles (Num/$m^3$)", fontsize=15, transform=plt.gcf().transFigure, rotation='vertical',ha='center', va='center')
plt.text(0.5, 0.07, "Distance (km)", fontsize=15,
         transform=plt.gcf().transFigure, ha='center', va='center')
#plt.legend(labels,bbox_to_anchor=(0.5, -0.18), loc='center',ncol=5, fontsize=15 )
plt.subplots_adjust(wspace=0.02, hspace=0.1)
handles, labels = axs[j, k].get_legend_handles_labels()
fig2.legend(handles, labels, bbox_to_anchor=(
    0.5, 0.03), loc='center', ncol=5, fontsize=12)
fig2_label = "ConcvsDist_Num_m3_Multiplot_" + composition + ".png"



"""SAVE RESULTS AND FIGURES"""

# Set current working directory
cwd = os.getcwd()

if record == "True":
    # create folder for the day if it doesnt already exists
    os.chdir(cwd+"/Results")

    # create folder for the day if it doesnt already exists
    path = cwd+"/Results"
    os.path.isdir(path)
    old_path = (daterun_label)
    new_path = os.path.isdir(old_path)
    if not new_path:
        os.makedirs(old_path)
        print("Created Folder : ", old_path)
    else:
        print(old_path, "folder already exists.")

    results_path = cwd+"/Results/"+old_path

    os.chdir(results_path)

    Fig_folder = "/Figures_SS"
    os.path.isdir(results_path+Fig_folder)

    new_path = os.path.isdir(results_path+Fig_folder)
    if not new_path:
        os.mkdir("Figures_SS")
        print("Created Folder : ", Fig_folder)
    else:
        print(Fig_folder, "folder already exists.")

    results_figures_path = results_path+Fig_folder

    # save
    # rate constants
    dfRC_filename = "RC_df" + runtitle + "_" + daterun_label + ".csv"
    RC_df_final.to_csv(dfRC_filename)

    # interactions dataframe to results folder
    interactions_filename = "interactionsdf_" + \
        runtitle + "_" + daterun_label + ".csv"
    interactions_df.to_csv(interactions_filename)

    # Results
    os.chdir(cwd)
    os.chdir(results_path)
    filename = "ConcFinal_Num_m3_" + runtitle + "_" + daterun_label+".csv"
    filename1 = "MassVStime_kg_" + runtitle + "_" + daterun_label+".csv"
    ConcFinal_num_m3.to_csv(filename)
    #ConcFinal_mg_m3.to_csv(filename1)
    RelativeAbun_MPtype_ss.to_csv(
        "Relative_abundance_MPtype_"+composition+daterun_label+".csv")
    RelativeAbun_Comp.to_csv(
        "Relative_abundance_Compartment_"+composition+daterun_label+".csv")

    # back to working directory
    os.chdir(cwd)

    # Save figures
    os.chdir(results_figures_path)
    fig2.savefig(fig2_label)
    for m in range(len(figuresHM)):
        hm = figuresHM[m]
        hm_label = HM_titles[m]
        hm.savefig(hm_label)
