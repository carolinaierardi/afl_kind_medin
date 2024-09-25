#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:23:00 2024

@author: carolinaierardi
"""

import os
import numpy as np 
from matplotlib import pyplot as plt
import pandas as pd

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data")


#%%Variables to change
Abeta_af = pd.read_csv("affinity.csv")            # affinities
inputs_order = pd.read_csv("inference_Abeta.csv") # order of molecules that were read
input_mol = pd.read_csv("input_Ab.csv")           # names of molecules with IDS

directory_mol_feat = "output_Ab_mol"

#%% Reorder and get corresponding smiles and affinities

#Obtain order of IDs
mols = list(inputs_order["mol_feat_file"])

#Define words to split string in
split_words = [directory_mol_feat + "/[","]"]

#split to obtain just the number and make into integer
split_mols = [i.split(split_words[0], 1)[1] for i in mols]
split_mols = [i.split(split_words[1],1)[0] for i in split_mols]
finalindex = [int(i) for i in split_mols]

#order the names of the molecules in the library original df
ordered_mols = list(input_mol["Compound Name"][finalindex])
ordered_smiles = list(input_mol["SMILES"][finalindex])

#obtain the names of the molecules and their corresponding affinity
affinity_names = pd.DataFrame({"Compound Name" : ordered_mols,
                               "SMILES" : ordered_smiles,
                               "Affinity" : Abeta_af["affinity"]})

affinity_names = affinity_names.sort_values("Affinity", ascending = False)
true_af = pd.read_csv("cayman_processed_cns.csv")


#%% Make figures 

hfont = {'fontname':'Arial'}  

fig, axs = plt.subplot_mosaic("AB",figsize=(10,4))          
                                                       #change font to arial  
axs["A"].hist(affinity_names["Affinity"], color = "#48D1CC", bins = 20, density=True, edgecolor='black', linewidth=1)
axs["A"].set_title(r" New model A$\beta$42 affinities", fontsize = 15) 
axs["A"].set_xlabel("Binding Affinity",**hfont, fontsize = 12);
axs["A"].set_ylabel("Frequency",**hfont, fontsize = 12);

axs["B"].hist(true_af["predicted affinity(ModelC)"], color = "#000080", bins = 20, density=True, edgecolor='black', linewidth=1)
axs["B"].set_title(r"Old model C A$\beta$42 affinities", fontsize = 15) 
axs["B"].set_xlabel("Binding Affinity",**hfont, fontsize = 12);
axs["B"].set_ylabel("Frequency",**hfont, fontsize = 12);
plt.tight_layout(h_pad = 2)  


sorted_new = affinity_names.sort_values(by=["Compound Name"], axis = 0)
sorted_old = true_af.sort_values(by=["cmpdname"], axis = 0)

b, a = np.polyfit(np.array(sorted_new["Affinity"]), 
            np.array(sorted_old["predicted affinity(ModelC)"]), deg=1)

pearsonr = np.corrcoef(np.array(sorted_new["Affinity"]), 
            np.array(sorted_old["predicted affinity(ModelC)"]))[0,1]

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(np.array(sorted_new["Affinity"]), 
            np.array(sorted_old["predicted affinity(ModelC)"]),
            c = "#48D1CC", ec = "black") 
ax.plot(sorted_new["Affinity"], a + b*np.array(sorted_new["Affinity"]), 
    color="k", lw=2.5)
ax.set_title(f"Pearson's R = {round(pearsonr,2)}", fontsize = 15) 
ax.set_xlabel("New model",**hfont, fontsize = 12);
ax.set_ylabel("Old model",**hfont, fontsize = 12);



#%% True affinities

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    
    b, a = np.polyfit(np.array(x), np.array(y), deg=1)

    pearsonr = np.corrcoef(np.array(x),np.array(y))[0,1]

    # the scatter plot:
    ax.scatter(x, y, c = "#48D1CC", ec = "black")
    ax.plot(np.array(x), a + b*np.array(x), color="k", lw=2.5)
    ax.set_xlabel("New model",**hfont, fontsize = 12);
    ax.set_ylabel("Old model",**hfont, fontsize = 12);

    # now determine nice limits by hand:
    binwidth = 0.25
    upp_lim = max(np.max(x), np.max(y))
    inf_lim = min(np.min(x), np.min(y))

    bins = np.arange(inf_lim, upp_lim, binwidth)
    ax_histx.hist(x, color = "#000080", bins = bins, density=True, 
                  edgecolor='black', linewidth=1)
    ax_histy.hist(y, orientation='horizontal',color = "#000080", bins = bins, 
                  density=True, edgecolor='black', linewidth=1)


#set(true_af["cmpdname"]), set(affinity_names["Compound Name"])

fig = plt.figure(figsize=(6, 6))
gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

scatter_hist(np.array(sorted_new["Affinity"]), 
            np.array(sorted_old["predicted affinity(ModelC)"]),
                     ax, ax_histx, ax_histy)



#%% Kinetics data


#Binding DB has only 240 values with koff and 30 with kon and koff
#https://www.bindingdb.org/rwd/bind/ByKI.jsp?specified=Kn

#this database has 3812 structures
#a lot of the targets are kinases

#https://pubs.rsc.org/en/content/articlelanding/2020/MD/d0md00178c#!divAbstract

#perhaps more available data
#https://kbbox.h-its.org/toolbox/data/

#https://www.csbj.org/action/showPdf?pii=S2001-0370%2824%2900068-0
#http://www.pdbbind.org.cn/index.php
#https://www.researchgate.net/publication/336555671_KOFFI_and_Anabel_20-a_new_binding_kinetics_database_and_its_integration_in_an_open-source_binding_analysis_software


