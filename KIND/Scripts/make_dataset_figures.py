#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:03:11 2024

@author: carolinaierardi
"""

import os 
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
    

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")


with open("preprocessed_kinetics_data.pkl", "rb") as f:
    input_data = pickle.load(f)
    
with open("data_splits.pkl", "rb") as f:
     data_splits = pickle.load(f)
     
plt.rcParams['font.family'] = 'Arial'
    
#%% Functions needed

def get_Morgan_fingerprint(smiles, rad = 2, bits = 124):
    
    """
    Parameters
    ----------
    smiles : str
        Smiles for a molecule.
    rad : int, optional
        Radius to use as parameters. The default is 2.
    bits : int, optional
        Bits to use as parameter. The default is 124.

    Returns
    -------
    vec : array
        array with Morgan Fingerprint for molecule.

    """
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, 
                                                radius=rad, nBits=bits)
    
    vec = np.array(fp)

    return vec


# def assign_p_symbol(p):
    
#     if 5.00e-02 < p <= 1.00e+00:
#         return "ns"
#     elif 1.00e-02 < p <= 5.00e-02:
#         return "*"
#     elif 1.00e-03 < p <= 1.00e-02:
#         return "**"
#     elif 1.00e-04 < p <= 1.00e-03:
#         return "***"
#     elif p <= 1.00e-04:
#         return "****"

# def add_stat_annotation(ax, data, x, y, box_pairs):
#     i = 0
#     for pair in box_pairs:
#         group1 = data[data[x] == pair[0]][y]
#         group2 = data[data[x] == pair[1]][y]
#         stat, p = ttest_ind(group1, group2)
#         y_max = max(group1.max(), group2.max())
#         ax.plot([0, 1], [y_max, y_max], color='black')
#         ax.text(0.5 + i*0.5, y_max, assign_p_symbol(p), ha='center')
#         i += 1


#%% Make plots

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data_preprocess/Figures")

#%% Distribution of values

n_prot = input_data["Protein_ID"].unique()
n_ligands = input_data["Ligand_ID"].unique()

print(f"Information provided: {list(input_data)}")
print(f"There are {len(n_prot)} unique proteins accounted for")
print(f"There are {len(input_data['target_class'].unique())} protein types")
print(f"There are {len(n_ligands)} unique ligands accounted for")


fig, axs = plt.subplot_mosaic("ABC",figsize=(12,4))                                   #get mosaic plot 

hfont = {'fontname':'Arial'}                                                         #change font to arial  
axs['A'].hist(input_data["pKD"], color="#479ad1",ec = "black", bins = 30)                            #first plot will be of mean FD
axs['A'].set_xlabel("pKD" ,**hfont, fontsize = 16)           #x-label
axs['A'].set_ylabel("Frequency",**hfont, fontsize = 16)                       #y-label
axs['A'].set_title("pKD",**hfont, weight = 'bold', fontsize = 16) #title for subplot
axs['A'].axvline(x = np.mean(input_data["pKD"]), color = 'k', linestyle = "dashed") #add dashed liine for excluded participants
axs['A'].tick_params(axis = 'both', labelsize=14)                                    #adjust size of the axis ticks
axs['A'].text(-0.1, 1.1, 'A', transform=axs['A'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot

                                                     
axs['B'].hist(input_data["pkon"], bins = 30, color="#479ad1", ec = "black")                   #second subplot will be of mean connectivity
axs['B'].set_xlabel("pKon",**hfont, fontsize = 16)             #x-label
axs['B'].set_ylabel("Frequency",**hfont, fontsize = 16)                         #y-label
axs['B'].set_title("pKon",**hfont, weight = 'bold', fontsize = 16)                                   #title for this subplot
axs['B'].axvline(x = np.mean(input_data["pkon"]) ,color = 'k', linestyle = "dashed")                                               #dashed line for excluded participants
axs['B'].tick_params(axis = 'both', labelsize=14)                                    #adjust size of the axis ticks
axs['B'].text(-0.1, 1.1, 'B', transform=axs['B'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot


axs['C'].hist(input_data["pkoff"], bins = 30, color="#479ad1", ec = "black")                   #second subplot will be of mean connectivity
axs['C'].set_xlabel("pKoff",**hfont, fontsize = 16)             #x-label
axs['C'].set_ylabel("Frequency",**hfont, fontsize = 16)                         #y-label
axs['C'].set_title("pKoff",**hfont, weight = 'bold', fontsize = 16)                                   #title for this subplot
axs['C'].axvline(x = np.mean(input_data["pkoff"]) ,color = 'k', linestyle = "dashed")                                               #dashed line for excluded participants
axs['C'].tick_params(axis = 'both', labelsize=14)                                    #adjust size of the axis ticks
axs['C'].text(-0.1, 1.1, 'C', transform=axs['C'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots


plt.savefig('Kinetics_distribution.png')

#%% Barplot per type of protein and assay

orders_type = list(input_data["target_class"].value_counts().index)
orders_assay = list(input_data["assay_method"].value_counts().index)

assays = ["FLBA","SPR","RL","Unknown","xCELLigence"]
target_types = ["Kinase", "GPCR", "HSP","Enzyme","Unknown","VGIC"]

fig, (ax1, ax2) = plt.subplots(2,1, figsize = (12,14))                                   #get mosaic plot 

ax1.bar(target_types,
        input_data["target_class"].value_counts(), 
        color = sns.color_palette("crest")[:5], ec = "k")
ax1.tick_params(axis = 'both', labelsize=20)   
ax1.set_title("Target class", **hfont, size = 20)


ax2.bar(assays,
        input_data["assay_method"].value_counts(), 
        color = sns.color_palette("magma")[:5], ec = "k")
ax2.tick_params(axis = 'both', labelsize=20) 
ax2.set_title("Assay method", **hfont, size = 20)

plt.savefig('Demographic_bars.png')


#%% Boxplots between types 

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 14))
fig.subplots_adjust(hspace=0.2)

# Boxplot for pKD by target_class
ax = sns.violinplot(x='target_class', y='pKD', order = orders_type,
                    data=input_data, ax=axes[0, 0],
                 palette='crest')
ax.set_xticklabels(target_types,**hfont, size = 12)
axes[0, 0].set_title('pKD',**hfont, size = 16)

# Boxplot for pKon by target_class
ax = sns.violinplot(x='target_class', y='pkon', order = orders_type,
                    data=input_data, ax=axes[0, 1], palette='crest')
# add_stat_annotation(ax, data=input_data, x='target_class', y='pkon',
#                     box_pairs=list(combinations(input_data["target_class"].unique(),2)))
ax.set_xticklabels(target_types,**hfont, size = 12)

axes[0, 1].set_title('pKon',**hfont, size = 16)

# Boxplot for pKoff by target_class
ax = sns.violinplot(x='target_class', y='pkoff', order = orders_type,
                    data=input_data, ax=axes[0, 2], palette='crest')
ax.set_xticklabels(target_types,**hfont, size = 12)
axes[0, 2].set_title('pKoff',**hfont, size = 16)


# Boxplot for pKD by assay_method
ax = sns.violinplot(x='assay_method', y='pKD', order = orders_assay,
                    data=input_data, ax=axes[1, 0], palette='magma')
ax.set_xticklabels(assays,**hfont, size = 12)
axes[1, 0].set_title('pKD',**hfont, size = 16)


# Boxplot for pKon by assay_method
ax = sns.violinplot(x='assay_method', y='pkon',  order = orders_assay,
                    data=input_data, ax=axes[1, 1], palette='magma')
ax.set_xticklabels(assays, **hfont, size = 12)
axes[1, 1].set_title('pKon', **hfont, size = 16)


# Boxplot for pKoff by assay_method
ax = sns.violinplot(x='assay_method', y='pkoff',  order = orders_assay,
                    data=input_data, ax=axes[1, 2], palette='magma')
ax.set_xticklabels(assays,**hfont, size = 12)
axes[1, 2].set_title('pKoff',**hfont, size = 16)


plt.savefig('Demographic_violin.png')


#%% Correlation between kinetics values


ax = sns.jointplot(data=input_data, x="pKD", y="pkoff", kind = "reg", 
              color = "#48D1CC", joint_kws = {'color':"#339b97"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pKD'], input_data['pkoff'])[1,0],3)}",
                **hfont, size=16)
ax.set_axis_labels('pKD', 'pKoff', **hfont, fontsize = 14)

plt.savefig('pKD_pKoff.png')


ax = sns.jointplot(data=input_data, x="pKD", y="pkon", kind = "reg", 
              color = "#d14747", joint_kws = {'color':"#9b3232"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pKD'], input_data['pkon'])[1,0],3)}",
                **hfont, size=16)
ax.set_axis_labels('pKD', 'pKon', **hfont, fontsize = 14)

plt.savefig('pKD_pKon.png')


ax = sns.jointplot(data=input_data, x="pkon", y="pkoff", kind = "reg", 
              color = "#6747d1", joint_kws = {'color':"#4b329b"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pkon'], input_data['pkoff'])[1,0],3)}",
                **hfont, size=16)
ax.set_axis_labels('pKon', 'pKoff', **hfont, fontsize = 14)

plt.savefig('pKon_pKoff.png')

#%% Dimensionality reduction for molecules

mols = input_data["canonical_smiles"].values

unique_mol_test = np.unique(data_splits["testing"]["canonical_smiles"].values)

fingerprints = [get_Morgan_fingerprint(sm, bits = 124) for sm in mols]
fingerprints = np.array(fingerprints)


# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(fingerprints)

# t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(fingerprints)

# UMAP
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_result = umap_model.fit_transform(fingerprints)

#To plot the train and test set in different colours
plotting_groups = np.zeros(len(input_data))
plotting_groups[np.array(data_splits["testing"].index)] = 1
plotting_groups = np.array([int(i) for i in plotting_groups])

colors = np.array(["lightblue", "green"])

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15,5))

axes[0].scatter(pca_result[:, 0], pca_result[:, 1], s=15, c=colors[plotting_groups])
axes[0].set_title("PCA", size = 20)
axes[0].set_xlabel('Component 1', size = 15)
axes[0].set_ylabel('Component 2', size = 15)

axes[1].scatter(tsne_result[:, 0], tsne_result[:, 1], s=15, c=colors[plotting_groups])
axes[1].set_title("t-SNE", **hfont, size = 20)
axes[1].set_xlabel('Component 1',  size = 15)
axes[1].set_ylabel('Component 2', size = 15)

axes[2].scatter(umap_result[:, 0], umap_result[:, 1], s=15, c=colors[plotting_groups])
axes[2].set_title("UMAP", **hfont, size = 20)
axes[2].set_xlabel('Component 1',  size = 15)
axes[2].set_ylabel('Component 2', size = 15)

fig.suptitle("Morgan Fingerprint (Train X Test)", **hfont, size = 25)

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots

plt.savefig('train_test_dimred.png')


#For the training validation splits

fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25,5))
colors = np.array(["lightblue", "green", "darkred"])


for i in range(len(data_splits["subtraining"])):
    
    plotting_groups = np.zeros(len(input_data))
    plotting_groups[np.array(data_splits["testing"].index)] = 1
    plotting_groups[np.array(data_splits["validation"][i].index)] = 2
    plotting_groups = np.array([int(i) for i in plotting_groups])


    axes[i].scatter(pca_result[:, 0], pca_result[:, 1], 
                    s=20, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1', **hfont, size = 15)
    axes[i].set_ylabel('Component 2', **hfont, size = 15)
    

fig.suptitle("PCA", **hfont, size = 25)

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_pca.png')


fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25,5))


for i in range(len(data_splits["subtraining"])):
    
    plotting_groups = np.zeros(len(input_data))
    plotting_groups[np.array(data_splits["testing"].index)] = 1
    plotting_groups[np.array(data_splits["validation"][i].index)] = 2
    plotting_groups = np.array([int(i) for i in plotting_groups])


    axes[i].scatter(tsne_result[:, 0], tsne_result[:, 1], s=15, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1', **hfont, size = 15)
    axes[i].set_ylabel('Component 2', **hfont, size = 15)


fig.suptitle("t-SNE", **hfont, size = 25)
fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_tsne.png')


fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25,5))

for i in range(len(data_splits["subtraining"])):
    
    plotting_groups = np.zeros(len(input_data))
    plotting_groups[np.array(data_splits["testing"].index)] = 1
    plotting_groups[np.array(data_splits["validation"][i].index)] = 2
    plotting_groups = np.array([int(i) for i in plotting_groups])


    axes[i].scatter(umap_result[:, 0], umap_result[:, 1], s=15, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1', **hfont, size = 15)
    axes[i].set_ylabel('Component 2', **hfont, size = 15)

fig.suptitle("UMAP", **hfont, size = 25)
fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_umap.png')








