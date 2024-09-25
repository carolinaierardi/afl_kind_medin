#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 12:12:20 2024

@author: carolinaierardi
"""

#This script will reproduce all the figures in the dissertation


#%% Import libraries

import os 
import pandas as pd
import numpy as np

#splits
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold

#rdkit
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


#plotting
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
    

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 20

    

#%% Import data
os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")

with open("preprocessed_kinetics_data.pkl", "rb") as f:
    meta = pickle.load(f)

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data")
input_data = pd.read_csv("kind.csv")


#%% sequences with more than 800 residues

def filter_csv(df):
    filtered_df = df[df['n_residue'] <= 800]
    return filtered_df

meta = filter_csv(meta)


#%% Split into train and test

stratify_column = meta["target_class"] + "_" + meta["assay_method"]

#adjust test size to have training and validation sets of the same size
test_set_size = int(np.round(np.floor(len(stratify_column)*0.1)))

training, testing = train_test_split(meta, 
                                     test_size=test_set_size-1, 
                                     shuffle = True,
                                     random_state = 42, 
                                     stratify=stratify_column) 


#%% Split between training and validation

n_folds = 5
k_split_stratify_column = training["target_class"] + "_" + training["assay_method"]

skfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42) #create stratified K-fold
subtraining = []    #empty list to store training indices for each loop
validating = []     #empty list to store validating indices for each loop

for train_index, test_index in skfold.split(training, k_split_stratify_column): #get the splits for our input data
    subtraining += [training.iloc[train_index]]                          #get the indices for the training set this fold  
    validating += [training.iloc[test_index]]                            #get the indices for the validating set this fold                                     



data_splits = {"training": training,"testing":testing,
               "subtraining":subtraining,"validation":validating}


with open("data_splits.pkl", "wb") as f:
     pickle.dump(data_splits, f)
     
#%% Make figures

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Figures")

fig, axs = plt.subplot_mosaic("ABC",figsize=(12,4))                                   #get mosaic plot 

hfont = {'fontname':'Arial'}                                                         #change font to arial  
axs['A'].hist(meta["pKD"], color="#479ad1",ec = "black", bins = 30)                            #first plot will be of mean FD
axs['A'].set_xlabel("pKD" )           #x-label
axs['A'].set_ylabel("Frequency")                       #y-label
axs['A'].set_title("pKD",weight = 'bold') #title for subplot
axs['A'].axvline(x = np.mean(meta["pKD"]), color = 'k', linestyle = "dashed") #add dashed liine for excluded participants
axs['A'].tick_params(axis = 'both')                                    #adjust size of the axis ticks
axs['A'].text(-0.1, 1.1, 'A', transform=axs['A'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot

                                                     
axs['B'].hist(meta["pkon"], bins = 30, color="#479ad1", ec = "black")                   #second subplot will be of mean connectivity
axs['B'].set_xlabel("pKon")             #x-label
axs['B'].set_ylabel("Frequency")                         #y-label
axs['B'].set_title("pKon", weight = 'bold')                                   #title for this subplot
axs['B'].axvline(x = np.mean(meta["pkon"]) ,color = 'k', linestyle = "dashed")                                               #dashed line for excluded participants
axs['B'].tick_params(axis = 'both')                                    #adjust size of the axis ticks
axs['B'].text(-0.1, 1.1, 'B', transform=axs['B'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot


axs['C'].hist(meta["pkoff"], bins = 30, color="#479ad1", ec = "black")                   #second subplot will be of mean connectivity
axs['C'].set_xlabel("pKoff")             #x-label
axs['C'].set_ylabel("Frequency")                         #y-label
axs['C'].set_title("pKoff", weight = 'bold')                                   #title for this subplot
axs['C'].axvline(x = np.mean(meta["pkoff"]) ,color = 'k', linestyle = "dashed")                                               #dashed line for excluded participants
axs['C'].tick_params(axis = 'both', labelsize=14)                                    #adjust size of the axis ticks
axs['C'].text(-0.1, 1.1, 'C', transform=axs['C'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots


plt.savefig('Kinetics_distribution.png')

#%% Barplot per type of protein and assay

orders_type = list(meta["target_class"].value_counts().index)
orders_assay = list(meta["assay_method"].value_counts().index)

assays = ["FLBA","SPR","RL","Unknown"]
target_types = ["Kinase", "GPCR", "HSP","Enz.","Unknown"]

fig, (ax1, ax2) = plt.subplots(2,1, figsize = (10,16))                                   #get mosaic plot 

ax1.bar(target_types,
        meta["target_class"].value_counts(), 
        color = sns.color_palette("crest")[:4], ec = "k")
ax1.tick_params(axis = 'both', labelsize = 25)   
ax1.set_title("Target class")


ax2.bar(assays,
        meta["assay_method"].value_counts(), 
        color = sns.color_palette("magma")[:5], ec = "k")
ax2.tick_params(axis = 'both', labelsize = 25) 
ax2.set_title("Assay method")

plt.savefig('Demographic_bars.png')

#%% Boxplots between types 

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 18))
fig.subplots_adjust(hspace=0.2)

# Boxplot for pKD by target_class
ax = sns.violinplot(x='target_class', y='pKD', order = orders_type,
                    data=meta, ax=axes[0, 0],
                 palette='crest')
ax.set_xticklabels(target_types, size = 20, rotation = 30)
axes[0, 0].set_title('pKD')


# Boxplot for pKon by target_class
ax = sns.violinplot(x='target_class', y='pkon', order = orders_type,
                    data=meta, ax=axes[0, 1], palette='crest')
# add_stat_annotation(ax, data=input_data, x='target_class', y='pkon',
#                     box_pairs=list(combinations(input_data["target_class"].unique(),2)))
ax.set_xticklabels(target_types, size = 20, rotation = 30)
axes[0, 1].set_title('pKon')

# Boxplot for pKoff by target_class
ax = sns.violinplot(x='target_class', y='pkoff', order = orders_type,
                    data=meta, ax=axes[0, 2], palette='crest')
ax.set_xticklabels(target_types, size = 20, rotation = 30)
axes[0, 2].set_title('pKoff')


# Boxplot for pKD by assay_method
ax = sns.violinplot(x='assay_method', y='pKD', order = orders_assay,
                    data=meta, ax=axes[1, 0], palette='magma')
ax.set_xticklabels(assays, size = 20, rotation = 30)
axes[1, 0].set_title('pKD')


# Boxplot for pKon by assay_method
ax = sns.violinplot(x='assay_method', y='pkon',  order = orders_assay,
                    data=meta, ax=axes[1, 1], palette='magma')
ax.set_xticklabels(assays, size = 20, rotation = 30)
axes[1, 1].set_title('pKon')


# Boxplot for pKoff by assay_method
ax = sns.violinplot(x='assay_method', y='pkoff',  order = orders_assay,
                    data=meta, ax=axes[1, 2], palette='magma')
ax.set_xticklabels(assays, size = 20, rotation = 30)
axes[1, 2].set_title('pKoff')


plt.savefig('Demographic_violin.png')

#%% Correlation between kinetics values


ax = sns.jointplot(data=input_data, x="pKD", y="pkoff", kind = "reg", 
              color = "#48D1CC", joint_kws = {'color':"#339b97"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pKD'], input_data['pkoff'])[1,0],3)}",
                **hfont, size=25)
ax.set_axis_labels('pKD', 'pKoff', **hfont, fontsize = 25)

ax.text(-0.1, 1.1, 'B', transform=axs['B'].transAxes, 
            size=20, weight='bold')                                                  #add the letter at the corner of the plot


plt.savefig('pKD_pKoff.png')


ax = sns.jointplot(data=input_data, x="pKD", y="pkon", kind = "reg", 
              color = "#d14747", joint_kws = {'color':"#9b3232"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pKD'], input_data['pkon'])[1,0],3)}",
                **hfont, size=25)
ax.set_axis_labels('pKD', 'pKon', **hfont, fontsize = 25)

plt.savefig('pKD_pKon.png')


ax = sns.jointplot(data=input_data, x="pkon", y="pkoff", kind = "reg", 
              color = "#6747d1", joint_kws = {'color':"#4b329b"})
ax.fig.suptitle(f"r = {round(np.corrcoef(input_data['pkon'], input_data['pkoff'])[1,0],3)}", size=25)
ax.set_axis_labels('pKon', 'pKoff', fontsize = 25)

plt.savefig('pKon_pKoff.png')


#%% Dimensionality reduction for molecules

def get_Morgan_fingerprint(smiles, rad = 2, bits = 124):
    
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, 
                                                radius=rad, nBits=bits)
    
    vec = np.array(fp)

    return vec


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

axes[0].scatter(pca_result[:, 0], pca_result[:, 1], 
                s=30, c=colors[plotting_groups], alpha = 0.7)
axes[0].set_title("PCA")
axes[0].set_xlabel('Component 1')
axes[0].set_ylabel('Component 2')

axes[1].scatter(tsne_result[:, 0], tsne_result[:, 1], 
                s=30, c=colors[plotting_groups], alpha = 0.7)
axes[1].set_title("t-SNE")
axes[1].set_xlabel('Component 1')
axes[1].set_ylabel('Component 2')

axes[2].scatter(umap_result[:, 0], umap_result[:, 1], 
                s=30, c=colors[plotting_groups], alpha = 0.7)
axes[2].set_title("UMAP")
axes[2].set_xlabel('Component 1')
axes[2].set_ylabel('Component 2')

fig.suptitle("Morgan Fingerprint (Train X Test)", size = 30)

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
                    s=30, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1')
    axes[i].set_ylabel('Component 2')
    

fig.suptitle("PCA",size = 30)

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_pca.png')


fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25,5))


for i in range(len(data_splits["subtraining"])):
    
    plotting_groups = np.zeros(len(input_data))
    plotting_groups[np.array(data_splits["testing"].index)] = 1
    plotting_groups[np.array(data_splits["validation"][i].index)] = 2
    plotting_groups = np.array([int(i) for i in plotting_groups])


    axes[i].scatter(tsne_result[:, 0], tsne_result[:, 1], s=30, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1')
    axes[i].set_ylabel('Component 2')


fig.suptitle("t-SNE",size = 30)

fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_tsne.png')


fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25,5))

for i in range(len(data_splits["subtraining"])):
    
    plotting_groups = np.zeros(len(input_data))
    plotting_groups[np.array(data_splits["testing"].index)] = 1
    plotting_groups[np.array(data_splits["validation"][i].index)] = 2
    plotting_groups = np.array([int(i) for i in plotting_groups])


    axes[i].scatter(umap_result[:, 0], umap_result[:, 1], s=30, c=colors[plotting_groups])
    axes[i].set_xlabel('Component 1')
    axes[i].set_ylabel('Component 2')

fig.suptitle("UMAP", size = 30)
fig.tight_layout(h_pad = 2)                                                          #tight layout so there is no overlay between plots
plt.savefig('train_val_umap.png')



medin_disorder = pd.read_excel("Medin_disorder_prediction.xlsx")

def add_secondary_xaxis(ax, positions, labels):
    secax = ax.secondary_xaxis('bottom')
    secax.set_xticks(positions)
    secax.set_xticklabels(labels, rotation='horizontal', fontsize=15)
    secax.spines['bottom'].set_position(('outward', 20))  # Adjust the position of the secondary x-axis
    secax.set_xlabel('Residue', size = 15)
    
fig, ax = plt.subplots(figsize=(25, 8))
ax.plot(medin_disorder["POS"], medin_disorder["IUPRED SCORE"],color="#5BC8AF", linewidth = 3, label = "Short disorder prediction")
ax.plot(medin_disorder["POS"], medin_disorder["ANCHOR SCORE"],color="#087F8C", linewidth = 3, label = "Disordered binding region")
ax.set_ylabel('Probability of Disorder')
ax.set_title('Medin Disorder Prediction')
ax.set_xticks(medin_disorder["POS"])
ax.set_xticklabels(range(1, len(medin_disorder)+1))
add_secondary_xaxis(ax, medin_disorder["POS"], medin_disorder["AMINO ACID"])
ax.legend()
plt.grid(False)

     


