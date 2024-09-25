#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:44:29 2024

@author: carolinaierardi
"""

import os
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle



#%% Import data

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")

with open("preprocessed_kinetics_data.pkl", "rb") as f:
    input_data = pickle.load(f)
    
with open("data_splits.pkl", "rb") as f: 
    data_splits = pickle.load(f)
    

# filtered_df = input_data[(input_data['target_id'].str.len() == 6) & (input_data['domains'] == '(full)')]

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data_preprocess/Figures")

#%% Functions needed 



final_data = input_data[["Complex_ID", "target", "Protein_ID","processed_sequence", "n_residue", 
            "Ligand_ID","canonical_smiles", 
            "pKD", "pkon", "pkoff"]]


with open("final_data.pkl", "wb") as f:
     pickle.dump(final_data, f)

training_data = final_data.loc[data_splits["training"].index]
testing_data = final_data.loc[data_splits["testing"].index]

#CV data

subtraining_data = [final_data.loc[data_splits["subtraining"][i].index] for i in range(len(data_splits["subtraining"]))]
validation_data = [final_data.loc[data_splits["validation"][i].index] for i in range(len(data_splits["validation"]))]

processed_splits = {"training": training_data,"testing":testing_data,
               "subtraining":subtraining_data,"validation":validation_data}

hfont = {"fontname":"Arial"}
plt.rcParams['font.family'] = 'Arial'


#Make figures
plt.figure()
plt.hist(training_data["n_residue"], ec = "k", color = "#479ad1", bins = 30, label = "training")
plt.hist(testing_data["n_residue"], ec = "k", color = "#a6cde9", bins = 30, label = "testing")
plt.xlabel("Protein n residue", **hfont, size = 12)
plt.ylabel("Frequency", **hfont, size = 12)
plt.legend(loc='upper right') 
plt.title("Number of residues", **hfont, size = 15)

plt.savefig("protein_lengths.png")


os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")

with open("processed_splits.pkl", "wb") as f:
     pickle.dump(processed_splits, f)
     
with open("processed_splits.pkl", "rb") as f:
     processed_splits = pickle.load(f)
     



     







