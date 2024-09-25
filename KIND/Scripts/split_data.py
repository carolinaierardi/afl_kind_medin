#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:01:03 2024

@author: carolinaierardi
"""


import os
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
    


#%% Import data
os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")


with open("preprocessed_kinetics_data.pkl", "rb") as f:
    input_data = pickle.load(f)
    


#%% Split into train and test

stratify_column = input_data["target_class"] + "_" + input_data["assay_method"]



test_set_size = int(np.round(np.floor(len(stratify_column)*0.1)))

training, testing = train_test_split(input_data, test_size=test_set_size-1, 
                                         shuffle = True,
                                     random_state = 42, 
                                     stratify=stratify_column) #take 20% of the data as testing set 


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


#%% Save to directory 

os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")

with open("data_splits.pkl", "wb") as f:
     pickle.dump(data_splits, f)
     
     
     
