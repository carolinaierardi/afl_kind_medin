#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 20:10:28 2024

@author: carolinaierardi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


#import data
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 20


os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kind1_CV")

r2_koff_stats = pd.read_csv("r2_kind_cv.csv")
r2_koff_stats_ = r2_koff_stats[["split_3 - validation/r2_score_koff", "split_5 - validation/r2_score_koff",
                      "split_4 - validation/r2_score_koff", "split_2 - validation/r2_score_koff",
                      "split_1 - validation/r2_score_koff"]]

r2_kon = pd.read_csv("r2_kon.csv")
r2_kon_stats_ = r2_kon[["split_3 - validation/r2_score_kon", "split_5 - validation/r2_score_kon",
                      "split_4 - validation/r2_score_kon", "split_2 - validation/r2_score_kon",
                      "split_1 - validation/r2_score_kon"]]


r2_aff = pd.read_csv("r2_aff.csv")
r2_aff_stats_ = r2_aff[["split_3 - validation/r2_score", "split_5 - validation/r2_score",
                      "split_4 - validation/r2_score", "split_2 - validation/r2_score",
                      "split_1 - validation/r2_score"]]


colours = ["#FB6107", "#F3DE2C", "#7CB518", "#5C8001", "#FBB02D"]

def make_plot(files, cols, metric): 
    
    k = ["pKD", "pKon", "pKoff"]
    
    for ii, file in enumerate(files): 
        df = pd.read_csv(file)
        df = df.iloc[:,cols]
        
        plt.figure(figsize=(8, 6))
        for i, column in enumerate(df.columns):
            plt.plot(df.index, df[column], label=f"Split {i + 1}", c = colours[i], linewidth = 2)
            
        plt.xlabel('Epochs')
        plt.title(metric + f" {k[ii]}")
        plt.ylabel(metric)
        plt.legend()
        plt.show()
        

    
files_r2 = ["r2_aff.csv", "r2_kon.csv", "r2_kind_cv.csv"]
files_pearson = ["pearson_aff.csv", "pearson_kon.csv", "pearson_koff.csv"]
files_mae = ["mae_aff.csv", "mae_kon.csv", "mae_koff.csv"]
cols = [4, 10, 16, 22, 28]


    
make_plot(files_r2, cols, r'$R^2$ score')
make_plot(files_pearson, cols, "Pearson's r")
make_plot(files_mae, cols,'MAE')






    




