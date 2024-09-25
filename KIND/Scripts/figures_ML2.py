#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:41:23 2024

@author: carolinaierardi
"""


import pandas as pd
import os
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools

def load_data(split):
    training_data = pd.read_pickle("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/kind/data/pickles/processed_splits.pkl")
    kd_data = pd.read_csv(f"/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/kind/data/CSV/Epoch108/kd_split{split}_epoch108.csv")
    kon_data = pd.read_csv(f"/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/kind/data/CSV/Epoch108/kon_split{split}_epoch108.csv")
    koff_data = pd.read_csv(f"/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/kind/data/CSV/Epoch108/koff_split{split}_epoch108.csv")
    assert len(kd_data) == len(kon_data) == len(koff_data), "Length of kd_data, kon_data, and koff_data must be the same."
    split5_data = pd.merge(pd.merge(kon_data, koff_data, on = "id"), kd_data, on = "id")
    split5_data = split5_data.drop_duplicates(subset=['id'], keep='first')
    return split5_data, training_data



def calculate_metrics(true_values, predicted_values):
    metrics = {}
    metrics['pearson_r'] = pearsonr(true_values, predicted_values)[0]
    metrics['spearman_rho'] = spearmanr(true_values, predicted_values)[0]
    metrics['r2'] = r2_score(true_values, predicted_values)
    metrics['mae'] = mean_absolute_error(true_values, predicted_values)
    metrics['mse'] = mean_squared_error(true_values, predicted_values)
    return metrics





output_data, training_data = load_data(2)
transformed_data = inverse_transform(output_data, training_data, 2)
transformed_data["Pred. Calc. Kd"] = transformed_data["Predicted Koff"] - transformed_data["Predicted Kon"]
transformed_data["True Calc. Kd"] = transformed_data["True Koff"] - transformed_data["True Kon"]

aff_data = transformed_data.iloc[:,-4:]

corr = np.corrcoef(aff_data.T)

df_corr = pd.DataFrame(corr, index=aff_data.columns, columns=aff_data.columns)
labels = ["Pred Kd", "True Kd", "Pred Calc Kd", "True Calc Kd"]

os.chdir("rds/project/rds-mphil-rkVIP0ZCj0k/carolina/kind/data/processing_scripts")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 20





pred_kon_KIND_recoverd = pred_kon_KIND*1.21-5.00
pred_koff_KIND_recoverd = pred_koff_KIND*0.88+1.64





os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/KIND_PDB/Epoch_14")
kaff_data = pd.read_csv("kaff_KIND.csv")
kd_data = pd.read_csv("aff_KIND.csv")
kon_data = pd.read_csv("kon_KIND.csv")
koff_data = pd.read_csv("koff_KIND.csv")
assert len(kd_data) == len(kon_data) == len(koff_data) == len(kaff_data), "Length of kd_data, kon_data, and koff_data must be the same."

aff_data = pd.DataFrame(data = {"True Affinity": kd_data["True Affinity KIND"],
                                "Predicted Affinity": kd_data["Predicted Affinity KIND"],
                                "Predicted Calc. Affinity": kaff_data["Predicted Kaff KIND"]})

columns = aff_data.columns
pairs = list(itertools.combinations(columns, 2))

# Set up the color palette
palette = sns.color_palette("crest", 2)
color = palette[1]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for (i, (col1, col2)) in enumerate(pairs):
    ax = axes[i]
    sns.scatterplot(x=aff_data[col1], y=aff_data[col2], color=color, ax=ax)
    ax.set_title(f'r = {aff_data[col1].corr(aff_data[col2]):.2f}')

plt.suptitle("KIND", size = 30)
plt.tight_layout()
plt.savefig("corr_scatter_KIND.png")


kaff_data = pd.read_csv("kaff_PDB.csv")
kd_data = pd.read_csv("aff_PDB.csv")
assert len(kd_data) == len(kon_data) == len(koff_data) == len(kaff_data), "Length of kd_data, kon_data, and koff_data must be the same."

aff_data = pd.DataFrame(data = {"True Affinity": kd_data["True Affinity"],
                                "Predicted Affinity": kd_data["Predicted Affinity"],
                                "Predicted Calc. Affinity": kaff_data["Predicted Kaff"]})


columns = aff_data.columns
pairs = list(itertools.combinations(columns, 2))

# Set up the color palette
palette = sns.color_palette("crest", 2)
color = palette[1]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for (i, (col1, col2)) in enumerate(pairs):
    ax = axes[i]
    sns.scatterplot(x=aff_data[col1], y=aff_data[col2], color=color, ax=ax)
    ax.set_title(f'r = {aff_data[col1].corr(aff_data[col2]):.2f}')

plt.suptitle("PDBbind", size = 30)
plt.tight_layout()
plt.savefig("corr_scatter_PDB.png")


os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/KIND_PDB")

colours = ["#FB6107", "#F3DE2C"]


def make_plot(files, metric): 
    
    k = ["pKon", "pKoff"]
    
    for ii, file in enumerate(files): 
        df = pd.read_csv(file)
        df = df.iloc[:,4]
        
        plt.figure(figsize=(8, 6))
        plt.plot(df.index, df, c = "#FB6107", linewidth = 2)
            
        plt.xlabel('Epochs')
        plt.title(metric + f" {k[ii]}")
        plt.ylabel(metric)
        plt.show()
        


files_r2 = ["r2_kon_KIND.csv", "r2_koff_KIND.csv"]
files_pearson = ["pearson_kon_KIND.csv", "pearson_koff_KIND.csv"]
files_mae = ["mae_kon_KIND.csv", "mae_koff_KIND.csv"]


make_plot(files_r2, r'$R^2$ score')
make_plot(files_pearson, "Pearson's r")
make_plot(files_mae,'MAE')

def make_plot(files, metric): 
        
    df = pd.read_csv(files[0])
    df_p = pd.read_csv(files[1])
    
    df = df.iloc[:,4]
    df_p = df_p.iloc[:,4]
    
    plt.figure(figsize=(8, 6))
    plt.plot(df.index, df, c = "#FB6107", linewidth = 2, label = "KIND")
    plt.plot(df.index, df_p, c = "#5C8001", linewidth = 2, label = "PDBbind")
        
    plt.xlabel('Epochs')
    plt.title(metric + " pKaff")
    plt.ylabel(metric)
    plt.legend()
    plt.show()


files_r2 = ["r2_aff_KIND.csv", "r2_PDB.csv"]
files_pearson = ["pearson_aff_KIND.csv", "pearson_PDB.csv"]
files_mae = ["mae_aff_KIND.csv", "mae_PDB.csv"]


make_plot(files_r2, r'$R^2$ score')
make_plot(files_pearson, "Pearson's r")
make_plot(files_mae,'MAE')


files_pearson = ["pearson_aff_KIND.csv", "pearson_kaff.csv"]
make_plot(files_pearson, "Pearson's r")

#get final metrics


os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/KIND_PDB/Epoch_14")
kaff_data = pd.read_csv("kaff_KIND.csv")
kd_data = pd.read_csv("aff_KIND.csv")
kon_data = pd.read_csv("kon_KIND.csv")
koff_data = pd.read_csv("koff_KIND.csv")


def calculate_metrics(predicted_values,true_values):
    metrics = {}
    metrics['pearson_r'] = pearsonr(true_values, predicted_values)[0]
    metrics['spearman_rho'] = spearmanr(true_values, predicted_values)[0]
    metrics['r2'] = r2_score(true_values, predicted_values)
    metrics['mae'] = mean_absolute_error(true_values, predicted_values)
    metrics['mse'] = mean_squared_error(true_values, predicted_values)
    return metrics

calculate_metrics(kaff_data.iloc[:,0], kaff_data.iloc[:,1])

kaff_data = pd.read_csv("kaff_PDB.csv")
kd_data = pd.read_csv("aff_PDB.csv")

calculate_metrics(kaff_data.iloc[:,0], kaff_data.iloc[:,1])





