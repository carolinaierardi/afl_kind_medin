import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20

# Load the CSV files
average_distances = pd.read_csv('/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/AF_medin/Medin_mean_0.csv', header=None, usecols = range(1,51), skiprows = 1)
stddev_distances = pd.read_csv('/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/AF_medin/Medin_std_0.csv', header=None, usecols = range(1,51), skiprows = 1)

# Normalize the standard deviations to use as alpha values
# Ensuring values are between 0 and 1 for transparency effect
stddev_norm = (stddev_distances - stddev_distances.min()) / (stddev_distances.max() - stddev_distances.min())

stddev_norm = 1 - stddev_norm
# # Create the heatmap for average distances
# plt.figure(figsize=(12, 10))
# sns.set(font_scale = 2)
# heatmap = sns.heatmap(average_distances, cmap='mako', annot=False)

# # Overlay the standard deviation as transparency
# # Iterate over the data to set alpha values
# for i in range(average_distances.shape[0]):
#     for j in range(average_distances.shape[1]):
#         heatmap.add_patch(
#             plt.Rectangle(
#                 (j, i), 1, 1, 
#                 fill=True, 
#                 color='white', 
#                 alpha=stddev_norm.iat[i, j]
#             )
#         )

# plt.title('Inter-Residue Distances (Å)')
# plt.xlabel('Residue Index')
# plt.ylabel('Residue Index')
# plt.savefig("Medin_distmat.png")


plt.figure(figsize=(8, 8))
dist = plt.imshow(average_distances, cmap='magma')
plt.colorbar(dist, shrink = 0.7)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title("Mean inter-residue distance")
plt.tight_layout()
plt.savefig("Medin_distmat.png")

plt.figure(figsize=(8, 8))
dist = plt.imshow(stddev_distances, cmap='magma')
plt.colorbar(dist, shrink = 0.7)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title("SD inter-residue distance")
plt.tight_layout()
plt.savefig("Medin_distmat_sd.png")