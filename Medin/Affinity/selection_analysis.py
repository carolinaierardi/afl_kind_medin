import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['axes.labelsize'] = 15

palette = sns.color_palette("flare", 3)
shade1, shade2, shade3 = palette

df = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/final_Medin_selection.csv")

df.index = [f"M{i}" for i in range(1, 31)]

# Plotting
fig, ax = plt.subplots(figsize=(7, 12))
num_molecules = len(df)
bar_width = 0.2
r1 = np.arange(num_molecules)
r2 = [x + bar_width for x in r1]
r3 = [x + bar_width for x in r2]

# Make the plot
ax.barh(r1, df['old affinity'], color=shade3, height=bar_width, edgecolor='k', label='Affinity Model 1')
ax.barh(r2, df['new affinity'], color=shade2, height=bar_width, edgecolor='k', label='Affinity Model 2')
ax.barh(r3, df['QED'], color=shade1, height=bar_width, edgecolor='k', label='QED')

# Add labels to the y-axis
plt.yticks([r + bar_width for r in range(num_molecules)], df.index)

# Add xlabel, ylabel, title and legend
ax.set_xlabel('Metrics')
ax.set_ylabel('Molecules')
ax.set_title('Selected molecule metrics')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, frameon=False)

# Show the plot
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/selection_data.png")
