import os
import pandas as pd
import re


# Define the input and output directories
input_directory = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/proteome_inf'
output_directory = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/selected_mol_proteome'

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Dictionary to hold data for the new CSV files
new_data = {}

def extract_number(value):
    match = re.search(r'\[(\d+(\.\d+)?)\]', value)
    if match:
        return float(match.group(1))
    return None


# Read each CSV file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.csv'):
        filepath = os.path.join(input_directory, filename)
        uniprot_id = os.path.splitext(filename)[0]
        # Read the CSV file without headers
        df = pd.read_csv(filepath, header=None)
        # Iterate through the rows of the dataframe
        for index, row in df.iterrows():
            original_filename_value = extract_number(row[0])
            affinity_value = extract_number(row[1])
            # Initialize the list for this new filename if it doesn't exist
            if original_filename_value is not None and affinity_value is not None:
                new_filename = f"{original_filename_value}.csv"
                if new_filename not in new_data:
                    new_data[new_filename] = []
            # Append the new row to the corresponding list
            new_data[new_filename].append({'Uniprot': uniprot_id, 'affinity': affinity_value})

# Write the new CSV files
for new_filename, rows in new_data.items():
    output_filepath = os.path.join(output_directory, new_filename)
    new_df = pd.DataFrame(rows)
    new_df.to_csv(output_filepath, index=False, columns=['Uniprot', 'affinity'])

#make plot of affintiies

import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 40
plt.rcParams['axes.labelsize'] = 30
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 30
plt.rcParams['legend.fontsize'] = 20

# Define the input and output directories
output_directory = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/selected_mol_proteome'
final_selection_file = '/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/final_Medin_selection.csv'

# Load the final selection affinities
final_selection_df = pd.read_csv(final_selection_file)

fig, axes = plt.subplots(6, 5, figsize=(30, 30), sharex = True)
axes = axes.flatten()

# Counter for the subplot index
plot_idx = 0

for _, row in final_selection_df.iterrows():
    base_filename = float(row['id'])
    final_affinity = row['new affinity']
    filepath = os.path.join(output_directory, f"{base_filename}.csv")
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
        final_row = pd.DataFrame({'Uniprot': ['final_selection'], 'affinity': [final_affinity]})
        df = pd.concat([df, final_row], ignore_index=True)
        df_sorted = df.sort_values(by='affinity').reset_index(drop=True)
        p25 = df_sorted['affinity'].quantile(0.25)
        p50 = df_sorted['affinity'].quantile(0.50)
        p75 = df_sorted['affinity'].quantile(0.75)
        ax = axes[plot_idx]
        ax.scatter(df_sorted.index, df_sorted['affinity'], label='Proteome', color='#70A288')
        if 'Q08431' in df_sorted['Uniprot'].values:
            idx = df_sorted[df_sorted['Uniprot'] == 'Q08431'].index[0]
            ax.scatter(idx, df_sorted.loc[idx, 'affinity'], color='#C95D63', label='Lactadherin',s = 500)
        idx = df_sorted[df_sorted['Uniprot'] == 'final_selection'].index[0]
        ax.scatter(idx, df_sorted.loc[idx, 'affinity'], color='#EE8434', label='Medin', s = 500)
        ax.axhline(y=p25, color='k', linestyle='--', linewidth=1, label='25th Percentile' if plot_idx == 0 else "")
        ax.axhline(y=p50, color='k', linestyle='--', linewidth=1, label='50th Percentile' if plot_idx == 0 else "")
        ax.axhline(y=p75, color='k', linestyle='--', linewidth=1, label='75th Percentile' if plot_idx == 0 else "")
        ax.set_title(f'M{plot_idx + 1}')
        ax.set_xlabel('Proteins')
        ax.set_ylabel('Affinity')
        plot_idx += 1

# Adjust layout and show plot
plt.legend()
plt.tight_layout()
plt.savefig("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/figures/selectivity.png")

from scipy import stats
import numpy as np

def calculate_p_value(affinity_value, distribution):
    # ecdf = stats.cumfreq(distribution, numbins=1000)
    # cdf = ecdf.cumcount / len(distribution)
    # bin_edges = ecdf.lowerlimit + np.linspace(0, ecdf.binsize * ecdf.cumcount.size, ecdf.cumcount.size)
    # idx = np.searchsorted(bin_edges, affinity_value, side='right')
    # p_value = 1 - cdf[idx-1] if idx > 0 else 1.0
    p_value = (1 + len(np.where(affinity_value < distribution)[0])) / (len(distribution) + 1)
    return p_value


fig, axes = plt.subplots(6, 5, figsize=(30, 30), sharey=True)
axes = axes.flatten()

# Counter for the subplot index
plot_idx = 0

# Process each CSV file in the order specified by final_selection_Medin.csv
for _, row in final_selection_df.iterrows():
    base_filename = float(row['id'])
    final_affinity = row['new affinity']
    filepath = os.path.join(output_directory, f"{base_filename}.csv")
    if os.path.exists(filepath):
        # Read the CSV file
        df = pd.read_csv(filepath)
        # Extract affinity values
        affinities = df['affinity'].astype(float)
        # Plot the histogram
        ax = axes[plot_idx]
        ax.hist(affinities, bins=30, alpha=0.7, label='Proteome', color='#70A288')
        # Plot the final selection affinity as a dashed line
        ax.axvline(x=final_affinity, color='#EE8434', linestyle='--', linewidth=4, label='Medin')
        # Highlight the affinity value for Uniprot Q08431 if it exists
        if 'Q08431' in df['Uniprot'].values:
            q08431_affinity = df.loc[df['Uniprot'] == 'Q08431', 'affinity'].values[0]
            ax.axvline(x=q08431_affinity, color='#C95D63', linestyle='--', linewidth=4, label='Lactadherin')
        # Calculate the p-value for the final selection affinity
        p_value = calculate_p_value(final_affinity, affinities)
        #ax.text(0.95, 0.95, f'p-value: {p_value:.2e}', transform=ax.transAxes, ha='right', va='top', fontsize=8)
        # Set plot title and labels
        ax.set_title(f'M{plot_idx + 1}, p = {round(p_value,2)}')
        ax.set_xlabel('Affinity')
        ax.set_ylabel('Frequency')
        # Increment the plot index
        plot_idx += 1
        # Add legend

plt.legend()
# Adjust layout and show plot
plt.tight_layout()
plt.savefig("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/figures/selectivity_hist.png")
