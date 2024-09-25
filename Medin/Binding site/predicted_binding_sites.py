import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import glob

# Path to the main directory containing the subdirectories
main_directory = "/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/inf_topA/Medin_outputs/plots/"

sequence = "RLDKQGNFNAWVAGSYGNDQWLQVDLGSSKEVTGIITQGARNFGSVQFVA"

positions = list(range(len(sequence)))
labels = list(sequence)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 20

# Lists to store the data for plotting
x_data_list = []
y_data_column2_list = []
y_data_column3_list = []

# Iterate through each subdirectory and process the .csv file
for sub_dir in os.listdir(main_directory):
    sub_dir_path = os.path.join(main_directory, sub_dir)
    
    if os.path.isdir(sub_dir_path):
        # Assume there's only one .csv file in each subdirectory
        csv_file_path = glob.glob(os.path.join(sub_dir_path, "*.csv"))[0]
        
        # Read the first three columns of the .csv file
        df = pd.read_csv(csv_file_path, usecols=[0, 1, 2])
        
        # Store the data
        x_data_list.append(df.iloc[:, 0])
        y_data_column2_list.append(df.iloc[:, 1])
        y_data_column3_list.append(df.iloc[:, 2])

x_data_df = pd.concat(x_data_list, axis=1)
y_data_column2_df = pd.concat(y_data_column2_list, axis=1)
y_data_column2_df["residue"] = x_data_list[0] + 1
y_data_column3_df = pd.concat(y_data_column3_list, axis=1)
y_data_column3_df["residue"] = x_data_list[0] + 1

average_probabilities = y_data_column3_df[y_data_column3_df.columns[:-1]].mean(axis=1)
average_prob = pd.DataFrame(data = {"residue": x_data_list[0] + 1, "average prob": average_probabilities})

binding_threshold = 0.8
required_count = len(y_data_column3_df.columns[1:]) / 2  # Half of the molecules
filtered_df = y_data_column3_df[y_data_column3_df.iloc[:, 1:].apply(lambda row: (row > binding_threshold).sum(), axis=1) >= required_count]

highest_probs = average_prob[average_prob["average prob"] > 0.8]

os.chdir("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin")
average_prob.to_csv("probabiltiy_per_res.csv", index = False)
highest_probs.to_csv("highest_prob.csv", index = False)

def add_secondary_xaxis(ax, positions, labels):
    secax = ax.secondary_xaxis('bottom')
    secax.set_xticks(positions)
    secax.set_xticklabels(labels, rotation='horizontal', fontsize=15)
    secax.spines['bottom'].set_position(('outward', 20))  # Adjust the position of the secondary x-axis
    secax.set_xlabel('Residue', size = 15)

fig, ax = plt.subplots(figsize=(25, 8))
for i, y_data in enumerate(y_data_column2_list):
    ax.plot(x_data_list[i], y_data,color="#5BC8AF", alpha = 0.5)
ax.set_ylabel('Residue distance to ligand (A)')
ax.set_title('Residue distance to ligand')
ax.set_xticks(positions)
ax.set_xticklabels(range(1, len(positions)+1))
add_secondary_xaxis(ax, positions, labels)
plt.grid(False)
plt.savefig("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/min_res_dist2.png")

fig, ax = plt.subplots(figsize=(25, 8))
for i, y_data in enumerate(y_data_column3_list):
    ax.plot(x_data_list[i], y_data,color="#25CED1", alpha = 0.5)
ax.set_ylabel('Probability of binding to ligand')
ax.set_title('Probability of binding per residue')
ax.set_xticks(positions)
ax.set_xticklabels(range(1, len(positions)+1))
add_secondary_xaxis(ax, positions, labels)
plt.grid(False)
plt.savefig("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/prob_contact2.png")
