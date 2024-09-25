import os 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import QED, Draw, AllChem, DataStructs
import numpy as np

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 40
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 20


def calculate_qed(smiles):
    mol = Chem.MolFromSmiles(smiles)
    qed = QED.qed(mol)
    return qed


def plot_molecule_grid(df, filepath):
    n_cols = 5  # Number of columns in the plot grid
    n_rows = len(df) // n_cols + (len(df) % n_cols > 0)  # Number of rows needed
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4 * n_rows))
    
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        img = Draw.MolToImage(mol, size=(200, 200))
        
        row_idx = idx // n_cols
        col_idx = idx % n_cols
        ax = axes[row_idx, col_idx] if n_rows > 1 else axes[col_idx]
        
        ax.imshow(img)
        ax.axis('off')
        title = f"ID: {row['id']}\nAffinity New: {row['new affinity']}\nAffinity Old: {row['old affinity']}\nQED: {round(row['QED'],3)}"
        ax.set_title(title, fontsize=20)
    
    # Remove empty subplots
    for ax in axes.flat:
        if not ax.has_data():
            fig.delaxes(ax)

    plt.tight_layout()
    plt.savefig(filepath)

clusters_all_nohype = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/v2/butina_meta.csv")

os.chdir("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/v2")

plt.figure(figsize=(20, 8))
sns.set(font_scale = 3)

# Create the scatter plot with adjustments
scatter_plot = sns.scatterplot(
    data=clusters_all_nohype,
    x='cluster',
    y='affinity',
    hue='affinity',  # Color by affinity
    size='QED',  # Size by scaled QED
    palette='flare',  # Color palette
    sizes=(20, 500),  # Size range for the points
    alpha=0.7,  # Transparency for overlapping points
    legend=False  # Do not show legend

)

sizes = [30, 60, 120, 150, 300]  # Adjust these sizes based on your QED_scaled values
legend_sizes = [size / 300 for size in sizes]  # Rescale for the legend
scatter_legend = scatter_plot.legend(
    handles=[
        plt.scatter([], [], s=sizes[size], 
                        label=f'{legend_sizes[size]}', 
                        color='black', alpha=0.5) for size in range(len(legend_sizes))
    ],
    title='QED',
    loc='upper left',  # Adjust location as needed
    bbox_to_anchor=(1.02, 1),  # Place legend outside the plot
    borderaxespad=0.
)


# Set plot labels and title
plt.xlabel('Clusters')
plt.ylabel('Affinity Score')
plt.title('Model 1')

# Adjust layout to make space for the colorbar
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("butina_50k.png")


#Molecule selection

#It was decided that choosing molecules based on their affinity alone is the best method

highest_affinity = clusters_all_nohype.loc[clusters_all_nohype.groupby('cluster')['affinity'].idxmax()]
highest_affinity.sort_values(by = "affinity", ascending = False, inplace = True)
selected_mol_old = highest_affinity.head(10)
selected_mol_old.to_csv("selected_mol_old_final.csv", index = False)

clusters_all_nohype['Selected'] = clusters_all_nohype['ID'].isin(selected_mol_old['ID'])


df_selected = clusters_all_nohype[clusters_all_nohype['Selected']]
df_non_selected = clusters_all_nohype[~clusters_all_nohype['Selected']]

#Method 1 - select the molecules with highest combination score

plt.figure(figsize=(15, 9))
sns.set(font_scale = 3)
# Plot all molecules in gray
scatter_plot = sns.scatterplot(
    data=df_non_selected,
    x='cluster',
    y='affinity',
    size='QED',  # Size by scaled QED
    color='gray',  # Color palette
    sizes=(20, 500),  # Size range for the points
    alpha=0.7,  # Transparency for overlapping points
    legend=False  # Do not show legend

)

scatter_plot = sns.scatterplot(
    data=df_selected,
    x='cluster',
    y='affinity',
    size='QED',  # Size by scaled QED
    color=(0.42355299, 0.16934709, 0.42581586),
    sizes=(20, 500),  # Size range for the points
    alpha=0.9,  # Transparency for overlapping points
    legend=False  # Do not show legend

)


plt.title('Molecule Selection')
plt.xlabel('Cluster')
plt.ylabel('Affinity Score')
plt.savefig("butina_selected_old_final.png")



#%% New model

clusters_all = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/hype/butina_hype_meta.csv")

os.chdir("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype")

plt.figure(figsize=(20, 8))
sns.set(font_scale = 3)

# Create the scatter plot with adjustments
scatter_plot = sns.scatterplot(
    data=clusters_all,
    x='cluster',
    y=' pred_affinity',
    hue=' pred_affinity',  # Color by affinity
    size='QED',  # Size by scaled QED
    palette='flare',  # Color palette
    sizes=(20, 500),  # Size range for the points
    alpha=0.7,  # Transparency for overlapping points
    legend=False  # Do not show legend

)

sizes = [30, 60, 120, 150, 300]  # Adjust these sizes based on your QED_scaled values
legend_sizes = [size / 300 for size in sizes]  # Rescale for the legend
scatter_legend = scatter_plot.legend(
    handles=[
        plt.scatter([], [], s=sizes[size], 
                        label=f'{legend_sizes[size]}', 
                        color='black', alpha=0.5) for size in range(len(legend_sizes))
    ],
    title='QED',
    loc='upper left',  # Adjust location as needed
    bbox_to_anchor=(1.02, 1),  # Place legend outside the plot
    borderaxespad=0.
)


# Set plot labels and title
plt.xlabel('Cluster')
plt.ylabel('Affinity Score')
plt.title('Model 2')

# Adjust layout to make space for the colorbar
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("butina_hype_50k.png")

#Model comparison

top1k_inter = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/hype/top_1k_intersection_ids.csv")
top10k_inter = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/hype/top_10k_intersection_ids.csv")

highest_affinity = clusters_all.loc[clusters_all.groupby('cluster')[' pred_affinity'].idxmax()]
highest_affinity.sort_values(by = " pred_affinity", ascending = False, inplace = True)


clusters_all['Combined_Score'] = clusters_all[' pred_affinity'] + clusters_all['QED']

highest_combined_score = clusters_all.loc[clusters_all.groupby('cluster')['Combined_Score'].idxmax()]
highest_combined_score.sort_values(by = "Combined_Score", ascending = False, inplace = True)


high_aff_10k = highest_affinity["id"].iloc[np.where(highest_affinity['id'].isin(top10k_inter["ID"]) == True)[0]]

palette = sns.color_palette("flare", 3)
shade1, shade2, shade3 = palette

def assign_color(molecule_id):
    if molecule_id in top1k_inter["ID"].values:
        return shade2
    elif molecule_id in high_aff_10k.values:
        return shade3
    elif molecule_id in highest_combined_score["id"].head(10).values:
        return (0.42355299, 0.16934709, 0.42581586)
    elif molecule_id in selected_mol_old["ID"]:
        print("Mol found")
        return "green"
    else:
        return 'gray'  # or any other default color

clusters_all['color'] = clusters_all['id'].apply(assign_color)

df_selected = clusters_all[clusters_all['color'] != 'gray']
df_non_selected = clusters_all[clusters_all['color'] == "gray"]

#Method 1 - select the molecules with highest combination score

plt.figure(figsize=(15, 9))
sns.set(font_scale = 3)

# Create the scatter plot with adjustments
scatter_plot = sns.scatterplot(
    data=df_non_selected,
    x='cluster',
    y=' pred_affinity',
    color = "gray",
    size='QED',  # Size by scaled QED
    sizes=(20, 500),  # Size range for the points
    alpha=0.7,  # Transparency for overlapping points
    legend=False  # Do not show legend

)

scatter_plot = sns.scatterplot(
    data=df_selected,
    x='cluster',
    y=' pred_affinity',
    hue=' pred_affinity',  # Color by affinity
    size='QED',  # Size by scaled QED
    palette='flare',  # Color palette
    sizes=(20, 500),  # Size range for the points
    alpha=0.9,  # Transparency for overlapping points
    legend=False  # Do not show legend

)

# Set plot labels and title
plt.xlabel('Cluster')
plt.ylabel('Affinity Score')
plt.title('Molecule Selection')

# Adjust layout to make space for the colorbar
plt.savefig("butina_hype_selected.png")

df_selected.to_csv("selected_mol_new.csv")

old_affinity = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/v2/finalv2_mapped_ids_affinities.csv")
new_affinity = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/data/processed/hype/final_hype_mapped_ids_affinities.csv")

all_selected_mol = pd.concat([selected_mol_old["ID"], df_selected["id"]], axis = 0, ignore_index = True).to_list()

old_affs = []
new_affs = []
smiles1 = []
smiles2 = []
qed = []

for i in all_selected_mol:
    smiles1 += [old_affinity["SMILES"].loc[old_affinity["ID"] == i].values[0]]
    smiles2 += [new_affinity[" smiles"].loc[new_affinity["id"] == i].values[0]]
    old_affs += [old_affinity["affinity"].loc[old_affinity["ID"] == i].values[0]]
    new_affs += [new_affinity[" pred_affinity"].loc[new_affinity["id"] == i].values[0]]
    qed += [calculate_qed(smiles2[-1])]



final_molecule_selection = pd.DataFrame(data = {"id": all_selected_mol, "smiles": smiles1, 
                                                "old affinity": old_affs,"new affinity": new_affs, 
                                                "QED": qed})


final_molecule_selection = pd.read_csv("/home/cm2231/rds/project/rds-fYhPa3It0hY/carolina/inferences/analysis/figures/figures/hype/final_Medin_selection.csv")
mols = final_molecule_selection['smiles'].apply(Chem.MolFromSmiles)

# Compute the fingerprints for each molecule
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in mols]

# Calculate the Tanimoto similarity matrix using BulkTanimotoSimilarity
num_molecules = len(fingerprints)
similarity_matrix = np.zeros((num_molecules, num_molecules))

for i in range(num_molecules):
    similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[i], fingerprints)
    similarity_matrix[i, :] = similarities

# Convert similarity to distance
similarity_df = pd.DataFrame(similarity_matrix, columns=final_molecule_selection.index, index=final_molecule_selection.index)



plt.figure(figsize=(4, 4))
plt.imshow(similarity_matrix, cmap = "flare")
plt.colorbar(shrink = 0.7)
plt.grid(False)
plt.title('Tanimoto Similarity Matrix')
plt.xlabel('Molecule Index')
plt.ylabel('Molecule Index')
plt.tight_layout()
plt.savefig("tanimoto_selection.png")


final_molecule_selection.to_csv("final_Medin_selection.csv")

plot_molecule_grid(final_molecule_selection, "final_selection.png")






