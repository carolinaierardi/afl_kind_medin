import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 30
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 20


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


pdb_file_path = "/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/scripts/Medin_MD_Gillian.pdb"

# Create a PDB parser
parser = Bio.PDB.PDBParser(QUIET=True)

# Parse the PDB file
structure = parser.get_structure("protein", pdb_file_path)

# Assuming you want to calculate the distance between all pairs of residues in the first model
model = structure[0]  # Get the first model
residues = [residue for residue in model.get_residues() if Bio.PDB.is_aa(residue)]  # Get all amino acid residues

distance_matrix = np.zeros((len(residues), len(residues)))

# Calculate distances
for i, residue_one in enumerate(residues):
    for j, residue_two in enumerate(residues):
        if i < j:  # Only calculate upper triangle of matrix to avoid redundant calculations
            distance = calc_residue_dist(residue_one, residue_two)
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # Mirror the distance


# Plot the distance matrix
fig, ax = plt.subplots(figsize=(8, 8))
cax = ax.matshow(distance_matrix, cmap='winter', alpha = 0.8)
cbar = fig.colorbar(cax, shrink = 0.7)
cbar.ax.tick_params(labelsize=15)  # Increase label size for the color bar
ax.grid(False)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.xlabel('Residue Position', fontsize=20)
plt.ylabel('Residue Position', fontsize=20)
plt.title('Inter-residue Distance Matrix', fontsize=25)

plt.tight_layout()  # Adjust subplots to fit into figure area.
plt.savefig("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/MD_interres_dist.png")
