# # Load the CSV file
# input_csv_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/fda_proteome/data/human_proteome_sorted.csv'
# df = pd.read_csv(input_csv_path)
# df.rename(columns={"Unnamed: 0": "id", "id":"uniprot_id","len": "resi_num"}, inplace = True)

# # Directory containing the .pkl files
# pkl_dir = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/mol_feat_Medin'
# pkl_files = [f for f in os.listdir(pkl_dir) if f.endswith('.pickle')]

# # Ensure we have exactly 30 pkl files
# if len(pkl_files) != 30:
#     raise ValueError("The directory must contain exactly 30 .pickle files.")

# # Sort the files to ensure consistent order
# pkl_files.sort()

# mapping = []

# os.chdir("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/proteome_input")
# # Create new CSV files with the added column
# for idx, pkl_file in enumerate(pkl_files):
#     pkl_file_path = os.path.join(pkl_dir, pkl_file)
#     df['mol_feat_file'] = pkl_file_path
#     output_csv_path = f'input_{idx + 1}.csv'
#     df.to_csv(output_csv_path, index=False)
#     print(f"Created {output_csv_path}")

#     pkl_number = os.path.splitext(pkl_file)[0].strip('[]')
#     mapping.append([idx + 1, pkl_number])

# mapping_df = pd.DataFrame(mapping, columns=['Input_Number', 'Pickle_ID'])
# mapping_csv_path = 'mapping.csv'
# mapping_df.to_csv(mapping_csv_path, index=False)
# print(f"Created mapping file: {mapping_csv_path}")

# print("All files created successfully.")

# print("All files created successfully.")

from rdkit import Chem
import pandas as pd
import os

df = pd.read_csv("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/final_Medin_selection.csv")

def count_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    hydrogen_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    non_hydrogen_count = mol.GetNumAtoms()  # RDKit's GetNumHeavyAtoms() counts non-hydrogen atoms
    return hydrogen_count + non_hydrogen_count

df['n_atoms'] = df['SMILES'].apply(count_atoms)
df.sort_values(by = "n_atoms", ascending = False, inplace = True)

df.to_csv("/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/sorted_mol_select.csv")