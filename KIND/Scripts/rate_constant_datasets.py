#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 10:28:21 2024

@author: carolinaierardi
"""

import os 
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import subprocess 
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from collections import Counter
import requests
import pypdb
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


os.chdir("/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/Kinetics_rate_constants")


fedorov = pd.read_csv("fedorov_et_al.csv")
kind = pd.read_excel("KIND.xlsx")


#%%% Functions needed

def reduce_dataset(original_data, binding_data):
    
    """

    Parameters
    ----------
    original_data : DataFrame
        Dataframe to be filtered.
    binding_data : Set of columns in the dataframe
        The columns in the dataframe you wish all rows to contain values in

    Returns
    -------
    DataFrame
        New filtered dataframe

    """
    
    df_copy = original_data.copy()
    incomplete = np.where(binding_data.isna().any(axis=1) == True)[0]
    
    if len(incomplete) == 0:
        return df_copy
    else: 
        df_copy_reduced = df_copy.drop(incomplete)
        return df_copy_reduced
    
def canonical_smiles(reduced_data, column_name):
    
    """
    Parameters
    ----------
    reduced_data : DataFrame
        Dataframe with some set of molecules.
    column_name : Column in dataframe
       Column with the description of those molecules.

    Returns
    -------
    reduced_data : TYPE
        Dataframe with added columns of canonical smiles and ID for smiles.

    """
    
    canonsmiles = [Chem.CanonSmiles(smiles) for smiles in reduced_data[column_name]]     
    reduced_data["canonical_smiles"] = canonsmiles
    reduced_data = reduced_data.assign(id=(reduced_data['canonical_smiles']).astype('category').cat.codes)
    
    
    return reduced_data

def get_protein_sequence(uniprot_id):
    
    """
    Parameters
    ----------
    uniprot_id : UNIPROT accession (str)
        Usually a 6 digit code corresponding to an entry on UNIPROT.

    Returns
    -------
    sequence : Entry sequence (str)
        Corresponding sequence to entry.

    """
    
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        sequence = ''.join(fasta_data.split('\n')[1:])
        return sequence
    else:
        return None
    
def count_hydrogens(mol):
    return sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

def add_natoms(df):

    # Convert SMILES to RDKit Mol objects
    df['molecules'] = df['canonical_smiles'].apply(Chem.MolFromSmiles)

    # Extract the number of atoms in each molecule
    df['num_atoms'] = df['molecules'].apply(lambda mol: mol.GetNumAtoms() if mol else None)
    df['num_hydro'] = df['molecules'].apply(lambda mol: count_hydrogens(mol) if mol else None)

    # Drop the temporary 'molecules' column
    df = df.drop(columns=['molecules'])

    # Reorder columns to place 'num_atoms' right after 'canonical_smiles'
    columns = list(df.columns)
    canonical_smiles_index = columns.index('canonical_smiles')
    columns.insert(canonical_smiles_index + 1, columns.pop(columns.index('num_atoms')))
    columns.insert(canonical_smiles_index + 2, columns.pop(columns.index('num_hydro')))
    df = df[columns]

    return df
    
def extract_sequence_from_pdb(file_path):
    
    """
    Parameters
    ----------
    file_path : path to .pdb file (str)
        where .pdb file is stored.

    Returns
    -------
    sequence : Corresponding sequence (str)
        Sequence for corresponding pdb file.

    """
    
    parser = PDBParser()
    structure = parser.get_structure('', file_path)
    sequence = []
    for model in structure:
        for chain in model:
            chain_seq = []
            for residue in chain:
                if 'CA' in residue:
                    chain_seq.append(residue.resname)
            # Convert three-letter amino acid codes to one-letter codes
            seq_str = ''.join([seq1(res) for res in chain_seq])
            sequence.append(SeqRecord(Seq(seq_str), id=chain.id))
    return sequence

# Function to convert three-letter amino acid codes to one-letter codes
def seq1(residue):
    
    """
    Parameters
    ----------
    residue : str
        residue found in file.

    Returns
    -------
    str
        corresponding letter for residue.

    """
    # Mapping of three-letter codes to one-letter codes
    d = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
         'GLN': 'Q', 'GLU': 'E', 'GLY': 'G','HIS': 'H', 'ILE': 'I', 
         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 
         'SER': 'S','THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    
    return d.get(residue, 'X')  # Return 'X' if residue not found in the dictionary


def remove_duplicates(entries):
    """
    

    Parameters
    ----------
    entries : list
        list of indices in existing dataframe where entries are repeated.

    Returns
    -------
    my_entry : DataFrame
        Single entry of DataFrame with the chosen segment.

    """
    
    #make smaller dataframe with the repeated entries
    dup_dat = [pd.DataFrame(new_dataframe.iloc[i] for i in entries)][0]
    
    #If there are no proteins with the PDB structures, I will chose the median values
    if not any(dup_dat["target_id"].str.len() == 4):
        
        #If there are more than two repeats
        if len(entries) > 2: 
            
            #get the median entry and use that 
            my_pKd = dup_dat["pKD"].median()
            my_entry = dup_dat[dup_dat["pKD"] == my_pKd]
        
        else:
            #use the one with the smallest pKD value
            my_pKd = dup_dat["pKD"].min()
            my_entry = dup_dat[dup_dat["pKD"] == my_pKd]
            
    else: #If there are proteins with PDB codes
        
    #Use that entry as it has more accurate 
        my_entry = dup_dat[dup_dat["target_id"].str.len() == 4]
        
    #But take some information about protein
        my_info = dup_dat[dup_dat["target_id"].str.len() != 4]
        
        my_entry.loc[my_entry.index, "target_class"] = my_info.loc[my_info.index[-1], "target_class"]
        my_entry.loc[my_entry.index, "assay_method"] = my_info.loc[my_info.index[-1], "assay_method"]

             
    return my_entry  

def update_unknowns(row):
    
    if row['assay_method'] == 'Unknown' or row['target_class'] == 'Unknown':
        target = row['target']
        if target in known_values_dict:
            if row['assay_method'] == 'Unknown':
                row['assay_method'] = known_values_dict[target]['assay_method']
            if row['target_class'] == 'Unknown':
                row['target_class'] = known_values_dict[target]['target_class']
    return row


def process_row(row):
    
    sequence = row['sequence']
    domain = row['domains']
    
    if domain == '(full)':
        seq_length = len(sequence)
        return pd.Series([sequence, seq_length])
    else:
        # Parse the range from the domain string
        range_str = domain.split(' ')[0]
        start, end = map(int, range_str.split('-'))
        
        new_sequence = sequence[start-1:end]
        new_seq_length = len(new_sequence)
        return pd.Series([new_sequence, new_seq_length])
    



#%% Preprocess data

#Transform values into log for fedorov
fedorov["pKD"] = -np.log10(fedorov["kd_nM"]/(10**9))
fedorov["pkon"] = -np.log10(fedorov["kon_s"])
fedorov["pkoff"] = -np.log10(fedorov["koff_s"])

#First we will obtain only the columns necessary in the data


kind_ = kind[["target", "smiles_washed", 
              "pKD","pkon","pkoff","assay_method", "target_class",]].copy()

fedorov_ = fedorov[["target", "smiles", "pKD", "pkon", "pkoff", "docked_pdb"]].copy()



#Now that we've selected only the columns we need, we will remove rows with incomplete data

fedorov_binding = fedorov[["pKD", "pkon", "pkoff"]]  #name data rows
kind_binding = kind[["pKD", 'pkon', "pkoff"]]            #name data rows

#remove rows with incomplete data
fedorov_filtered = reduce_dataset(fedorov_, fedorov_binding)
kind_filtered = reduce_dataset(kind_, kind_binding)


#%%This section is to find all sequences for KIND 

kind_info = pd.read_excel("KIND_sequence_info.xlsx")
kind_proteins = kind_info["Protein Name"].to_list()
kind_uniprotid = kind_info["UniProt ID"].to_list()
kind_domains = kind_info["Sequence"].to_list()

kind_mapping = dict(zip(kind_proteins, kind_uniprotid))
kind_used_seq = dict(zip(kind_proteins, kind_domains))

kind_filtered["target_id"]  = kind_filtered["target"].map(kind_mapping)
kind_filtered["domains"] = kind_filtered["target"].map(kind_used_seq)


#%% Sequences for targets

#Next, we will obtain sequences for the targets in the data

#For KIND, we use the manual curation information to get sequences from each accession code

# Fetch sequences for each ID and store in a new column
kind_filtered['sequence'] = kind_filtered['target_id'].apply(get_protein_sequence)

# For Fedorov et al., we use the PDB files available 
# Dictionary to store PDB IDs and their corresponding sequences
sequences = {}

# Directory containing subdirectories with PDB files
root_dir = '/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data/5206636/ci0c00450_si_001/structures'

# Traverse directories and find PDB files
for subdir, _, files in os.walk(root_dir):
    for file in files:
        if file.endswith('.pdb'):
            pdb_id = os.path.splitext(file)[0]
            file_path = os.path.join(subdir, file)
            sequences[pdb_id] = extract_sequence_from_pdb(file_path)


# If you need to flatten the sequences into a single sequence per PDB file, you can concatenate them
flattened_sequences = {pdb_id: ''.join([str(seq_record.seq) for seq_record in seq_list])
                       for pdb_id, seq_list in sequences.items()}

shortened_id_to_sequence = {pdb_id[:4]: sequence for pdb_id, sequence in flattened_sequences.items()}

fedorov_filtered['sequence'] = fedorov_filtered['docked_pdb'].map(shortened_id_to_sequence)


fedorov_filtered["domains"] = "(full)"
fedorov_filtered["target_class"] = "Unknown"
fedorov_filtered["assay_method"] = "Unknown"
fedorov_filtered = fedorov_filtered.rename(columns = {"docked_pdb": "target_id", "smiles":"smiles_washed"})


#Concatenate two dataframes
new_dataframe = pd.concat([kind_filtered, fedorov_filtered], ignore_index=True)
#There are some entries that are different but correspond to the same proteins
#We will change those to detect where the changes occur 
new_dataframe['target'] = new_dataframe['target'].replace('A2A', 'A2a')
new_dataframe['target'] = new_dataframe['target'].replace('p38 MAP kinase', 'MAP38')
new_dataframe['target'] = new_dataframe['target'].replace('hCAII', 'hCA II')
new_dataframe['target'] = new_dataframe['target'].replace('Histamine H1 receptor', 'H1')
new_dataframe['target'] = new_dataframe['target'].replace('CDK2 Kinase', 'CDK2')
new_dataframe['target'] = new_dataframe['target'].replace('Muscarinic M3', 'M3')


#%% Canonical smiles

new_dataframe = canonical_smiles(new_dataframe, "smiles_washed")


#%% Remove duplicates

#create unique complex ID

entry_id = np.array(new_dataframe["target"] + "_" + new_dataframe["id"].astype(str))

entry_id_l = [i.lower() for i in entry_id]

unique_ids = set()
duplicate_ids = [x for x in entry_id_l if x in unique_ids or unique_ids.add(x)]  
duplicate_ids = np.unique(duplicate_ids)

duplicate_locs = [list(np.where(np.array(entry_id_l) == j)[0]) for j in duplicate_ids]
                          
decided_dat = [remove_duplicates(i) for i in duplicate_locs]
decided_dat = pd.concat(decided_dat, ignore_index=True)

unique_new = new_dataframe.drop(sum(duplicate_locs, []))
unique_new = pd.concat([unique_new, decided_dat], ignore_index=True)


known_values = unique_new[unique_new['assay_method'] != 'Unknown'].groupby('target').agg({
    'assay_method': 'first',
    'target_class': 'first'
}).reset_index()

known_values_dict = known_values.set_index('target').to_dict(orient='index')

# Apply the function to the DataFrame
unique_new = unique_new.apply(update_unknowns, axis=1)

unique_new[['processed_sequence', 'n_residue']] = unique_new.apply(process_row, axis=1)


#%% Add Ids to rows 

#create ids for the proteins and ligands 
unique_protein_names = unique_new['processed_sequence'].unique()
protein_id_dict = {protein: f'P{idx+1}' for idx, protein in enumerate(unique_protein_names)}

# Map the protein names in the DataFrame to their unique IDs
unique_new['Protein_ID'] = unique_new['processed_sequence'].map(protein_id_dict)
unique_new["Ligand_ID"] = "L" + unique_new['id'].astype(str)


# Create a new column 'Complex' by combining 'Protein' and 'Ligand'
unique_new['Complex'] = unique_new['Protein_ID'] + '_' + unique_new["Ligand_ID"]
unique_complexes = unique_new['Complex'].unique()
complex_id_dict = {complex_: f'C{idx+1}' for idx, complex_ in enumerate(unique_complexes)}

unique_new['Complex_ID'] = unique_new['Complex'].map(complex_id_dict)
unique_new.drop(columns=['Complex'], inplace=True)

unique_new.drop_duplicates(subset=['Complex_ID'], inplace=True)

input_data = unique_new[["Complex_ID", "target", "Protein_ID", "processed_sequence","n_residue",
                         "canonical_smiles", "Ligand_ID",
                         "pKD","pkon",'pkoff',"assay_method", "target_class"]]

input_data = add_natoms(input_data)

input_data.reset_index(drop = True, inplace = True)

with open("preprocessed_kinetics_data.pkl", "wb") as f:
     pickle.dump(input_data, f)
     






