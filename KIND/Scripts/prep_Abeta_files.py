#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:15:09 2024

@author: carolinaierardi
"""

import pandas as pd
import numpy as np 
import os
from Bio import SeqIO



#variables to change
current_wd = "/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data"

#for molecular feature extraction
screening_lib = "zinc-cayman_lib.xlsx"
output_csv = "input_Ab.csv"

#for protein feature extraction
input_fasta = "input.fasta"
output_fasta = "input_Ab.fasta"

#for inference
directory_prot_feat = "test_Ab/output_Ab/"
directory_mol_feat = "/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data"



#make molecular and protein features files
os.chdir(current_wd)                              #change wd
my_lib = pd.read_excel(screening_lib)             #import screening library
my_lib_red = my_lib.drop(["QED", 
                          "cns_guacamol","complexity_rdkit"], 
                         axis = 1)                #remove unneccessary columns
my_lib_red.insert(0, "id",my_lib_red.index, True) #add index as in the example code
my_lib_red.to_csv(output_csv)                     #save as a csv file


fasta_sequences = list(SeqIO.parse(input_fasta,'fasta'))

#CHANGE THIS DEPENDING ON WHICH SEQUENCE IS NEEDED
target_seq = fasta_sequences[-1]
SeqIO.write(target_seq, output_fasta, "fasta")


#write input .csv file for inference

#inference_df = pd.DataFrame(columns = ["id", "seq","resi_num",
 #                       "protein_feat_prefix", "mol_feat_file"])


#https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
mol_feat = []
for file in os.listdir(directory_mol_feat):
    filename = os.fsdecode(file)
    if filename.endswith(".pickle"): 
        mol_feat += [os.path.join(directory_mol_feat, filename)]
        continue
    else:
        continue

prot_seq = str(target_seq.seq)
num_resi = len(prot_seq)
prot_feat_prefix = str(directory_prot_feat + target_seq.id)

