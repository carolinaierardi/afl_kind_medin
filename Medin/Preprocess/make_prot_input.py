import os
import pandas as pd
from Bio import SeqIO

def process_fasta(fasta_file_path, output_csv_path, protein_feat_prefix):
    # Initialize lists to store the sequences and their lengths
    sequences = []
    lengths = []
    prefixes = []

    # Parse the FASTA file and extract sequences and their lengths
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        sequences.append(str(record.seq))
        lengths.append(len(record.seq))
        prefixes.append(protein_feat_prefix)

    # Create a DataFrame from the extracted data
    df = pd.DataFrame({
        'seq': sequences,
        'resi_num': lengths,
        'protein_feat_prefix': prefixes
    })

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv_path, index=False)

    print(f"Processing complete. Output saved to {output_csv_path}")

# Define the input FASTA file path, output CSV file path, and protein feature prefix
fasta_file_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/medin.fasta'
output_csv_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/Medin.csv'
protein_feat_prefix = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/output_Medin_prot/Medin'

# Process the FASTA file and generate the CSV
process_fasta(fasta_file_path, output_csv_path, protein_feat_prefix)
