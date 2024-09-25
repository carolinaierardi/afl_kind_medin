""" import os
import pandas as pd

# Define the input directory containing the original .csv files
input_directory = '/home/cm2231/rds/project/rds-a1NGKrlJtrw/zinc/split'
# Define the output directory where the modified .csv files will be saved
output_directory = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/split_zinc'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Loop through each file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.csv') and filename.startswith('all_clean_'):
        # Extract the index number from the filename
        index = filename.split('_')[-1].split('.')[0]

        # Construct the full file path
        file_path = os.path.join(input_directory, filename)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Keep only the first three columns
        df = df.iloc[:, [1, 2, 0]]
        
        # Rename the columns
        df.columns = ['id', 'Name', 'SMILES']
        
        # Construct the new filename with just the index
        new_filename = f'{index}.csv'
        new_file_path = os.path.join(output_directory, new_filename)
        
        # Save the modified DataFrame to the new file
        df.to_csv(new_file_path, index=False)

print("Processing complete. Check the output directory for the modified files.")

 """
import os
import pandas as pd
import numpy as np

def split_csv(file_path, output_dir, start_index):
    # Read the CSV file and keep only the first three columns
    df = pd.read_csv(file_path, usecols=[0, 1, 2])

    # Reorder the columns
    df = df[[df.columns[1], df.columns[2], df.columns[0]]]
    
    # Rename the columns
    df.columns = ["id", "Name", "SMILES"]
    
    # Determine the size of each split
    num_splits = 10
    rows_per_split = int(np.ceil(len(df) / num_splits))
    
    # Split the dataframe and save each part
    for i in range(num_splits):
        start_row = i * rows_per_split
        end_row = start_row + rows_per_split
        split_df = df.iloc[start_row:end_row]
        
        # Construct the output file name
        output_file_name = f"{start_index + i}.csv"
        output_file_path = os.path.join(output_dir, output_file_name)
        
        # Save the split dataframe to a CSV file
        split_df.to_csv(output_file_path, index=False)

def process_directory(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize the index for the output files
    current_index = 0
    
    # Loop through the files in the input directory
    for file_name in sorted(os.listdir(input_dir)):
        if file_name.startswith("all_clean_") and file_name.endswith(".csv"):
            file_path = os.path.join(input_dir, file_name)
            split_csv(file_path, output_dir, current_index)
            current_index += 10

if __name__ == "__main__":
    # Specify the input and output directories
    input_directory = "/home/cm2231/rds/project/rds-a1NGKrlJtrw/zinc/split"
    output_directory = "/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/split_zinc"
    
    # Process the directory
    process_directory(input_directory, output_directory)

print("Processing complete. Check the output directory for the modified files.")
