import os
import csv
import pandas as pd
import re


# Define the paths
existing_csv_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/Medin.csv'
directory_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/mol_feat_Medin'
output_csv_path = '/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/input_inf_top.csv'

with open(existing_csv_path, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)
    single_row = next(reader)

# Step 2: Loop through the directory to find all .pickle files and extract IDs
pickle_files = []
ids = []

# Regex pattern to extract number from filename
pattern = re.compile(r'\[(\d+)\]\.pickle$')

for root, dirs, files in os.walk(directory_path):
    for file in files:
        if file.endswith('.pickle'):
            full_path = os.path.join(root, file)
            pickle_files.append(full_path)
            
            # Extract number from filename
            match = pattern.search(file)
            if match:
                ids.append(int(match.group(1)))
            else:
                ids.append(None)  # or raise an exception if necessary

# Step 3: Create a new DataFrame with the required structure
# Repeat the single row values for the original columns
data = {col: [single_row[i]] * len(pickle_files) for i, col in enumerate(header)}

# Add the id and mol_feat_file columns
data['id'] = ids
data['mol_feat_file'] = pickle_files

# Create the DataFrame
df = pd.DataFrame(data)

# Reorder columns to make 'id' the first column and 'mol_feat_file' the last column
column_order = ['id'] + header + ['mol_feat_file']
df = df[column_order]

# Save the DataFrame to a new CSV file
df.to_csv(output_csv_path, index=False)