import os
import pandas as pd

def check_unique_ids(directory):
    all_ids = []
    
    # Loop through all CSV files in the directory
    for file_name in sorted(os.listdir(directory)):
        if file_name.endswith(".csv"):
            file_path = os.path.join(directory, file_name)
            df = pd.read_csv(file_path, usecols=["id"])
            
            # Append all ids to the list
            all_ids.extend(df["id"].tolist())
    
    # Convert list to a set to find unique ids
    unique_ids = set(all_ids)
    
    # Check for duplicates
    if len(all_ids) == len(unique_ids):
        print("All 'id' entries are unique.")
    else:
        print("There are duplicate 'id' entries.")
    
    # Print the number of unique ids
    print(f"Number of unique ids: {len(unique_ids)}")

if __name__ == "__main__":
    # Specify the directory containing the .csv files
    directory = "/home/cm2231/rds/project/rds-mphil-rkVIP0ZCj0k/carolina/test_Medin/split_zinc"
    
    # Check for unique ids
    check_unique_ids(directory)
