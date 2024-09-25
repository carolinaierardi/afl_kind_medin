import pandas as pd
import glob
import os
import argparse
from csvsort import csvsort

def main(path_to_processed_csv, path_to_index_files, chunk_size, output_csv_path, top_10000_path):
    # Paths to the files
    index_files = glob.glob(os.path.join(path_to_index_files, "*.csv"))

    # Step 1: Create a dictionary to map ID to NAME and SMILES from INDEX.csv files
    id_info = {}
    for index_file in index_files:
        index_df = pd.read_csv(index_file)  # The file already has headers
        for _, row in index_df.iterrows():
            id_info[str(row['id'])] = {'smiles': row['smiles']}

    # Step 2: Process the large CSV in chunks
    chunks = pd.read_csv(path_to_processed_csv, chunksize=chunk_size)

    with open(output_csv_path, 'w') as f_out:
        # Write the header
        f_out.write("id, pred_affinity, smiles\n")

        for chunk in chunks:
            # Ensure ID is treated as string for accurate mapping
            chunk['id'] = chunk['id'].astype(str)

            # Map the IDs to NAME and SMILES
            chunk['smiles'] = chunk['id'].apply(lambda x: id_info.get(x, {}).get('smiles', ''))

            # Write to the output file without header
            chunk.to_csv(f_out, index=False, header=False)

    # Print the number of rows in the final DataFrame
    csvsort(output_csv_path, [1])
    num_rows = sum(1 for _ in open(output_csv_path)) - 1  # Subtract header row
    top_x = pd.read_csv(output_csv_path, skiprows = num_rows-10000)
    top_x.columns = ["id","pred_affinity","smiles"]
    top_x.sort_values(by = "pred_affinity", ascending = False, inplace = True)
    top_x.to_csv(top_10000_path, index = False)

    print(f"The final DataFrame has {num_rows} rows.")
    print(f"Final DataFrame with IDs, affinities, NAME, and SMILES has been saved to '{output_csv_path}'.")
    print(f"Top 10000 DataFrame with IDs, affinities, NAME, and SMILES has been saved to '{top_10000_path}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CSV files.')
    parser.add_argument('--path_to_processed_csv', type=str, required=True, help='Path to the processed CSV file')
    parser.add_argument('--path_to_index_files', type=str, required=True, help='Path to the index files directory')
    parser.add_argument('--chunk_size', type=int, required=True, help='Chunk size for processing the CSV')
    parser.add_argument('--output_csv_path', type=str, required=True, help='Output path for the final CSV file')
    parser.add_argument('--top_10000_path', type=str, required=True, help='Output path for the top 10000 CSV file')
    args = parser.parse_args()

    main(args.path_to_processed_csv, args.path_to_index_files, args.chunk_size, args.output_csv_path, args.top_10000_path)
