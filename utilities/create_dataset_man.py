import pandas as pd
import os

# Directory containing the files
data_dir = "C:/Users/walth/OneDrive - University of Strathclyde/PROJECTS/QuantumSSA/files/man_dataset/"

# Get a list of all files in the directory
all_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith('.txt')]
all_files = sorted(all_files, key=lambda f: int(''.join(filter(str.isdigit, os.path.basename(f)))))

# Initialize an empty list to store DataFrames
data_frames = []

# Loop through files and load each one into a DataFrame
for file in all_files:
    # Load the file into a DataFrame
    # Assuming the data is space-separated and the header row is missing
    df = pd.read_csv(file, delim_whitespace=True, header=None)
    data_frames.append(df)

# Concatenate all DataFrames into a single DataFrame
full_data = pd.concat(data_frames, ignore_index=True)

# Print the combined DataFrame
print(full_data)