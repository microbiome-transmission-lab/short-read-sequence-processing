#!/bin/bash

# Check if a directory is provided as an argument
if [ $# -eq 0 ]; then
  echo "Error: Please provide a directory path"
  exit 1
fi

# Store the provided directory path
input_dir=$1
output_dir=$2
koslicki_ko_sketch=$3
kmer_size=$4
scale_size=$5

# Specify the name of your Python script
python_script="$HOME/funprofiler/funcprofiler.py"  # Replace with the actual name of your script

# Note: Zenodo source for FracMiniHash sketches of Orthologs in KEGG database from the KoslickiLab: https://zenodo.org/records/10045253
# koslicki_ko_sketch="$HOME/funprofiler/demo/KOs_sketched_scaled_1000.sig.zip"

# Loop through all files in the directory
for file in "$input_dir"/*; do
  if [ -f "$file" ]; then  # Make sure it's a regular file
    filename=$(basename "$file")
    output_file_ko="${filename%.*}_funprofiler_ko_profiles.csv"  # Example output naming
    output_file_gather="${filename%.*}_funprofiler_gather_out.csv" # Second output name

    # Run the Python script, passing input and output filenames
    python "$python_script" "$file" "$koslicki_ko_sketch" "$kmer_size" "$scale_size" "$output_dir/$output_file_ko" -g "$output_dir/$output_file_gather"
  fi
done


