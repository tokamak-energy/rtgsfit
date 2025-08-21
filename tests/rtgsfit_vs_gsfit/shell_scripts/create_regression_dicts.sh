#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Directory containing the data files (relative to the script location)
DATA_DIR="$SCRIPT_DIR/../data"

# Loop over files starting with rtgsfit_output_dict_ in the data directory
for file in "$DATA_DIR"/rtgsfit_output_dict_*.npy; do
    # Skip if no files match
    [ -e "$file" ] || continue

    # Get the filename without the directory path
    filename=$(basename "$file")

    # Create the new filename by inserting '_regression' before '.npy'
    new_filename="${filename%.npy}_regression.npy"

    # Copy the file
    cp "$file" "$DATA_DIR/$new_filename"

    echo "Copied $file to $DATA_DIR/$new_filename"
done
