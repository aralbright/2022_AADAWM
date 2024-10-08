#!/bin/bash

# Source directories
source_dir1="/Volumes/albright_postdoc/tubulinRNAseq/trimmed"
source_dir2="/Volumes/albright_postdoc/tubulinRNAseq/trimmed_20230809"

# Destination directory for concatenated files
output_dir="/Volumes/albright_postdoc/tubulinRNAseq/merge"

# Loop through files in the first source directory
for file1 in "$source_dir1"/*.fastq.gz; do
    # Extract the filename without the path
    filename=$(basename "$file1")
    
    # Form the corresponding filename in the second source directory
    file2="$source_dir2/$filename"
    
    # Check if the corresponding file exists in the second source directory
    if [ -f "$file2" ]; then
        # Concatenate the matching files and save to the output directory
        cat "$file1" "$file2" > "$output_dir/$filename"
        echo "Concatenated $filename"
    else
        echo "No matching file found for $filename in $source_dir2"
    fi
done