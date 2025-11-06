#!/bin/bash

# Pull files from prokka output into a directory needed for BLAST
# Better version of copy_gbk_faa.sh

input_dir=$1

if [[ -z "$input_dir" ]]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

mkdir -p "${input_dir}/genomes"

shopt -s nullglob
for f in "$input_dir"/*; do
    if [[ -d "$f" && "$(basename "$f")" != "genomes" ]]; then
        folder_name=$(basename "$f")
        echo "Processing: $folder_name"

        for file in "$f"/*.ffn "$f"/*.faa "$f"/*.gbk ; do # change target file extension to match need
            if [[ -f "$file" ]]; then
                ext="${file##*.}"  # Get file extension
                new_name="${folder_name}.${ext}"
                cp "$file" "${input_dir}/genomes/${new_name}"
                echo "Copied and renamed: $file -> ${new_name}"
            fi
        done
    fi
done
shopt -u nullglob

#for file in "$f"/*.gbff "$f"/*.faa "$f"/*.ffn "$f"/*.gbk; do
#             genbank,    protein,     nucleotide,    other genbank