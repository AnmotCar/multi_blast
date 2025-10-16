#!/bin/bash

# --- HELP MESSAGE ---
show_help() {
    cat << EOF
Copy files from prokka output into a directory needed for BLAST (genomes)

Usage: <input_directory>

Arguments:
  input_directory  ->  Directory containing subdirectories with Prokka output files.

Output:
  A 'genomes' directory inside the input_directory with renamed files.

EOF
}

# --- PARSE ARGUMENTS ---
if [[ "$1" == "-h" || "$1" == "--help" || $# -eq 0 ]]; then
    show_help
    exit 0
fi


# Check for correct number of arguments
input_dir=$1

if [[ -z "$input_dir" ]]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

mkdir -p "${input_dir}/genomes"

shopt -s nullglob
for f in "$input_dir"/*; do
    if [[ -d "$f" ]]; then
        folder_name=$(basename "$f")
        echo "Processing: $folder_name"

        for file in "$f"/*.faa "$f"/*.gbk ; do # change target file extension to match need
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
echo "All files copied and renamed in ${input_dir}/genomes"