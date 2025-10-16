from Bio import SeqIO
import sys
import glob
import os

# Same as gene_cooridnates_from_blast.py but works with prokka annotated files (becuase they are different from PGAP)

# Get command-line arguments
locus_tags_file = sys.argv[1]  # Path to the file containing locus tags
input_dir = sys.argv[2]        # Path to the directory containing GenBank files

# Optional third argument: custom output directory
if len(sys.argv) > 3:
    output_dir = sys.argv[3]
else:
    output_dir = input_dir  # default: write directly into input_dir

# Define output file path
output_file = os.path.join(output_dir, "gene_coordinates.tsv")

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Read the locus tags
with open(locus_tags_file, "r") as f:
    locus_tags = set(line.strip() for line in f)

# Debugging: Check if locus tags were read correctly
print(f"Locus tags loaded ({len(locus_tags)} total): {locus_tags}")

# Check if any GenBank files exist in the directory
gbk_files = glob.glob(os.path.join(input_dir, "*.gbk"))
print(f"GenBank files found ({len(gbk_files)} total): {gbk_files}")

if not gbk_files:
    print("No GenBank files found! Check your input directory.")
    sys.exit(1)

# Open the output file once for writing
with open(output_file, "w") as out:
    # Write the header only once
    out.write("Locus_Tag\tProtein_ID\tStart\tEnd\tStrand\tGene\tProduct\tOrganism\tSource_File\n")

    # Iterate over all .gbk files in the input directory
    for genbank_file in gbk_files:
        print(f"\nProcessing file: {genbank_file}")
        base_name = os.path.basename(genbank_file)  # Get the file name

        # Parse the GenBank file and extract gene information
        with open(genbank_file, "r") as gbk:
            for record in SeqIO.parse(gbk, "genbank"):
                organism = record.annotations.get("organism", "Unknown_Organism")  # Extract organism name
                print(f"Processing organism: {organism}")

                # Loop through features in the GenBank file
                for feature in record.features:
                    if feature.type == "CDS":  # Look for CDS annotations
                        if "locus_tag" in feature.qualifiers:
                            locus_tag = feature.qualifiers["locus_tag"][0]
                            print(f"Found locus_tag: {locus_tag}")

                            if locus_tag in locus_tags:
                                print(f"Matching locus_tag found: {locus_tag}")

                                protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                                start = feature.location.start + 1  # Convert 0-based to 1-based indexing
                                end = feature.location.end
                                strand = "+" if feature.location.strand == 1 else "-"
                                gene = feature.qualifiers.get("gene", ["N/A"])[0]
                                product = feature.qualifiers.get("product", ["N/A"])[0]

                                # Write extracted data to the output file
                                out.write(f"{locus_tag}\t{protein_id}\t{start}\t{end}\t{strand}\t{gene}\t{product}\t{organism}\t{base_name}\n")
                                print(f"Written to output: {locus_tag}, {start}-{end}, {strand}")

print(f"\nResults written to: {output_file}")
