from Bio import SeqIO
import sys
import glob
import os

# Command-line arguments
locus_tags_file = sys.argv[1]  # Path to locus tag list
input_dir = sys.argv[2]        # Directory with GenBank files
output_dir = sys.argv[3] if len(sys.argv) > 3 else input_dir

# Define output file path
output_file = os.path.join(output_dir, "gene_coordinates.tsv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Read the locus tags
with open(locus_tags_file, "r") as f:
    locus_tags = set(line.strip() for line in f)

print(f"Locus tags loaded ({len(locus_tags)}): {list(locus_tags)[:5]}...")

# Find all GenBank files
gbk_files = glob.glob(os.path.join(input_dir, "*.gbk"))
if not gbk_files:
    print("No GenBank files found! Check your input directory.")
    sys.exit(1)

# Open output file
with open(output_file, "w") as out:
    out.write("Locus_Tag\tProtein_ID\tContig\tStart\tEnd\tStrand\tGene\tProduct\tOrganism\tSource_File\n")

    for genbank_file in gbk_files:
        print(f"\nProcessing file: {genbank_file}")
        base_name = os.path.basename(genbank_file)

        for record in SeqIO.parse(genbank_file, "genbank"):
            organism = record.annotations.get("organism", "Unknown_Organism")
            contig = record.id  #this is the contig name

            for feature in record.features:
                if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers["locus_tag"][0]

                    if locus_tag in locus_tags:
                        protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]
                        start = int(feature.location.start) + 1  # 1-based index
                        end = int(feature.location.end)
                        strand = "+" if feature.location.strand == 1 else "-"
                        gene = feature.qualifiers.get("gene", ["N/A"])[0]
                        product = feature.qualifiers.get("product", ["N/A"])[0]

                        # Include contig in output
                        out.write(f"{locus_tag}\t{protein_id}\t{contig}\t{start}\t{end}\t{strand}\t{gene}\t{product}\t{organism}\t{base_name}\n")

                        #print(f"✔ {locus_tag} ({contig}) {start}-{end} ({strand})")

print(f"\n Results written to: {output_file}")
