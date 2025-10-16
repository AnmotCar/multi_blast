# BLAST Multi-Synteny Pipeline

A Bash-based workflow for running multi-genome BLAST comparisons and extracting regions around query genes.

## Overview ##

This pipeline automates:
- Database creation for all `.faa` files
- BLAST of multiple query `.fasta` sequences
- Extraction of surrounding genes (customizable)
- Retrieval of locus tags and contextual regions
- Optional Python post-processing


## Requirements ##

- BLAST+ tools (`makeblastdb`, `blastp`)
- Biopyton (if using the final step)
- GNU core utilities (`awk`, `sed`, etc.)
	-bioawk used here


## File structure layout ##
- protein containing .faa file and GenBank (.gbk) need to be in the same directory and have matching names
- Script optimized to work with prokka output
- Since prokka output is required an included script "rename_pull.sh" is included -> will copy required files from a prokka output and name based on folder name

- Results folder will contain
	- genome1_vs_quiry.fasta -> the blast hit 
	- genome1_vs_quiry.txt -> locus lag and e-value
	- genome1_vs_quiry_cluster.fasta -> if using the neighbourhood search, sequences of all hits
	- all_locus_tags.txt -> list of all locus tags in the analysis
	- gene_coordinates.csv -> file containing hits with "Locus_Tag/Protein_ID/Start/End/Strand/Gene/Product/Organism/Source_File"

Directory structure example:
"
~/path/you/have/ 
|    ├── genomes/ 
|    ├── genome1.faa 
|    ├── genome1.gbk 
|    ├── genome2.faa 
|    ├── genome1.gbk 
|    └── output_runname/ 
|           ├── database 
|           ├── results 
|           └── run.log 
└── input/ 
 └── gene1-3.fasta
"
   
## Usage ##

Run the main script:

---
bash /path/to/script/multi_blast_complete.sh

---

Usage: /path/to/script/multi_blast_complete.sh [STEP] [DRY_RUN] [--blast-only]

This script performs BLAST searches and retrieves surrounding genes.

Arguments:
STEP        Optional. Step to start from. Default is 'all'.
              Valid steps: make_databases, run_blast, extract_sequences,
                           retrieve_context_genes, collect_locus_tags, run_python
              You can also use 'all' to run the entire pipeline.
  
DRY_RUN     Optional. 'true' to perform a dry run (commands are printed but not executed).

--blast-only Optional flag to run only up to BLAST and sequence extraction skipping neighbourhood search.

Examples:
  Run the full pipeline:
      /path/to/script/multi_blast_complete.sh all

  Run only BLAST + sequence extraction:
      /path/to/script/multi_blast_complete.sh --blast-only

  Dry run of BLAST-only:
      /path/to/script/multi_blast_complete.sh --blast-only true

  Show help:
      /path/to/script/multi_blast_complete.sh -h



