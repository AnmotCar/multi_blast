# Multi-BLAST Pipeline with Gene Neighbourhood

A Bash-based workflow for running multi-genome BLAST comparisons and extracting regions around query genes.

## Overview ##

This pipeline automates:
- Database creation for all `.faa` or `.fna` files
- BLAST of multiple query `.fasta` sequences
- Extraction of surrounding genes (customizable)
- Retrieval of locus tags and contextual regions
- Creation of table with data summury as a `.csv`
- BLASTp is default, can switch to BLASTn


## Requirements ##

- BLAST+ tools (`makeblastdb`, `blastp`, `blastn`)
- Biopyton (if using the final step)
- GNU core utilities (`awk`, `sed`, etc.)
	-bioawk used here


## File structure layout ##
- protein containing .faa or nucleotide containing .fna files and GenBank (.gbk) need to be in the same directory and have matching names
- Script optimized to work with prokka output
- Since prokka output is required an included script "rename_pull.sh" is included -> will copy required files from a prokka output and name based on folder name

- Results folder will contain
	- genome1_vs_quiry.fasta -> the blast hit 
	- genome1_vs_quiry.txt -> locus lag and e-value
	- genome1_vs_quiry_cluster.fasta -> if using the neighbourhood search, sequences of all hits
	- all_locus_tags.txt -> list of all locus tags in the analysis
	- gene_coordinates.csv -> file containing hits with "Locus_Tag/Protein_ID/Contig/Start/End/Strand/Gene/Product/Organism/Source_File"

Directory structure example:

<pre> ```
~/path/you/have/ 
|    └── genomes/ 
|         ├── genome1.faa 
|         ├── genome1.fna 
|         ├── genome1.gbk 
|         ├── genome2.faa
|         ├── genome2.faa 
|         ├── genome2.gbk 
|         └── output_runname/ 
|           ├── database 
|           ├── results 
|           └── run.log 
└── input/ 
 └── gene1-3.fasta
``` </pre>
   
## Usage ##

Run the main script:

---
bash /path/to/script/multi_blast_complete.sh

---

Usage: /path/to/script/multi_blast_complete.sh [STEP] [DRY_RUN] [--blast-only] [--blastn]

When script is run prompts will be
- Path to .faa or .fna genomes directory:
- Path to query .fasta files:
- Run name (blank = timestamp):
- How many surrounding genes? (Will take # up and downstream): note if running blast mode only just hit enter and ignore


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
      /path/to/script/multi_blast_complete.sh

  Run only BLAST + sequence extraction:
      /path/to/script/multi_blast_complete.sh --blast-only
    
  Run only nucleotide BLAST:
      /path/to/script/multi_blast_complete.sh --blastn

  Dry run of BLAST-only:
      /path/to/script/multi_blast_complete.sh --blast-only true

  Show help:
      /path/to/script/multi_blast_complete.sh -h

Note:
- This script will full all hits matching the blast quiry until the e-value cut off (1e-10) modify script on line 194 to change.
- Can also make this script only return the first best hit by adding `-max_target_seqs 1 -max_hsps 1` on line 194 as well. 

