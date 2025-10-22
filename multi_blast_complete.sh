#!/bin/bash
   # blast_multi_syntany_improved.sh
   # Description: Master pipeline for multi-genome BLAST synteny analysis
   # Author: Anna Motnenko.
   # Version: 1.0
   # Dependencies: BLAST+, Python 3, GNU utils

# --- HELP MESSAGE ---
show_help() {
    cat << EOF
Usage: $0 [STEP] [DRY_RUN] [--blast-only]

This script performs BLAST searches and retrieves surrounding genes.

Arguments:
  STEP        Optional. Step to start from. Default is 'all'.
              Valid steps: make_databases, run_blast, extract_sequences,
                           retrieve_context_genes, collect_locus_tags, run_python
              You can also use 'all' to run the entire pipeline.
  DRY_RUN     Optional. 'true' to perform a dry run (commands are printed but not executed).
  --blast-only Optional flag to run only up to BLAST and sequence extraction.

Examples:
  Run the full pipeline:
      $0 all

  Run only BLAST + sequence extraction:
      $0 --blast-only

  Dry run of BLAST-only:
      $0 --blast-only true

  Show help:
      $0 -h
      $0 --help
EOF
}

############################################
# PARSE FLAGS
############################################
STEP="all"    # default step
DRY_RUN=false
BLAST_ONLY=false

# Temporary array for remaining positional args
ARGS=()

for arg in "$@"; do
    case "$arg" in
        -h|--help)
            show_help
            exit 0
            ;;
        --blast-only)
            BLAST_ONLY=true
            ;;
        true|false)
            DRY_RUN="$arg"
            ;;
        *)
            STEP="$arg"
            ;;
    esac
done

set -euo pipefail

############################################
# SUPPORTING SCRIPTS
############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

############################################
# GLOBALS
############################################
declare -A summary  # associative array to store status for each query/db combination

############################################
# FUNCTIONS
############################################

log_status() {
    local query="$1"
    local db="$2"
    local step="$3"
    local status="$4"
    summary["$query|$db|$step"]="$status"
}

check_tools() {
    echo "Checking required tools..."
    for tool in makeblastdb blastp blastdbcmd python3; do
        if ! command -v "$tool" &>/dev/null; then
            echo "Error: $tool not found in PATH."
            exit 1
        fi
    done
    echo "All required tools are available."
}

get_inputs() {
    echo "Path to .faa genomes directory:"
    read input_dir
    echo "Path to query .fasta files:"
    read input_query
    echo "Run name (blank = timestamp):"
    read run_name
    if [ -z "$run_name" ]; then
        run_name=$(date +"%Y%m%d_%H%M%S")
    fi
    echo "How many surrounding genes? (Will take # up and downstream):"
    read context
    export input_dir input_query context run_name
}

validate_inputs() {
    [[ -d "$input_dir" && "$(ls -A "$input_dir"/*.faa 2>/dev/null)" ]] || { echo "Error: input_dir invalid"; exit 1; }
    [[ -d "$input_query" && "$(ls -A "$input_query"/*.fasta 2>/dev/null)" ]] || { echo "Error: input_query invalid"; exit 1; }
}

setup_output() {
    RUN_DIR="$input_dir/output_${run_name}"
    DB_DIR="$RUN_DIR/database"
    RESULTS_DIR="$RUN_DIR/results"
    mkdir -p "$DB_DIR" "$RESULTS_DIR"
    LOGFILE="$RUN_DIR/run.log"
    exec > >(tee -a "$LOGFILE") 2>&1
    echo "Run directory: $RUN_DIR"
    echo "Log file: $LOGFILE"
    export RUN_DIR DB_DIR RESULTS_DIR LOGFILE
}

make_databases() {
    echo "STEP 1: Building BLAST databases..."
    database_list="$DB_DIR/database_list.txt"
    > "$database_list"
    for f in "$input_dir"/*.faa; do
        filename=$(basename "${f%.*}")
        db_path="$DB_DIR/${filename}_db/${filename}_db"
        if [[ -f "$db_path.psq" ]]; then
            echo "Database $filename exists, skipping."
            log_status "ALL" "$filename" "db" "skipped"
        else
            mkdir -p "$(dirname "$db_path")"
            makeblastdb -in "$f" -dbtype prot -out "$db_path" -parse_seqids
            echo "Database $filename created."
            log_status "ALL" "$filename" "db" "done"
        fi
        echo "$db_path" >> "$database_list"
    done
    export database_list
}

run_blast() {
    echo "STEP 2: Running BLAST queries (multi-copy aware)..."
    for query_file in "$input_query"/*.fasta; do
        query_base=$(basename "${query_file%.*}")
        while read -r db_path; do
            db_name=$(basename "$db_path")
            out_file="$RESULTS_DIR/${db_name}_vs_${query_base}.txt"

            if [[ -f "$out_file" ]]; then
                echo "BLAST $query_base vs $db_name exists, skipping."
                log_status "$query_base" "$db_name" "blast" "skipped"
            else
                # Remove -max_target_seqs 1 and -culling_limit 1 to keep all hits
                blastp -query "$query_file" -db "$db_path" \
                       -out "$out_file" -outfmt '6 sseqid evalue' \
                       -max_hsps 1 -num_threads 4 -evalue 1e-10

                echo "BLAST $query_base vs $db_name done."
                log_status "$query_base" "$db_name" "blast" "done"
            fi
        done < "$database_list"
    done
}

extract_sequences() {
    echo "STEP 3: Extracting sequences..."
    for f in "$RESULTS_DIR"/*.txt; do
        filename=$(basename "${f%.*}")
        fasta_out="$RESULTS_DIR/${filename}.fasta"
        if [[ -f "$fasta_out" ]]; then
            echo "Sequences $filename exist, skipping."
            log_status "$filename" "ALL" "extract" "skipped"
            continue
        fi
        while read -r db_path; do
            db_name=$(basename "$db_path")
            if [[ "$filename" == "$db_name"* ]]; then
                blastdbcmd -db "$db_path" -entry_batch "$f" -out "$fasta_out" -dbtype prot
                echo "Extracted sequences $filename"
                log_status "$filename" "$db_name" "extract" "done"
            fi
        done < "$database_list"
    done
}
retrieve_context_genes() {
    echo "STEP 4: Retrieving context genes (multi-copy aware)..."

    # associative array to track how many times each locus tag has been processed
    declare -A tag_counts

    # loop over all extracted BLAST sequence files
    for blast_result in "$RESULTS_DIR"/*.fasta; do
        [[ "$blast_result" == *_cluster.faa ]] && continue  # skip already processed files

        base_name=$(basename "${blast_result%.*}")

        # loop over all genome files
        for genome_file in "$input_dir"/*.faa; do
            genome_base=$(basename "${genome_file%.*}")

            # only process if the BLAST result matches this genome
            if [[ "$base_name" == ${genome_base}_db* ]]; then
                echo "Processing BLAST result: $blast_result  against genome: $genome_file"

                # read all locus tags from the BLAST result into a bash array
                mapfile -t locus_tags < <(grep "^>" "$blast_result" | cut -d' ' -f1 | tr -d '>')

                # process each locus tag
                for locus_tag in "${locus_tags[@]}"; do
                    # initialize or increment count safely
                    if [[ -z "${tag_counts[$locus_tag]+_}" ]]; then
                        tag_counts[$locus_tag]=1
                    else
                        tag_counts[$locus_tag]=$(( tag_counts[$locus_tag] + 1 ))
                    fi
                    copy_num=${tag_counts[$locus_tag]}

                    # create a unique output filename for multiple copies
                    if (( copy_num > 1 )); then
                        annotation_output="$RESULTS_DIR/${base_name}_${locus_tag}_copy${copy_num}_cluster.faa"
                    else
                        annotation_output="$RESULTS_DIR/${base_name}_${locus_tag}_cluster.faa"
                    fi

                    # extract surrounding genes using awk
                    awk -v target="$locus_tag" -v range="$context" '
                        /^>/ { headers[++h]=NR; tags[h]=substr($1,2) }
                        { lines[NR]=$0 }
                        END {
                            for (i=1; i<=h; i++) {
                                if (tags[i]==target) {
                                    start=(i-range>1)?headers[i-range]:1
                                    end=(i+range<=h)?headers[i+range]:NR+1
                                    for (j=start; j<end; j++) print lines[j]
                                }
                            }
                        }
                    ' "$genome_file" >> "$annotation_output"

                    echo "Context for $locus_tag (copy $copy_num) written to $annotation_output"
                done
            fi
        done
    done
}


collect_locus_tags() {
    echo "STEP 5: Collecting locus tags (including multi-copy clusters)..."
    tag_output="$RESULTS_DIR/all_locus_tags.txt"
    > "$tag_output"

    total_tags=0
    cluster_count=0

    # Loop over all _cluster.faa files, including *_copyN_cluster.faa
    for tag_file in "$RESULTS_DIR"/*_cluster.faa; do
        [[ ! -s "$tag_file" ]] && continue  # skip empty files
        cluster_count=$((cluster_count + 1))
        tag_count=$(grep -c "^>" "$tag_file" || true)
        total_tags=$((total_tags + tag_count))
        grep "^>" "$tag_file" | cut -d' ' -f1 | tr -d '>' >> "$tag_output"
    done

    echo "✅ Processed $cluster_count cluster files, collected $total_tags locus tags total."
    echo "Locus tags written to: $tag_output"
}
run_python() {
    echo "STEP 6: Running Python script for gene coordinates..."
    tag_output="$RESULTS_DIR/all_locus_tags.txt"
    python3 "$SCRIPT_DIR/gene_coordinates.py" "$tag_output" "$input_dir" "$RESULTS_DIR"
    echo "Python script completed."
}

print_summary() {
    echo
    echo "######## PIPELINE SUMMARY ########"
    printf "%-25s %-25s %-10s %-10s\n" "QUERY" "DB" "STEP" "STATUS"
    for key in "${!summary[@]}"; do
        IFS="|" read -r q db s <<< "$key"
        printf "%-25s %-25s %-10s %-10s\n" "$q" "$db" "$s" "${summary[$key]}"
    done
    echo "##################################"
}

############################################
# MAIN
############################################
main() {
    local step="$1"
    local dry_run="$2"
    local blast_only="$3"

    check_tools
    get_inputs
    validate_inputs
    setup_output

    # Full pipeline steps
    local steps=(make_databases run_blast extract_sequences retrieve_context_genes collect_locus_tags run_python)

    # If --blast-only is set, keep only first three steps
    if [[ "$blast_only" == true ]]; then
        steps=(make_databases run_blast extract_sequences)
    fi

    # Map for resume logic
    declare -A step_map
    for i in "${!steps[@]}"; do
        step_map[${steps[i]}]=$i
    done

    # Determine starting index
    local start_idx=0
    if [[ "$step" != "all" ]]; then
        if [[ -n "${step_map[$step]+_}" ]]; then
            start_idx=${step_map[$step]}
        else
            echo "Unknown step: $step"
            echo "Valid steps: all, ${!step_map[@]}"
            exit 1
        fi
    fi

    echo "Dry run: $dry_run"
    echo "Steps to execute (resume from $step): ${steps[@]:start_idx}"

    # Execute steps
    for ((i=start_idx; i<${#steps[@]}; i++)); do
        if [[ "$dry_run" == "true" ]]; then
            echo "[DRY RUN] Would execute step: ${steps[i]}"
        else
            "${steps[i]}"
        fi
    done

    print_summary
    echo "Pipeline finished. Results in: $RUN_DIR"
}

main "$STEP" "$DRY_RUN" "$BLAST_ONLY"

