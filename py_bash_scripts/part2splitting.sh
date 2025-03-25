#!/bin/bash
# Usage: ./part2splitting.sh <fasta_file> <number_of_parts>
# Run from the base dir i.e. above bins/ from part1chunking.sh
# File read in below was calculated from y = 4.4562e-07x² + 7.0125e-04x + 1.3098e-03 
# where x is avg length of seq per chunk (+2 in case of 0_99 as curve overestimates speed here)
# while IFS=, read -r file num_splits; do
#    time /lisc/project/dome/protein_embeddings/py_bash_scripts/part2splitting.sh "bins/$file" "$num_splits"
# done < stats/files_by_how_many_to_split_excluding_1001AAmin.csv
# Old way below that doesn't account for variation in #of seqs per chunk
# for fasta in bins/head500.part_00*_0_99.fasta; do time /lisc/project/dome/protein_embeddings/py_bash_scripts/part2splitting.sh "$fasta" 2 ; done
set -euo pipefail

# Check that exactly two arguments are provided.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <number_of_parts>"
    exit 1
fi

INFILE="$1" # assumed to be in output/bins
NUM_PARTS="$2"
# Determine base directory
BINS_DIR=$(dirname "$INFILE")
BASE_OUTPUT=$(dirname "$BINS_DIR")
PREFIX=$(basename "$INFILE" .fasta)

SPLITS_DIR="${BASE_OUTPUT}/splits"
STATS_DIR="${BASE_OUTPUT}/stats"

mkdir -p "$SPLITS_DIR" "$STATS_DIR"

# Split the FASTA file into the specified number of parts.
seqkit split -p "$NUM_PARTS" -O "$SPLITS_DIR" "$INFILE"

# Define the stats file for this bin (appending stats for each split)
STAT_FILE="${STATS_DIR}/${PREFIX}_stats.txt"

# echo "INFILE=${INFILE}"
# echo "BINS_DIR=${BINS_DIR}"
# echo "BASE_OUTPUT=${BASE_OUTPUT}"
# echo "PREFIX=${PREFIX}"

# -b prints only the basename
seqkit stats -b "$SPLITS_DIR/${PREFIX}"*.fasta >> "$STAT_FILE" 

