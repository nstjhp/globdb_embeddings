#!/bin/bash
# Usage: ./part2splitting.sh <fasta_file> <number_of_parts>

# Check that exactly two arguments are provided.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <number_of_parts>"
    exit 1
fi

INFILE="$1"
NUM_PARTS="$2"

# Split the FASTA file into the specified number of parts.
seqkit split -p "$NUM_PARTS" "$INFILE"

# The split files are saved in a directory named "$INFILE.split"
SPLIT_DIR="${INFILE}.split"

# Check if the split directory exists.
if [ ! -d "$SPLIT_DIR" ]; then
    echo "Error: Split directory $SPLIT_DIR not found."
    exit 1
fi

# Run seqkit stats on each FASTA file in the split directory.
for file in "$SPLIT_DIR"/*.fasta; do
    echo "Stats for $file:"
    seqkit stats "$file"
    echo "-------------------------------------"
done
