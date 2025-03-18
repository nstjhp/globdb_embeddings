#!/bin/bash
# This script processes a large FASTA file by:
# 1. Filtering out proteins longer than 1000 AAs (those >= 1001 go into a their own 'bin')
# 2. Creating a tab-delimited table of sequence IDs and lengths.
# 3. Splitting the table into bins (by length range) without doing a full sort.
# 4. Extracting the sequences for each bin into separate FASTA files.
# 5. Running seqkit stats on each split file.
#
# To loop over a dir of fastas
# for fasta in *.fasta; do time /lisc/project/dome/protein_embeddings/py_bash_scripts/part1chunking.sh -i "$fasta" ; done
# 
# Once you've decided by how many to split these FASTAs, you are ready for part2splitting.sh
set -euo pipefail

usage() {
    echo "Usage: $0 -i <input_fasta> [-o <output_directory>]"
    exit 1
}

# Check if seqkit is installed
if ! command -v seqkit >/dev/null 2>&1; then
    echo "Error: seqkit is not installed. Please install seqkit and try again."
    exit 1
fi

# Parse command-line options
OUTPUT_DIR="output"  # default base output directory
while getopts "i:o:" opt; do
    case ${opt} in
        i)
            INPUT_FASTA="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${INPUT_FASTA:-}" ]; then
    echo "Error: Input FASTA file is required."
    usage
fi

# Create a structured output directory
FILTERED_DIR="${OUTPUT_DIR}/filtered"
LENGTHS_DIR="${OUTPUT_DIR}/lengths"
BINS_DIR="${OUTPUT_DIR}/bins"
STATS_DIR="${OUTPUT_DIR}/stats"

mkdir -p "$FILTERED_DIR" "$LENGTHS_DIR" "$BINS_DIR" "$STATS_DIR"

# Get the base name (without path and extension) to serve as a unique prefix
BASE=$(basename "$INPUT_FASTA")
PREFIX="${BASE%.*}"

echo "Processing file: $INPUT_FASTA"
echo "Using prefix: $PREFIX"
echo "Output directory: $OUTPUT_DIR"

# Step 1: Filter sequences based on protein length
OUT_MAX="${FILTERED_DIR}/${PREFIX}_filtered1000AAmax.fasta"
OUT_MIN="${FILTERED_DIR}/${PREFIX}_filtered1001AAmin.fasta"

echo "Filtering proteins with length <= 1000 AAs into $OUT_MAX..."
seqkit seq -M 1000 "$INPUT_FASTA" > "$OUT_MAX"

echo "Filtering proteins with length >= 1001 AAs into $OUT_MIN..."
seqkit seq -m 1001 "$INPUT_FASTA" > "$OUT_MIN"

# Initialise stats file
MAIN_STATS="${STATS_DIR}/${PREFIX}_stats.txt"
#echo "# Main seqkit stats for $PREFIX" > "$MAIN_STATS"

# Step 2: Create a tab-delimited table (ID and length) for proteins <= 1000 AAs
LENGTHS_FILE="${LENGTHS_DIR}/${PREFIX}_lengths.tsv"
echo "Generating lengths table: $LENGTHS_FILE..."
seqkit fx2tab -l -n -i "$OUT_MAX" > "$LENGTHS_FILE"

# Step 3: Bin sequences by length (bins: 0-99, 100-199, …, 900-1000)
for bin in {0..9}; do
    lower=$(( bin * 100 ))
    if [ $bin -eq 9 ]; then
        upper=1000
    else
        upper=$(( lower + 99 ))
    fi
    BIN_IDS="${LENGTHS_DIR}/${PREFIX}_bin_${lower}_${upper}.tsv"
    echo "Creating bin for lengths ${lower}-${upper} into $BIN_IDS..."
    awk -v low="$lower" -v high="$upper" '$2 >= low && $2 <= high {print $1}' "$LENGTHS_FILE" > "$BIN_IDS"

    BIN_FASTA="${BINS_DIR}/${PREFIX}_filtered1000AAmax_sorted_${lower}_${upper}.fasta"
    echo "Extracting sequences for bin ${lower}-${upper} into $BIN_FASTA..."
    seqkit grep -f "$BIN_IDS" "$OUT_MAX" > "$BIN_FASTA"
done

# Step 4: Run seqkit stats on each binned FASTA file and the 1001AAmin fasta
seqkit stats -b "$OUT_MAX" "${BINS_DIR}/${PREFIX}_filtered1000AAmax_sorted_"*.fasta "$OUT_MIN" >> "$MAIN_STATS" 

# Clean up the large intermediate file (filtered1000AAmax) as its data is now in the bins
rm -f "$OUT_MAX"
# Remove binned ID files as they can be recreated from the lengths file
rm -f "${LENGTHS_DIR}/${PREFIX}_bin_"*.tsv

# Move the filtered1001AAmin file to the bins directory (as it represents the >=1001 bin)
mv "$OUT_MIN" "${BINS_DIR}/"
# This is now empty so remove
rmdir "$FILTERED_DIR"

echo "Processing complete. All output files are organized under the directory: $OUTPUT_DIR"

