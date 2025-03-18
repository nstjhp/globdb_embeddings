#!/bin/bash
# This script processes a large FASTA file by:
# 1. Filtering out proteins longer than 1000 AAs.
# 2. Creating a tab-delimited table of sequence IDs and lengths.
# 3. Splitting the table into bins (by length range) without doing a full sort.
# 4. Extracting the sequences for each bin into separate FASTA files.
# 5. Running seqkit stats on each split file.
#
# Once you've decided by how many to split these FASTAs, you are ready for part2splitting.sh
# Case 1) In 00:04:47 I processed 6827 proteins that had average length of 50 AAs.
# Case 2) In 00:20:35 I processed 6828 proteins that had average length of 239 AAs.
# Case 3) In 01:33:56 I processed 6828 proteins that had average length of 806 AAs.

# Input big FASTA file
BIG_FASTA="chloroflexi_test1000.faa"

# 1. Filter out proteins longer than 1000 AAs
seqkit seq -M 1000 "$BIG_FASTA" > filtered1000AAmax.fasta

# 2. Create a tab-delimited table (ID and length)
seqkit fx2tab -l -n -i filtered1000AAmax.fasta > lengths.tsv

# 3. Define bins by length range.
# For bins 0-8, we'll have: 0-99, 100-199, ..., 800-899.
# For bin 9, we'll have: 900-1000 (including 1000).
for bin in {0..9}; do
    lower=$(( bin * 100 ))
    if [ $bin -eq 9 ]; then
        upper=1000
    else
        upper=$(( lower + 99 ))
    fi
    bin_ids="bin_${lower}_${upper}.tsv"
    # Filter the lengths.tsv file. Each line is expected to have two columns: ID and length.
    # We assume no leading ">" because we used the -n flag.
    awk -v low="$lower" -v high="$upper" '$2 >= low && $2 <= high {print $1}' lengths.tsv > "$bin_ids"

    bin_fasta="filtered1000AAmax_sorted_${lower}_${upper}.fasta"
    # Extract sequences for the current bin into a FASTA file.
    seqkit grep -f "$bin_ids" filtered1000AAmax.fasta > "$bin_fasta"
done

# 4. Run seqkit stats on each generated FASTA file.
for fasta in filtered1000AAmax_sorted_*.fasta; do
    echo "Stats for $fasta:"
    seqkit stats "$fasta"
    echo "----------------------------------------------------"
done
