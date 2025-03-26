#!/bin/bash

OUTPUT_FILE="/lisc/scratch/dome/pullen/GlobDB/linclust/valid_tasks4slurm.txt"
NUM_SPLITS_FILE="/lisc/scratch/dome/pullen/GlobDB/linclust/stats/files_by_how_many_to_split_excluding_1001AAmin.csv"
# Find max SPLIT_ID from NUM_SPLITS_FILE (second column has split counts)
# NR>1 means ignore the first line i.e. the header
MAX_SPLIT_ID=$(awk -F',' 'NR>1 {if($2>max){max=$2}} END {print max}' "$NUM_SPLITS_FILE")
echo "MAX_SPLIT_ID = ${MAX_SPLIT_ID}"
# Clear previous output file
> ${OUTPUT_FILE}

# Define chunks in Bash array
CHUNKS=("0_99" "100_199" "200_299" "300_399" "400_499" \
        "500_599" "600_699" "700_799" "800_899" "900_1000")

# Iterate through CHUNKS, SPLIT_IDs (zero-padded), and PART_IDs (also zero-padded)
for PART in {1..10}; do
    PART_PADDED=$(printf "%03d" "$PART")  # Zero-pad to 3 digits
    for CHUNK in "${CHUNKS[@]}"; do
        for (( SPLIT_ID=1; SPLIT_ID<=MAX_SPLIT_ID; SPLIT_ID++ )); do
            SPLIT_ID_PADDED=$(printf "%03d" "$SPLIT_ID")  # Zero-pad to 3 digits
            FILE="/lisc/scratch/dome/pullen/GlobDB/linclust/splits/clusters_more_than1.part_${PART_PADDED}_filtered1000AAmax_sorted_${CHUNK}.part_${SPLIT_ID_PADDED}.fasta"
            if [ -f "${FILE}" ]; then
                echo "${PART_PADDED},${CHUNK},${SPLIT_ID_PADDED}" >> ${OUTPUT_FILE}
            fi
        done
    done
done

