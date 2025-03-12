#!/bin/bash
#SBATCH --job-name=merge_checkpoints
#SBATCH --output=merge_checkpoints_%A.out
#SBATCH --error=merge_checkpoints_%A.err
#SBATCH --mail-type=END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000M
#SBATCH --time=00:20:00
#SBATCH --partition=basic,short

# Directory where individual checkpoint files are stored
CHECKPOINT_DIR="/lisc/scratch/dome/pullen/GlobDB"
MASTER_FILE="${CHECKPOINT_DIR}/all_embeddings.h5"

# ***********************
CHUNK="200_299" # change here
# ***********************
# CHUNK="100_199" # change here

# List all individual checkpoint files matching your pattern.
# Adjust the pattern as needed.
temp_files=$(find "${CHECKPOINT_DIR}" -maxdepth 1 -name "embeddings_${CHUNK}_*.h5")

echo "Merging checkpoint files:"
echo "${temp_files}"

# Create (or open) the master file and merge datasets from individual files.
python <<EOF
import h5py, glob, os, fcntl
from pathlib import Path

master_path = "${MASTER_FILE}"
checkpoint_dir = "${CHECKPOINT_DIR}"
# Get list of temp checkpoint files.
temp_files = glob.glob(os.path.join("${CHECKPOINT_DIR}", "embeddings_${CHUNK}_*.h5"))

try:
    master_file = Path(master_path)
    # Open or create master file in append mode to obtain a descriptor for locking
    with open(master_file, 'ab+') as f:
        print("Acquiring exclusive lock...")
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            # We open it a second time, now with h5py to interact with the HDF5 file
            with h5py.File(f, "a") as master_hf:
                pre_merge = len(master_hf.keys())
                print("Pre-merging master file keys count:", pre_merge)
                total_added = 0
                total_skipped = 0
                for temp_file in temp_files:
                    print("Merging", temp_file)
                    with h5py.File(temp_file, "r") as temp_hf:
                        for key in temp_hf.keys():
                            if key not in master_hf:
                                # Use h5py.h5o.copy to copy the dataset directly without loading into memory.
                                h5py.h5o.copy(temp_hf.id, key.encode("utf-8"), master_hf.id, key.encode("utf-8"))
                                total_added += 1
                            else:
                                print(f"Protein {key} already exists in master; skipping.")
                                total_skipped += 1
                post_merge = len(master_hf.keys())
                print("Total proteins added:", total_added)
                print("Total proteins skipped:", total_skipped)
                print("Post-merging master file keys count:", post_merge)
        finally:
            print("Releasing exclusive lock...")
            fcntl.flock(f, fcntl.LOCK_UN)
except Exception as e:
    print(f"Error merging {temp_file} into master embeddings file - {e}")
    raise
EOF

echo "Merge complete. Master file: ${MASTER_FILE}"

