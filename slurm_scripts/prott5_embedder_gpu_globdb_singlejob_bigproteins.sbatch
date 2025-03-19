#!/bin/bash
#SBATCH --job-name=prott5_embed_redo_fails
# #SBATCH --output=job%A_GPU.out
# #SBATCH --error=job%A_GPU.err
#SBATCH --mail-type=END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU
#SBATCH --mem=8000M
#SBATCH --time=0-02:00:00
#SBATCH --partition=gpu # needs to use a GPU with more RAM for big proteins
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:l40s:1

# Exit the slurm script if a command fails
set -e

# ***********************
SPLIT_ID="REDO_FAILS"
# ***********************

# Construct filenames dynamically
ID_STRING="${SLURM_JOB_ID}_${SLURMD_NODENAME}_GPU"
SLURM_OUTPUT_FILENAME="${ID_STRING}.out"
SLURM_ERROR_FILENAME="${ID_STRING}.err"
SLURM_OUTPUT_FILEPATH="/lisc/scratch/dome/pullen/GlobDB/${SLURM_OUTPUT_FILENAME}"
SLURM_ERROR_FILEPATH="/lisc/scratch/dome/pullen/GlobDB/${SLURM_ERROR_FILENAME}"

# Redirect output and error to my generated filenames
exec 1>${SLURM_OUTPUT_FILEPATH} 2>${SLURM_ERROR_FILEPATH}

# Test these job commands
echo "Job ID: ${SLURM_JOB_ID}"
echo "TMPDIR: ${TMPDIR}"
echo "Node: ${SLURMD_NODENAME}"

# Try the unsplit >=1001 AAs fasta and hope that the previous embeddings can be skipped appropriately 
# saves doing the python in a loop over every split
FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/filtered1001AAmin_sorted.fasta"
# FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/filtered1001AAmin_sorted.fasta.split/filtered1001AAmin_sorted.part_${SPLIT_ID}.fasta"

seqkit stats -b ${FASTA_INPUT_FILE}

MAX_RESIDUES=4000
MAX_SEQ_LEN=1000
MAX_BATCH=200

# For use in the Python logging, here we set MY_SLURM_PROCESS_ID
# and base it on whether this is an array job (-n checks if a string is non-empty)
if [[ -n "${SLURM_ARRAY_JOB_ID}" && -n "${SLURM_ARRAY_TASK_ID}" ]]; then
    # This is an array job
    export MY_SLURM_PROCESS_ID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    # This is not an array job; use SLURM_JOB_ID instead
    export MY_SLURM_PROCESS_ID="${SLURM_JOB_ID}"
fi
echo "MY_SLURM_PROCESS_ID: ${MY_SLURM_PROCESS_ID}"

# Run the Python script
python /lisc/project/dome/protein_embeddings/py_bash_scripts/prott5_embedder_nick.py \
  --input ${FASTA_INPUT_FILE} \
  --output $TMPDIR/embeddings_1001AAmin_${SPLIT_ID}_${SLURM_JOB_ID}.h5 \
  --log $TMPDIR/${SLURM_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log \
  --max_residues ${MAX_RESIDUES} --max_seq_len ${MAX_SEQ_LEN} --max_batch ${MAX_BATCH} \
  --master_embedding_file /lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/embeddings/master_embeddings.h5

if [ -f "$TMPDIR/embeddings_1001AAmin_${SPLIT_ID}_${SLURM_JOB_ID}.h5" ]; then
    cp "$TMPDIR/embeddings_1001AAmin_${SPLIT_ID}_${SLURM_JOB_ID}.h5" /lisc/scratch/dome/pullen/GlobDB/
else
    echo "File $TMPDIR/embeddings_1001AAmin_${SPLIT_ID}_${SLURM_JOB_ID}.h5 not found; skipping copy."
fi

if [ -f "$TMPDIR/${SLURM_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log" ]; then
    cp "$TMPDIR/${SLURM_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log" /lisc/scratch/dome/pullen/GlobDB/
else
    echo "File $TMPDIR/${SLURM_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log not found; skipping copy."
fi

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

