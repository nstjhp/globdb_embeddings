#!/bin/bash
#SBATCH --job-name=prott5_embed_array
# #SBATCH --output=job%A_%a_GPU.out
# #SBATCH --error=job%A_%a_GPU.err
#SBATCH --mail-type=ARRAY_TASKS,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU and more shards
#SBATCH --mem=16000M
#SBATCH --time=0-02:00:00
#SBATCH --partition=basic#,gpu 
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:t4:1
#SBATCH --exclude=node-c[01-02] # trick to exclude c nodes that for unknown reasons use 2 threads
# #SBATCH --threads-per-core=1 # docs not clear if this is min 1, or max 1, but test shows it does not force 1 so a bit pointless
#SBATCH --array=11-14

# Exit the slurm script if a command fails
set -e

# Construct filenames dynamically
ID_STRING="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURMD_NODENAME}_GPU"
SLURM_OUTPUT_FILENAME="${ID_STRING}.out"
SLURM_ERROR_FILENAME="${ID_STRING}.err"
SLURM_OUTPUT_FILEPATH="/lisc/scratch/dome/pullen/GlobDB/${SLURM_OUTPUT_FILENAME}"
SLURM_ERROR_FILEPATH="/lisc/scratch/dome/pullen/GlobDB/${SLURM_ERROR_FILENAME}"

# Redirect output and error to my generated filenames
exec 1>${SLURM_OUTPUT_FILEPATH} 2>${SLURM_ERROR_FILEPATH}

# SPLIT_ID=$SLURM_ARRAY_TASK_ID
# instead of above make a zero padded SPLIT_ID here e.g. 1 -> 001
# to match the seqkit split format
# printf prints the formatted string not to standard output but instead assigns it to the shell variable specified by -v
printf -v SPLIT_ID "%03d" "$SLURM_ARRAY_TASK_ID"
FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/filtered1001AAmin_sorted.fasta.split/filtered1001AAmin_sorted.part_${SPLIT_ID}.fasta"

seqkit stats ${FASTA_INPUT_FILE}

MAX_SEQ_LEN=2000
MAX_RESIDUES=16000
MAX_BATCH=200

# JOB_PARAM_STRING="${MAX_RESIDUES}_${MAX_SEQ_LEN}_${MAX_BATCH}_${SPLIT_ID}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Array ID: ${SLURM_ARRAY_JOB_ID}"
echo "TMPDIR: ${TMPDIR}"
echo "Node: ${SLURMD_NODENAME}"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo "Running task $IDX with parameters:"
echo "  MAX_RESIDUES: $MAX_RESIDUES"
echo "  MAX_SEQ_LEN:  $MAX_SEQ_LEN"
echo "  MAX_BATCH:    $MAX_BATCH"

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
  --output $TMPDIR/embeddings_1001AAmin_${SPLIT_ID}.h5 \
  --log $TMPDIR/${SLURM_ARRAY_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log \
  --max_residues ${MAX_RESIDUES} --max_seq_len ${MAX_SEQ_LEN} --max_batch ${MAX_BATCH} \
  --master_embedding_file /lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/embeddings/chloroflexi_test100_embeddings.h5

if [ -f "$TMPDIR/embeddings_1001AAmin_${SPLIT_ID}.h5" ]; then
    cp "$TMPDIR/embeddings_1001AAmin_${SPLIT_ID}.h5" /lisc/scratch/dome/pullen/GlobDB/
else
    echo "File $TMPDIR/embeddings_1001AAmin_${SPLIT_ID}.h5 not found; skipping copy."
fi

if [ -f "$TMPDIR/${SLURM_ARRAY_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log" ]; then
    cp "$TMPDIR/${SLURM_ARRAY_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log" /lisc/scratch/dome/pullen/GlobDB/
else
    echo "File $TMPDIR/${SLURM_ARRAY_JOB_ID}_chloroflexi_test1000_gpu_splits_1001AAmin_${SPLIT_ID}.log not found; skipping copy."
fi

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR
