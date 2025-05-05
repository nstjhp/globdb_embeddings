#!/bin/bash
#SBATCH --job-name=prott5_embed_big002 #redo_fails
#SBATCH --output=/lisc/scratch/dome/pullen/GlobDB/outfiles/job%A_%N_GPU.out
#SBATCH --error=/lisc/scratch/dome/pullen/GlobDB/outfiles/job%A_%N_GPU.err
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU
#SBATCH --mem=8000M
#SBATCH --time=0-03:00:00
#SBATCH --partition=gpu # needs to use a GPU with more RAM for big proteins
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --exclude=node-c[01-02] # trick to exclude c nodes that for unknown reasons use 2 threads
#SBATCH --gres=gpu:l40s:1

# Exit the slurm script if a command fails
set -e

# ***********************
# SPLIT_ID="REDO_FAILS"
SPLIT_ID="part_002"
CHUNK="min1001"
# ***********************
JOB_PARAM_STRING="${SLURM_JOB_ID}_${CHUNK}_${SPLIT_ID}"

# Test these job commands
echo "Job ID: ${SLURM_JOB_ID}"
echo "TMPDIR: ${TMPDIR}"
echo "Node: ${SLURMD_NODENAME}"

# Try the unsplit >=1001 AAs fasta
FASTA_INPUT_FILE="/lisc/scratch/dome/pullen/GlobDB/fastas/bins/globdb_r226_all_prot_${SPLIT_ID}_filtered1001AAmin.fasta"

seqkit stats -b ${FASTA_INPUT_FILE}

MAX_RESIDUES=4000
MAX_SEQ_LEN=1000
MAX_BATCH=200
echo "  MAX_RESIDUES: $MAX_RESIDUES"
echo "  MAX_SEQ_LEN:  $MAX_SEQ_LEN"
echo "  MAX_BATCH:    $MAX_BATCH"
echo "JOB_PARAM_STRING: $JOB_PARAM_STRING"

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
python /lisc/project/dome/protein_embeddings/py_bash_scripts/prott5_embedder_globdb.py \
  --input ${FASTA_INPUT_FILE} \
  --output $TMPDIR/embed_${JOB_PARAM_STRING}.h5 \
  --log $TMPDIR/${JOB_PARAM_STRING}.log \
  --max_residues ${MAX_RESIDUES} --max_seq_len ${MAX_SEQ_LEN} --max_batch ${MAX_BATCH} \
  --master_embedding_file /lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/embeddings/master_embeddings.h5

if [ -f "$TMPDIR/embed_${JOB_PARAM_STRING}.h5" ]; then
    cp "$TMPDIR/embed_${JOB_PARAM_STRING}.h5" /lisc/scratch/dome/pullen/GlobDB/embeddings
else
    echo "File $TMPDIR/embed_${JOB_PARAM_STRING}.h5 not found; skipping copy."
fi

if [ -f "$TMPDIR/${JOB_PARAM_STRING}.log" ]; then
    cp "$TMPDIR/${JOB_PARAM_STRING}.log" /lisc/scratch/dome/pullen/GlobDB/logs
else
    echo "File $TMPDIR/${JOB_PARAM_STRING}.log not found; skipping copy."
fi

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

