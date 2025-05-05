#!/bin/bash
#SBATCH --job-name=prott5_embed_single
# #SBATCH --output=job%A_GPU.out
# #SBATCH --error=job%A_GPU.err
#SBATCH --mail-type=END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU and more shards
#SBATCH --mem=10000M
#SBATCH --time=0-01:45:00
#SBATCH --partition=basic#,gpu 
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:t4:1

# Exit the slurm script if a command fails
set -e

# ***********************
CHUNK="100_199" # change here
SPLIT_ID="012"
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
echo "Node: ${SLURMD_NODENAME}"

FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/filtered1000AAmax_sorted_${CHUNK}.fasta.split/filtered1000AAmax_sorted_${CHUNK}.part_${SPLIT_ID}.fasta"

MAX_RESIDUES=16000
MAX_SEQ_LEN=2000
MAX_BATCH=200
JOB_PARAM_STRING="${MAX_RESIDUES}_${MAX_SEQ_LEN}_${MAX_BATCH}_${CHUNK}_${SPLIT_ID}"

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
  --output /lisc/scratch/dome/pullen/GlobDB/embeddings_${JOB_PARAM_STRING}.h5 \
  --log /lisc/scratch/dome/pullen/GlobDB/${JOB_PARAM_STRING}_chloroflexi_test1000_gpu_splits.log \
  --max_residues ${MAX_RESIDUES} --max_seq_len ${MAX_SEQ_LEN} --max_batch ${MAX_BATCH} \
  --master_embedding_file /lisc/scratch/dome/pullen/GlobDB/backup_embeddings/embeddings_100_199_012.h5
#  --master_embedding_file /lisc/scratch/dome/pullen/GlobDB/all_embeddings.h5

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

