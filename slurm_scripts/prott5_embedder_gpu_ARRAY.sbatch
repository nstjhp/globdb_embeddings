#!/bin/bash
#SBATCH --job-name=prott5_embed_array
# #SBATCH --output=job%A_%a_GPU.out
# #SBATCH --error=job%A_%a_GPU.err
#SBATCH --mail-type=ARRAY_TASKS,END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU and more shards
#SBATCH --mem=16000M
#SBATCH --time=0-00:50:00
#SBATCH --partition=basic#,gpu 
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:t4:1
#SBATCH --array=0-5

# Exit the slurm script if a command fails
set -e

# Construct filenames dynamically
ID_STRING="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURMD_NODENAME}_GPU"
SLURM_OUTPUT_FILENAME="${ID_STRING}.out"
SLURM_ERROR_FILENAME="${ID_STRING}.err"
SLURM_OUTPUT_FILEPATH="/lisc/scratch/dome/pullen/Human/${SLURM_OUTPUT_FILENAME}"
SLURM_ERROR_FILEPATH="/lisc/scratch/dome/pullen/Human/${SLURM_ERROR_FILENAME}"

# Redirect output and error to my generated filenames
exec 1>${SLURM_OUTPUT_FILEPATH} 2>${SLURM_ERROR_FILEPATH}

# FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/Human/uniprotkb_taxonomy_id_9606_AND_existence_2025_02_25.fasta"
SPLIT_ID="005"
FASTA_INPUT_FILE="/lisc/project/dome/protein_embeddings/Human/human_proteome_filtered1000AAmax_sorted.fasta.split/human_proteome_filtered1000AAmax_sorted.part_${SPLIT_ID}.fasta"

MAX_SEQ_LEN=2000

MAX_RESIDUES_LIST=(4000 4000 8000 8000 16000 16000)
MAX_BATCH_LIST=(100 200 100 200 100 200)

IDX=$SLURM_ARRAY_TASK_ID
MAX_RESIDUES=${MAX_RESIDUES_LIST[$IDX]}
MAX_BATCH=${MAX_BATCH_LIST[$IDX]}

JOB_PARAM_STRING="${MAX_RESIDUES}_${MAX_SEQ_LEN}_${MAX_BATCH}_${SPLIT_ID}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo "Running task $IDX with parameters:"
echo "  MAX_RESIDUES: $MAX_RESIDUES"
echo "  MAX_SEQ_LEN:  $MAX_SEQ_LEN"
echo "  MAX_BATCH:    $MAX_BATCH"

# Run the Python script
python /lisc/project/dome/protein_embeddings/prott5_embedder_nick.py \
  --input ${FASTA_INPUT_FILE} \
  --output /lisc/scratch/dome/pullen/Human/embeddings_GPU_${JOB_PARAM_STRING}.h5 \
  --log /lisc/scratch/dome/pullen/Human/human_gpu_splits_${JOB_PARAM_STRING}.log \
  --max_residues ${MAX_RESIDUES} --max_seq_len ${MAX_SEQ_LEN} --max_batch ${MAX_BATCH}

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

