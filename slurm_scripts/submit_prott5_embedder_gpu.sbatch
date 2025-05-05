#!/bin/bash
#SBATCH --job-name=prott5_embed
# #SBATCH --output=job%A_GPU.out
# #SBATCH --error=job%A_GPU.err
#SBATCH --mail-type=END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU and more shards
#SBATCH --mem=16000M
#SBATCH --time=0-00:30:00
#SBATCH --partition=basic#,gpu 
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:t4:1

# Exit the slurm script if a command fails
set -e

# Construct filenames dynamically
ID_STRING="${SLURM_JOB_ID}_${SLURMD_NODENAME}_GPU"
SLURM_OUTPUT_FILENAME="${ID_STRING}.out"
SLURM_ERROR_FILENAME="${ID_STRING}.err"
SLURM_OUTPUT_FILEPATH="/lisc/scratch/cube/pullen/Ecoli/${SLURM_OUTPUT_FILENAME}"
SLURM_ERROR_FILEPATH="/lisc/scratch/cube/pullen/Ecoli/${SLURM_ERROR_FILENAME}"

# Redirect output and error to my generated filenames
exec 1>${SLURM_OUTPUT_FILEPATH} 2>${SLURM_ERROR_FILEPATH}

# Test these job commands
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"

FASTA_INPUT_FILE="/lisc/project/cube/protein_embeddings/Ecoli/4998lines.fasta"
# FASTA_INPUT_FILE="/lisc/project/cube/protein_embeddings/Ecoli/uniprotkb_taxonomy_id_562_AND_existence_2025_02_21.fasta"
# FASTA_INPUT_FILE="/lisc/project/cube/protein_embeddings/Human/uniprotkb_taxonomy_id_9606_AND_existence_2025_02_25.fasta"
# FASTA_INPUT_FILE="/lisc/project/cube/protein_embeddings/GlobDB/chloroflexi_test100.faa"

# Run the Python script
python  -m cProfile -s cumulative /lisc/project/cube/protein_embeddings/prott5_embedder_globdb.py \
  --input ${FASTA_INPUT_FILE} \
  --output /lisc/scratch/cube/pullen/Ecoli/embeddings_confirmAllOK.h5 \
  --log /lisc/scratch/cube/pullen/testing_ecoli_gpu.log
#  --output /lisc/scratch/cube/pullen/Human/embeddings_GPU_maxresid8k_seqlen2k.h5 \

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "=========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

