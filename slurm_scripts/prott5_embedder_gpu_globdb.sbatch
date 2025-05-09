#!/bin/bash
#SBATCH --job-name=prott5_embed_array
#SBATCH --output=/lisc/scratch/dome/pullen/GlobDB/outfiles/job%A_%a_%N_GPU.out
#SBATCH --error=/lisc/scratch/dome/pullen/GlobDB/outfiles/job%A_%a_%N_GPU.err
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # for GPU usage, only need 1 CPU and more shards
#SBATCH --mem=8000M
#SBATCH --time=0-01:00:00
#SBATCH --partition=basic #,gpu 
# #SBATCH --gres=shard:1 # If you don't know the optimal CPU/GPU ratio of your job, we recommend to request as many shard resources as CPU cores.
#SBATCH --gres=gpu:t4:1
#SBATCH --exclude=node-c[01-02] # trick to exclude c nodes that for unknown reasons use 2 threads
# #SBATCH --threads-per-core=1 # docs not clear if this is min 1, or max 1, but test shows it does not force 1 so a bit pointless
#SBATCH --array=1-7

# Exit the slurm script if a command fails
set -e

# ***********************
CHUNK="900_1000" # change here and --array above
# CHUNK="0_99" # change here and --array above
PART="002"
# ***********************

# SPLIT_ID=$SLURM_ARRAY_TASK_ID
# instead of above make a zero padded SPLIT_ID here e.g. 1 -> 001
# to match the seqkit split format
# printf prints the formatted string not to standard output but instead assigns it to the shell variable specified by -v
printf -v SPLIT_ID "%03d" "$SLURM_ARRAY_TASK_ID"
FASTA_INPUT_FILE="/lisc/scratch/dome/pullen/GlobDB/fastas/splits/globdb_r226_all_prot_part_${PART}_filtered1000AAmax_sorted_${CHUNK}.part_${SPLIT_ID}.fasta"
# globdb_r226_all_prot_part_*_filtered1000AAmax_sorted_*[_1234]99.fasta
# globdb_r226_all_prot_part_001_filtered1000AAmax_sorted_0_99.part_001.fasta

seqkit stats -b ${FASTA_INPUT_FILE}

# Construct filenames dynamically
# ID_STRING="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURMD_NODENAME}_GPU"
JOB_PARAM_STRING="${SLURM_ARRAY_JOB_ID}_${PART}_${CHUNK}_${SPLIT_ID}"

MAX_SEQ_LEN=2000
MAX_RESIDUES=16000
MAX_BATCH=200

echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Array ID: ${SLURM_ARRAY_JOB_ID}"
echo "TMPDIR: ${TMPDIR}"
echo "Node: ${SLURMD_NODENAME}"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo "Running task $IDX with parameters:"
echo "  MAX_RESIDUES: $MAX_RESIDUES"
echo "  MAX_SEQ_LEN:  $MAX_SEQ_LEN"
echo "  MAX_BATCH:    $MAX_BATCH"
echo "JOB_PARAM_STRING: $JOB_PARAM_STRING"

# For use in the Python logging, here we set MY_SLURM_PROCESS_ID 
# and base it on whether this is an array job (-n checks if a string is non-empty)
if [[ -n "${SLURM_ARRAY_JOB_ID}" && -n "${SLURM_ARRAY_TASK_ID}" ]]; then
    # This is an array job
    export MY_SLURM_PROCESS_ID="${JOB_PARAM_STRING}"
#    export MY_SLURM_PROCESS_ID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
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
