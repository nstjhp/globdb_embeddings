#!/bin/bash
#SBATCH --job-name=merge_h5_part001
#SBATCH --output=job%A_%x_%N.out
#SBATCH --error=job%A_%x_%N.err
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16000M
#SBATCH --time=00:30:00
#SBATCH --partition=basic,short

# Exit the slurm script if a command fails
set -e

ID="part001"
INPUT_FILES="/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_462*,/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_46330*"
OUTPUT_FILE="${ID}_embeddings.h5"

# Small to inter
/usr/bin/time python /lisc/project/dome/protein_embeddings/py_bash_scripts/merge_h5_small_to_inter.py --input-pattern "${INPUT_FILES}" --output-dir $TMPDIR --num-groups 16

echo "Small to intermediates DONE"

# Inter to 1 big master
/usr/bin/time python /lisc/project/dome/protein_embeddings/py_bash_scripts/merge_h5_inter_to_big.py --input-dir $TMPDIR --output-file $TMPDIR/$OUTPUT_FILE

if [ -f "$TMPDIR/$OUTPUT_FILE" ]; then
    cp $TMPDIR/$OUTPUT_FILE /lisc/scratch/dome/pullen/GlobDB/embeddings
else
    echo "File  $TMPDIR/$OUTPUT_FILE not found, copying intermediate results"
    cp -r $TMPDIR /lisc/scratch/dome/pullen/GlobDB/embeddings
fi

# Append the contents of this script to the output file
echo "=== Job Script Contents ==="
cat $0
echo "==========================="

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

