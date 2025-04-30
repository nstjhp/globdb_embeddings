#!/bin/bash
#SBATCH --job-name=merge_h5_1001AAmin
#SBATCH --output=job%A_%x_%N.out
#SBATCH --error=job%A_%x_%N.err
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=18000M
#SBATCH --time=1:00:00

# Exit the slurm script if a command fails
set -e

INPUT_FILES="/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_4899843*,/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_4912031*"
OUTPUT_FILE="globdb_linclust_embeddings_1001AAmin.h5"

# Small to inter
/usr/bin/time python /lisc/project/dome/protein_embeddings/py_bash_scripts/merge_h5_small_to_inter.py --input-patterns "${INPUT_FILES}" --output-dir $TMPDIR --num-groups 5

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

