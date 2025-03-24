#!/bin/bash
#SBATCH --job-name=globdb_linclust
#SBATCH --output=job%A_globdb_linclust.out
#SBATCH --error=job%A_globdb_linclust.err
#SBATCH --mail-type=END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=1400000M
#SBATCH --time=2-00:00:00
#SBATCH --partition=himem

# Exit the slurm script if a command fails
set -e

module load mmseqs2

/usr/bin/time cp /lisc/scratch/dome/speth/resources_scratch/globdb/globdb_r226_all_prot.faa $TMPDIR
echo "TIME globdb_r226_all_prot.faa copied to $TMPDIR"

mkdir $TMPDIR/tmp

/usr/bin/time mmseqs createdb $TMPDIR/globdb_r226_all_prot.faa $TMPDIR/globdb
echo "TIME Created DB"
/usr/bin/time mmseqs linclust $TMPDIR/globdb $TMPDIR/globdb_clusters $TMPDIR/tmp --min-seq-id 0.4 --threads 100 --split-memory-limit 1300G --remove-tmp-files 1
echo "TIME Linclustered"
/usr/bin/time mmseqs createtsv $TMPDIR/globdb $TMPDIR/globdb $TMPDIR/globdb_clusters $TMPDIR/globdb_clusters.tsv
echo "TIME TSV created"
/usr/bin/time mmseqs createsubdb $TMPDIR/globdb_clusters $TMPDIR/globdb $TMPDIR/globdb_clusters_rep
echo "TIME Representative sequences subDB made"
/usr/bin/time mmseqs convert2fasta $TMPDIR/globdb_clusters_rep $TMPDIR/globdb_clusters_rep.fasta
echo "TIME Representative sequences extracted to fasta"

rm $TMPDIR/globdb_r226_all_prot.faa

/usr/bin/time cp -r $TMPDIR /lisc/scratch/dome/pullen/GlobDB/linclust
echo "TIME All files copied from $TMPDIR to /lisc/scratch/dome/pullen/GlobDB/linclust"

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

