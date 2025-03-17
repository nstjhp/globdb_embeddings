#!/bin/bash
#SBATCH --job-name=split_big_fasta
#SBATCH --output=job%A.out
#SBATCH --error=job%A.err
#SBATCH --mail-type=END,TIME_LIMIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=16000M # only for CPU memory
#SBATCH --mem=6000M # only for CPU memory
#SBATCH --time=0-00:30:00
#SBATCH --partition=basic,short#,himem

# Exit the slurm script if a command fails
set -e

# These commands will go into stdout
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "TMPDIR: ${TMPDIR}"

python split_fasta_from_index.py

# Compared to the seqkit way my nice script work pretty fast and with <3 GB memory requirement
# seff 4550358
# Job ID: 4550358
# Cluster: lisc
# User/Group: pullen/login
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:03:02
# CPU Efficiency: 16.67% of 00:18:12 core-walltime
# Job Wall-clock time: 00:18:12
# Memory Utilized: 2.91 GB
# Memory Efficiency: 49.65% of 5.86 GB (5.86 GB/node)

# Seqkit way OOMed even with 256 GB memory
# -2: two-pass mode read files twice to lower memory usage
# -U: For the two-pass mode (-2/--two-pass), The flag -U/--update-faidx is recommended to ensure the .fai file matches the FASTA file.
# -s: split sequences into multi parts with N sequences (here ~1/250 of GlobDB's fasta)
# seqkit split /lisc/scratch/dome/pullen/GlobDB/fastas/globdb_r226_all_prot.faa -2 -U -s 3354461 -O $TMPDIR

# cp -r $TMPDIR /lisc/scratch/dome/pullen/GlobDB/fastas

# If we reached this point, we succeeded. We clean up resources.
rm -rf $TMPDIR

