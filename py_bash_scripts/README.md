# globdb_embeddings

A collection of Python, Bash, and Slurm scripts for generating and managing protein embeddings of GlobDB (https://globdb.org/home).

## Table of Contents

- [Overview](#overview)  
- [Methods](#methods)  
- [Results](#results)  

## Overview

`globdb_embeddings` is a set of scripts to:
- use linclust to cluster sequences in the entire GlobDB fasta file
- process large fasta files into smaller, more manageable files
- generate protein embeddings from the fasta files
- merge the embeddings into h5 files
- give file information and timings

## Methods

We use the ProtT5-XL-U50 protein large language model (see https://github.com/agemagician/ProtTrans).
The script they kindly provide was refactored to make it more friendly for running on our life sciences cluster [LiSC](https://lisc.univie.ac.at/).
For each protein its embedding is a vector of 1024 floating point numbers.
During testing we found that to embed the entire GlobDB would take over a year to run with our settings on the available GPUs.
You can see the results of this extensive testing in this [graph](https://github.com/nstjhp/globdb_embeddings/blob/master/docs/length_vs_time_interactive_plot.html).
We therefore clustered GlobDB at 40% identity using `linclust` to provide representative clusters and aimed to embed a representative for each cluster of size 2 or more.
A subset fasta of only the cluster representatives was made using seqkit.
To efficiently use HPC resources, which on LiSC have an optimal runtime of 1 hour per task, we first broke the subsetted fasta into 10 parts, then chunked the representative proteins into bins by size (length of amino acid sequence), and then split these bins into '1 T4 GPU hour'-sized fasta files using the quadratic-fit formula shown in the graph mentioned above.
Thus, we had 5649 tasks of 1 T4 GPU hour run in a SLURM array.
Proteins ≥ 1001 amino acids were run on a L40s GPU which has higher vRAM and is faster than a T4 GPU and were not split further.
Specifically that meant 10 tasks (1 per part) each of which ran in around 12 hours.

## Results

The clusters represent 660303898 proteins (approx 80% of GlobDB), and we were able to generate embeddings for almost all cluster representatives.
The ones we could not embed were generally proteins longer than 9000 amino acids as our GPUs do not have enough vRAM to successfully calculate the embeddings, but they were only a few hundred proteins overall.
In total we were able to calculate 82972511 embeddings from GlobDB, which when merged together are 321 Gb on disk.
