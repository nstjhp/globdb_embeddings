#!/bin/bash
# Submit the array job and capture its job ID.
array_job_id=$(sbatch --parsable prott5_embedder_gpu_globdb.sbatch)
echo "Submitted array job with ID: $array_job_id"

# Submit the merge job with a dependency on the array job.
merge_job_id=$(sbatch --dependency=afterok:$array_job_id --parsable merge_checkpoints.sh)
echo "Submitted merge job with ID: $merge_job_id"

