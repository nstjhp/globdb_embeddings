import h5py
import glob
import numpy as np
import os
import argparse
from multiprocessing import Pool
import sys

def process_group(args):
    """
    Process a group of HDF5 source files.
    For each file, iterate over its keys (each embedding),
    read the 1024-d float array, and store it along with its key.
    Finally, write an intermediate master file.
    """
    group_files, output_filename = args

    # First pass: count total embeddings in this group
    total_embeddings = 0
    for fname in group_files:
        with h5py.File(fname, 'r') as f:
            total_embeddings += len(f.keys())

    print(f"Processing {len(group_files)} files into {output_filename} with {total_embeddings} embeddings")
    sys.stdout.flush()

    # Preallocate numpy arrays:
    embeddings = np.empty((total_embeddings, 1024), dtype=np.float32)
    keys = []  # Will be converted to an array of strings later

    # Second pass: read each embedding and key
    idx = 0
    for fname in group_files:
        with h5py.File(fname, 'r') as f:
            for key in f.keys():
                data = f[key][:]
                if data.shape != (1024,):
                    raise ValueError(f"Unexpected shape for key {key} in file {fname}: {data.shape}")
                embeddings[idx] = data
                keys.append(key)
                idx += 1

    # Write the aggregated data to an intermediate file.
    with h5py.File(output_filename, 'w') as f_out:
        f_out.create_dataset('embeddings', data=embeddings)
        # Create a dataset for keys using a variable-length UTF-8 string type.
        dt = h5py.string_dtype(encoding='utf-8')
        keys_array = np.array(keys, dtype=object)
        f_out.create_dataset('keys', data=keys_array, dtype=dt)

    print(f"Finished writing {output_filename} with {idx} embeddings")
    sys.stdout.flush()
    return output_filename

def main():
    parser = argparse.ArgumentParser(description="Process and aggregate HDF5 files into intermediate master files.")
    parser.add_argument('--input-patterns', type=str, required=True,
                        help="Comma-separated glob patterns for source HDF5 files (e.g., 'embed_462*,embed_46330*')")
    parser.add_argument('--output-dir', type=str, required=True,
                        help="Directory where intermediate master files will be stored")
    parser.add_argument('--num-groups', type=int, default=16,
                        help="Number of groups/intermediate files to create (default: 16)")
    args = parser.parse_args()

    # Split the comma-separated patterns and combine the results.
    patterns = [p.strip() for p in args.input_patterns.split(',')]
    source_files = []
    for pattern in patterns:
        source_files.extend(glob.glob(pattern))

    if not source_files:
        print("No source files found. Exiting.")
        return
    else:
        print(f"Found {len(source_files)} files")
        print("Files:", source_files)

    # Divide the list of files into groups.
    groups = [source_files[i::args.num_groups] for i in range(args.num_groups)]

    # Prepare tasks: each task is a tuple (group_files, output_filename)
    tasks = []
    for i, group in enumerate(groups):
        output_filename = os.path.join(args.output_dir, f'intermediate_{i}.h5')
        tasks.append((group, output_filename))

    # Process groups in parallel using multiprocessing.
    with Pool(processes=args.num_groups) as pool:
        results = pool.map(process_group, tasks)

    print("Created intermediate files:", results)

if __name__ == '__main__':
    main()

# master_filename = "master_embeddings.h5"

# input_pattern = '/lisc/project/dome/protein_embeddings/GlobDB/chloroflexi_test1000/embeddings/*.h5'

# python /lisc/project/dome/protein_embeddings/py_bash_scripts/merge_h5_small_to_inter.py --input-patterns "/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_462*,/lisc/scratch/dome/pullen/GlobDB/embeddings/embed_46330*" --output-dir $TMPDIR --num-groups 16
