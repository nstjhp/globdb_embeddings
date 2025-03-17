import h5py
import glob
import numpy as np
import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Merge intermediate HDF5 files into one master file."
    )
    parser.add_argument(
        '--input-dir', type=str, required=True,
        help="Directory containing the intermediate HDF5 files (e.g., 'path/to/intermediate_dir')"
    )
    parser.add_argument(
        '--output-file', type=str, required=True,
        help="Path for the output master HDF5 file (e.g., 'master.h5')"
    )
    parser.add_argument(
        '--total-embeddings', type=int, default=3052886,
        help="Total number of embeddings to expect (default: 3052886)"
    )
    args = parser.parse_args()

    # Find all intermediate files in the input directory.
    intermediate_files = sorted(glob.glob(os.path.join(args.input_dir, "intermediate_*.h5")))
    if not intermediate_files:
        print(f"No intermediate files found in {args.input_dir}")
        return

    # Open one intermediate file to get the embedding width.
    with h5py.File(intermediate_files[0], 'r') as f:
        emb_shape = f['embeddings'].shape  # e.g., (n, 1024)
    n_cols = emb_shape[1]
    total_embeddings = args.total_embeddings

    print(f"Merging {len(intermediate_files)} files into {args.output_file}")
    print(f"Preallocating a dataset of shape ({total_embeddings}, {n_cols})")

    # Create the master file and preallocate datasets.
    with h5py.File(args.output_file, 'w') as master:
        master_emb = master.create_dataset(
            'embeddings', shape=(total_embeddings, n_cols), dtype='float32'
        )
        # Create dataset for keys using a variable-length UTF-8 string dtype.
        dt = h5py.string_dtype(encoding='utf-8')
        master_keys = master.create_dataset(
            'keys', shape=(total_embeddings,), dtype=dt
        )

        current_index = 0
        for file in intermediate_files:
            with h5py.File(file, 'r') as f:
                emb_block = f['embeddings'][:]
                keys_block = f['keys'][:]
                n_block = emb_block.shape[0]

                if current_index + n_block > total_embeddings:
                    raise ValueError(
                        f"Block from {file} exceeds allocated size: current index {current_index} + block size {n_block} > {total_embeddings}"
                    )

                master_emb[current_index:current_index + n_block, :] = emb_block
                master_keys[current_index:current_index + n_block] = keys_block
                current_index += n_block
                print(f"Copied {n_block} embeddings from {file}")
                sys.stdout.flush()

        print(f"Finished merging. Total embeddings copied: {current_index}")
        if current_index != total_embeddings:
            print(f"Warning: Total embeddings copied ({current_index}) does not match expected ({total_embeddings}).")

if __name__ == '__main__':
    main()

