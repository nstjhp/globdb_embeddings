import h5py
import glob
import numpy as np
import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Merge intermediate HDF5 files into one master file with unlimited rows."
    )
    parser.add_argument(
        '--input-dir', type=str, required=True,
        help="Directory containing the intermediate HDF5 files (e.g., 'path/to/intermediate_dir')"
    )
    parser.add_argument(
        '--output-file', type=str, required=True,
        help="Path for the output master HDF5 file (e.g., 'master.h5')"
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

    print(f"Merging {len(intermediate_files)} files into {args.output_file}")
    print(f"Creating master dataset with unlimited rows and {n_cols} columns")

    # Create the master file and initialize datasets with 0 rows, but with an unlimited maxshape.
    with h5py.File(args.output_file, 'w') as master:
        master_emb = master.create_dataset(
            'embeddings', shape=(0, n_cols), maxshape=(None, n_cols),
            dtype='float32', chunks=(10000, n_cols)# chunks are memory storage size
        )
        # Create dataset for keys using a variable-length UTF-8 string dtype.
        dt = h5py.string_dtype(encoding='utf-8')
        master_keys = master.create_dataset(
            'keys', shape=(0,), maxshape=(None,), dtype=dt, chunks=(1_000_000,)
        )

        current_index = 0
        for file in intermediate_files:
            with h5py.File(file, 'r') as f:
                emb_block = f['embeddings'][:]
                keys_block = f['keys'][:]
                n_block = emb_block.shape[0]

                # Extend the master datasets to accommodate the new block.
                new_total = current_index + n_block
                master_emb.resize((new_total, n_cols))
                master_keys.resize((new_total,))
                
                master_emb[current_index:new_total, :] = emb_block
                master_keys[current_index:new_total] = keys_block
                current_index = new_total

                print(f"Copied {n_block} embeddings from {file}")
                sys.stdout.flush()

        print(f"Finished merging. Total embeddings copied: {current_index}")

if __name__ == '__main__':
    main()

# python /lisc/project/dome/protein_embeddings/py_bash_scripts/merge_h5_inter_to_big.py --input-dir $TMPDIR --output-file $TMPDIR/master_embeddings.h5
