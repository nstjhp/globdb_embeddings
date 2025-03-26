import h5py
import numpy as np
import glob
import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Merge large HDF5 files into one master file."
    )
    parser.add_argument('--input-files', type=str, required=True,
                        help="Comma-separated glob patterns for source HDF5 files (e.g., 'embed_462*,embed_46330*')")
    parser.add_argument('--output-file', type=str, required=True,
                        help="Path to the final output master HDF5 file.")
    parser.add_argument('--block-size', type=int, default=100000,
                        help="Number of rows to process per block (default: 100000).")
    args = parser.parse_args()

    # Split the comma-separated patterns and combine the results.
    patterns = [p.strip() for p in args.input_files.split(',')]
    input_files = []
    for pattern in patterns:
        input_files.extend(glob.glob(pattern))

    if not input_files:
        print("No input files found. Exiting.")
        return
    else:
        print(f"Found {len(input_files)} files")
        print("Files:", input_files)

    # Open one input file to determine the embedding width (assumed same for all)
    with h5py.File(input_files[0], 'r') as f:
        emb_shape = f['embeddings'].shape  # e.g., (n, 1024)
    n_cols = emb_shape[1]
    print(f"Found {len(input_files)} input files. Embedding width taken to be {n_cols}.")

    # Create the master file with unlimited rows.
    with h5py.File(args.output_file, 'w') as master:
        # Create the embeddings dataset.
        emb_ds = master.create_dataset(
            'embeddings',
            shape=(0, n_cols),
            maxshape=(None, n_cols),
            dtype='float32',
            chunks=(args.block_size, n_cols)  # Using block size for appending.
        )
        # Create the keys dataset with a chunk size of 1_000_000 for efficient reads.
        dt = h5py.string_dtype(encoding='utf-8')
        keys_ds = master.create_dataset(
            'keys',
            shape=(0,),
            maxshape=(None,),
            dtype=dt,
            chunks=(1_000_000,)
        )

        current_index = 0

        # Process each input file.
        for infile in input_files:
            print(f"Processing file: {infile}")
            with h5py.File(infile, 'r') as f:
                total_rows = f['embeddings'].shape[0]
                # Process the file in blocks for memory efficiency.
                for start in range(0, total_rows, args.block_size):
                    end = min(start + args.block_size, total_rows)
                    # Read keys in this block.
                    keys_block = f['keys'][start:end]
                    # Build a valid mask: exclude keys that are b'' or empty string.
                    valid_mask = np.array([k != b'' and k != "" for k in keys_block])
                    if not np.any(valid_mask):
                        # Skip block if no valid entries.
                        continue
                    valid_indices = np.where(valid_mask)[0]
                    # Read the corresponding embeddings.
                    emb_block = f['embeddings'][start:end]
                    valid_emb = emb_block[valid_indices, :]
                    valid_keys = keys_block[valid_indices]

                    n_valid = valid_emb.shape[0]
                    # Extend the master datasets.
                    new_total = current_index + n_valid
                    emb_ds.resize((new_total, n_cols))
                    keys_ds.resize((new_total,))
                    # Append the valid entries.
                    emb_ds[current_index:new_total, :] = valid_emb
                    keys_ds[current_index:new_total] = valid_keys
                    current_index = new_total

                    print(f"  Processed rows {start} to {end}: added {n_valid} valid entries (total so far: {current_index})")
                    sys.stdout.flush()

        print(f"Finished merging. Total valid embeddings in master file: {current_index}")

if __name__ == '__main__':
    main()

# /usr/bin/time python merge_h5_big_to_final.py --input-files "/lisc/project/dome/protein_embeddings/GlobDB/embeddings/chlor_plus_part001.h5,/lisc/project/dome/protein_embeddings/GlobDB/embeddings/part002_embeddings.h5" --output-file /lisc/project/dome/protein_embeddings/GlobDB/embeddings/chlor_plus_part001and002.h5
