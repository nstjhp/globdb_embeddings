#!/usr/bin/env python3
import os
import glob
import h5py
import argparse

def print_h5_info(file_path):
    # Get file size on disk in bytes.
    file_size = os.path.getsize(file_path)
    print(f"File: {file_path}")
    print(f"  File size on disk: {file_size} bytes")
    
    with h5py.File(file_path, 'r') as f:
        # Embeddings dataset info
        if 'embeddings' in f:
            ds = f['embeddings']
            print("  Dataset 'embeddings':")
            print(f"    Shape: {ds.shape}")
            print(f"    Data type: {ds.dtype}")
            # Estimate the in-memory size (uncompressed)
            est_size = ds.size * ds.dtype.itemsize
            print(f"    Estimated in-memory size: {est_size} bytes")
        else:
            print("  Dataset 'embeddings' not found!")
        
        # Keys dataset info
        if 'keys' in f:
            ds = f['keys']
            print("  Dataset 'keys':")
            print(f"    Shape: {ds.shape}")
            print(f"    Data type: {ds.dtype}")
        else:
            print("  Dataset 'keys' not found!")
    print("-" * 50)

def main():
    parser = argparse.ArgumentParser(
        description="Compare HDF5 files by printing size, shape, and datatype info."
    )
    parser.add_argument(
        "--input-pattern", type=str, required=True,
        help="Glob pattern to match HDF5 files (e.g. '/path/to/*.h5')"
    )
    args = parser.parse_args()

    # Expand the glob pattern into a list of file paths.
    files = sorted(glob.glob(args.input_pattern))
    if not files:
        print(f"No files matched the pattern: {args.input_pattern}")
        return

    for file_path in files:
        print_h5_info(file_path)

if __name__ == "__main__":
    main()

# /usr/bin/time python print_h5_info.py --input-pattern "../GlobDB/embeddings/*.h5"
