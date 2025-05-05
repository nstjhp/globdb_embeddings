#!/usr/bin/env python3
import h5py
import argparse

def save_keys(h5_path: str, out_txt: str, block_size: int = 10_000_000):
    """
    Read the 'keys' dataset from h5_path in blocks of block_size rows,
    decode any byte-strings to UTF-8, and write one key per line to out_txt.
    """
    with h5py.File(h5_path, 'r') as f, open(out_txt, 'w', encoding='utf-8') as outf:
        keys_ds = f['keys']
        total = keys_ds.shape[0]
        for start in range(0, total, block_size):
            end = min(start + block_size, total)
            block = keys_ds[start:end]  # array of bytes or str
            for k in block:
                # decode bytes, leave str alone
                if isinstance(k, (bytes, bytearray)):
                    outf.write(k.decode('utf-8') + '\n')
                else:
                    outf.write(k + '\n')
            print(f"  Wrote keys {start}–{end} / {total}", end='\r')
    print(f"\nDone—wrote {total} keys to {out_txt}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Dump HDF5 'keys' → text file.")
    p.add_argument("--h5",     required=True, help="Master HDF5 file path")
    p.add_argument("--out",    required=True, help="Output .txt path (one key per line)")
    p.add_argument("--block-size", type=int, default=10_000_000,
                   help="Rows to read at a time (default: 1e7)")
    args = p.parse_args()
    save_keys(args.h5, args.out, args.block_size)

