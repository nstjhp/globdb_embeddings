from Bio import SeqIO
from collections import defaultdict
import sys
import hashlib

def sequence_hash(seq_str):
        return hashlib.md5(seq_str.encode('utf-8')).hexdigest()

def find_exact_duplicates_total(fasta_files):
    # Dictionary to map each sequence to the number of times it appears.
    seq_counts = defaultdict(int)
    total_records = 0
    
    # Process each FASTA file provided
    for fasta_path in fasta_files:
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_records += 1
            key = sequence_hash(str(record.seq))
            seq_counts[key] += 1
            if total_records % 10000000 == 0:
                print(f"Total records so far: {total_records}, key is {key}, seq_counts[key] is {seq_counts[key]}")
    
    # Count how many unique sequences have duplicates and sum the extra copies.
    num_proteins_with_duplicates = sum(1 for count in seq_counts.values() if count > 1)
    total_duplicate_occurrences = sum(count - 1 for count in seq_counts.values() if count > 1)
    
    return total_records, num_proteins_with_duplicates, total_duplicate_occurrences

def find_exact_duplicates_seqs(fasta_path):
    # Dictionary to map each sequence (as a string) to a list of record IDs
    seq_to_ids = defaultdict(list)
    total_records = 0
    
    # Iterate through each record in the FASTA file
    for record in SeqIO.parse(fasta_path, "fasta"):
        total_records += 1
        # Convert the sequence to a string to use as a dictionary key
        seq_str = str(record.seq)
        seq_to_ids[seq_str].append(record.id)
    
    # Filter out unique sequences; only keep those that appear more than once
    duplicate_sequences = {seq: ids for seq, ids in seq_to_ids.items() if len(ids) > 1}
    
    # Calculate the total number of duplicate occurrences (beyond the first instance)
    total_duplicate_occurrences = sum(len(ids) - 1 for ids in duplicate_sequences.values())
    
    return total_records, duplicate_sequences, total_duplicate_occurrences

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py input.fasta")
        sys.exit(1)
    
    fasta_files = sys.argv[1:]
    total_records, dup_proteins, dup_occurrences = find_exact_duplicates_total(fasta_files)

    print(f"Total number of FASTA records processed: {total_records}")
    print(f"Number of unique proteins with duplicates: {dup_proteins}")
    print(f"Total non-unique copies (extra occurrences): {dup_occurrences}")

#    fasta_file = sys.argv[1]
#    duplicates, total_dups = find_exact_duplicates(fasta_file)
#    
#    if duplicates:
#        print("Found duplicate sequences:")
#        for seq, ids in duplicates.items():
#            print(f"Sequence present in {len(ids)} records: {', '.join(ids)}")
#    else:
#        print("No duplicate sequences found.")
#    
#    print(f"Total duplicate occurrences (non-unique copies): {total_dups}")
#
