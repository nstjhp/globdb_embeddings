#!/usr/bin/env python3
import csv
import sys
import itertools

def process_clusters_sorted(input_file, output_file):
    """
    Process a sorted TSV file where the first column is the cluster representative.
    For each group, if there is more than one member, write a row to the output file with:
      - The representative (key)
      - The count (number of occurrences)
      - A comma-separated list of all members (from the second column)
      
    This method only holds one group in memory at a time.
    """
    with open(input_file, "r") as f_in, open(output_file, "w", newline='') as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        
        # Group rows by the representative (first column)
        for rep, group in itertools.groupby(reader, key=lambda row: row[0]):
            # Create a list of members for this cluster
            members = [row[1] for row in group]
            count = len(members)
            # Only output clusters with more than one member
            if count > 1:
                writer.writerow([rep, count, ",".join(members)])

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python get_clusters_more_than_1.py slurm-4625318/globdb_clusters.tsv slurm-4625318/clusters_more_than1_lists.tsv")
        sys.exit(1)
    process_clusters_sorted(sys.argv[1], sys.argv[2])

