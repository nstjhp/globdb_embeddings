# first go into the stats dir made from running part1chunking.sh
# and check the stats txt files have the columns from seqkit stats:
# particularly 1 is file, 4 is num_seqs, 7 is avg_len
# then do something like this:
# for file in clusters_more_than1.part_*.txt; do     tail -n 11 "$file" | awk '{print $1, $4, $7}' | sed 's/,//g'; done | awk 'BEGIN {print "file,num_seqs,avg_len"} {print $1 "," $2 "," $3}' > numseqs_avglen_per_part.csv

import pandas as pd
import numpy as np
import re

# Read the CSV file
df = pd.read_csv("numseqs_avglen_per_part.csv")

# Define a function to extract the "part" and "chunk" using regex
def extract_part_and_chunk(file_name):
    # Extract the 3 digits after "part_"
    part_match = re.search(r"part_(\d{3})", file_name)
    part = part_match.group(1) if part_match else None

    # Extract the "chunk" (either number_number or 1001AAmin)
    chunk_match = re.search(r"_(\d+_\d+\.fasta|filtered(1001AAmin)\.fasta)$", file_name)
    if chunk_match:
        chunk = chunk_match.group(2) if chunk_match.group(2) else chunk_match.group(1).replace(".fasta", "")
    else:
        chunk = None

    return part, chunk

# Apply the function to the "file" column and create new columns
df[["part", "chunk"]] = df["file"].apply(lambda x: pd.Series(extract_part_and_chunk(x)))

# Calculate time_per_avg_protein from our fit formula
df["time_per_avg_protein"] = (
    4.4562e-07 * df["avg_len"]**2 + 7.0125e-04 * df["avg_len"] + 1.3098e-03
)

# Reorder columns
df = df[["file", "part", "chunk", "num_seqs", "avg_len", "time_per_avg_protein"]]

# Make 'chunk' a categorical variable with the specified order
chunk_order = [
    "0_99", "100_199", "200_299", "300_399", "400_499",
    "500_599", "600_699", "700_799", "800_899", "900_1000", "1001AAmin"
]
df["chunk"] = pd.Categorical(df["chunk"], categories=chunk_order, ordered=True)

# floor should now be OK because we already have 1% protein embeddings which 
# will be skipped when we run the code
df["num_splits"] = np.floor(df['num_seqs']/(3600/df['time_per_avg_protein']))
# Add 2 to 'num_splits' for rows where 'chunk' is "0_99" as the formula is too optimistic here
df.loc[df["chunk"] == "0_99", "num_splits"] += 2
# Display the updated DataFrame
# print(df)

df[["file", "num_splits"]].to_csv("files_by_how_many_to_split.csv")

# with this file it's useful to:
# cp files_by_how_many_to_split.csv files_by_how_many_to_split_excluding_1001AAmin.csv
# vi files_by_how_many_to_split_excluding_1001AAmin.csv
# :g/min/d to remove the lines with the 1001 min as we don't split them - they will go
# directly to the L40s which has the RAM to deal with big proteins
# :%s/\.0$// to remove the decimal part which won't work with seqkit split
# e.g. invalid argument "82.0" for "-p, --by-part" flag: strconv.ParseInt: parsing "82.0": invalid syntax
# Our file without the 1001 min will now go to part2splitting.sh
