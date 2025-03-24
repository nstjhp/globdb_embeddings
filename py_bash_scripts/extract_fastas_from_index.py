from itertools import islice
import sys
from typing import List, Set, Tuple
import os

def extract_seqs(fasta_file: str, fai_file: str, IDs_to_extract: str, output_file: str) -> None:
    """
    Extract desired sequences from a large FASTA file using its .fai index file.

    The .fai (FASTA index) file is a tab-separated, 5-column file with the following fields:
    1. NAME: Name of the reference sequence.
    2. LENGTH: Total length of the sequence, in bases.
    3. OFFSET: Byte offset in the FASTA file of the sequence's first base.
    4. LINEBASES: Number of bases per line in the sequence.
    5. LINEWIDTH: Number of bytes per line, including the newline character.

    The .fai file can be generated using tools like `seqkit faidx` or `samtools faidx`.

    Example of a .fai file:
    ```
    AMXMAG_0088___0 78      17      78      79
    AMXMAG_0088___1 689     113     689     690
    AMXMAG_0088___2 182     820     182     183
    ```

    We can exploit the knowledge of which sequence starts at which byte (column 3) and the standardised
    format of the fasta to calculate the byte at which the ID of the sequence starts too. This makes
    splitting the file into parts much more efficient than reading the whole fasta into memory.
    A useful command to help is `od -A d -t c -w16 file.fasta` to see the byte positions and characters.
    In the example above we have ID strings of 15 bytes, plus a `>` and a newline. So from the byte 
    offset we subtract 15+2 to get the start of the sequence ID (which should = the previous sequence's
    offset plus linewidth, as we have the sequence on 1 line).

    Key Steps:
    0. Make the .fai file e.g. `seqkit faidx file.fasta`.
    1. Use `islice` to read every matching ID's line from the .fai file which avoids the need to
       read the entire file into memory. See https://docs.python.org/3/library/itertools.html#itertools.islice
       and https://stackoverflow.com/a/27108718
    2. For each desired sequence:
       - Extract the sequence ID and calculate its length in bytes.
       - Add 2 bytes (for '>' and newline) to account for the FASTA header format.
       - Adjust the byte offset to include the ID line.
    3. Write each sequence to a separate output file, ensuring valid FASTA format.

    Args:
        fasta_file (str): Path to the input FASTA file.
        fai_file (str): Path to the corresponding .fai index file.
        IDs_to_extract (str): Path to a text file containing sequence IDs to extract (one per line).
        output_file (str): Output filename

    Output:
        Generates one FASTA file containing only our desired sequences.

    Notes:
        - The first part always starts at byte 0 to include the header of the first sequence.
        - The script assumes the FASTA file is well-formatted, with each sequence starting
          with a '>' followed by the sequence ID and a newline.
        - Having non-ASCII characters in the ID string might cause problems.
    """
    # Input validation
    assert isinstance(fasta_file, str) and fasta_file.endswith(('.fasta', '.faa')), "fasta_file must be a valid .fasta or .faa file path."
    assert isinstance(fai_file, str) and fai_file.endswith('.fai'), "fai_file must be a valid .fai file path."
    assert isinstance(IDs_to_extract, str), "IDs_to_extract must be a valid file path as a string."
    assert isinstance(output_file, str) and output_file.endswith(('.fasta', '.faa')), "output_file must be a valid .fasta or .faa file path."
    if not os.path.exists(IDs_to_extract):
        raise FileNotFoundError(f"IDs_to_extract file '{IDs_to_extract}' not found.")

    ids_set: Set[str] = set()
    with open(IDs_to_extract, 'r') as id_file:
        for line in id_file:
            # One ID per line, stripping any whitespace
            seq_id = line.strip()
            if seq_id:
                ids_set.add(seq_id)

    # Get the total size of the FASTA file (needed for the last sequence)
    fasta_size = os.stat(fasta_file).st_size

    # Helper: given a line from the .fai file, compute the header start offset and return also the seq ID.
    def parse_fai_line(line: str) -> Tuple[int, str]:
        cols = line.rstrip("\n").split('\t')
        if len(cols) != 5:
            raise ValueError(f"Invalid .fai file format in line: {line}")
        seq_id = cols[0] # Sequence ID
        seq_offset = int(cols[2]) # Byte offset of the sequence
        # Calculate the byte position of the ID line: 
        # subtract header length (ID encoded in bytes + 2 for '>' and newline)
        header_start = seq_offset - (len(seq_id.encode('utf-8')) + 2)
        return header_start, seq_id

    line_count = 0
    extracted_count = 0
    # Scan the .fai file and record intervals (start, end) for matching sequences.
    # If consecutive intervals are contiguous, merge them.
    matching_intervals = []  # List of tuples (start_offset, end_offset)
    with open(fai_file, 'r') as fai:
        prev_line = None
        prev_header_start = None
        prev_seq_id = None
        
        for line in fai:
            line_count += 1
            if line_count % 10000000 == 0:
                print(f"Processed {line_count} lines of the .fai file so far...")
                print(f"Extracted {extracted_count} sequences so far...")
                sys.stdout.flush()
            if prev_line is None:
                prev_header_start, prev_seq_id = parse_fai_line(line)
                prev_line = line
                continue
            
            # For current line, compute its header start offset.
            current_header_start, _ = parse_fai_line(line)
            
            # If the previous sequence is one of the IDs to extract,
            # record its interval as [prev_header_start, current_header_start)
            if prev_seq_id in ids_set:
                extracted_count += 1
                current_interval = (prev_header_start, current_header_start)
                # Merge with the previous interval if contiguous.
                if matching_intervals and matching_intervals[-1][1] == current_interval[0]:
                    matching_intervals[-1] = (matching_intervals[-1][0], current_interval[1])
                else:
                    matching_intervals.append(current_interval)
            
            # Advance to next line.
            prev_header_start, prev_seq_id = current_header_start, parse_fai_line(line)[1]
        
        # Process the last line in the .fai file.
        if prev_seq_id in ids_set:
            extracted_count += 1
            last_interval = (prev_header_start, fasta_size)
            if matching_intervals and matching_intervals[-1][1] == last_interval[0]:
                matching_intervals[-1] = (matching_intervals[-1][0], last_interval[1])
            else:
                matching_intervals.append(last_interval)
    
    # Now open the FASTA file once and write out each matching interval.
    with open(fasta_file, 'rb') as fasta, open(output_file, 'wb') as out_f:
        for start, end in matching_intervals:
            # Seek to the beginning of the interval.
            fasta.seek(start)
            # Read exactly the number of bytes for this interval.
            data = fasta.read(end - start)
            out_f.write(data)

    print(f"Extraction complete. {len(matching_intervals)} intervals were written to {output_file}.")
    print(f"Processed {line_count} lines of the .fai file in total.")
    print(f"Total sequences extracted: {extracted_count}")

#    with open("/lisc/scratch/dome/pullen/GlobDB/fastas/offsets_for_splitting.txt", 'w') as offsets_file:
#        for offset in offsets:
#            offsets_file.write(f"{offset}\n")
#    print(f"len offset: {len(offsets)}")

def main():
    # Parameters
    fasta_file = "/lisc/scratch/dome/pullen/GlobDB/linclust/slurm-4625318/globdb_clusters_rep.fasta"
    fai_file = "/lisc/scratch/dome/pullen/GlobDB/linclust/slurm-4625318/globdb_clusters_rep.fasta.fai"
    output_file = "/lisc/scratch/dome/pullen/GlobDB/linclust/slurm-4625318/clusters_more_than1.fasta"
    IDs_to_extract = "/lisc/scratch/dome/pullen/GlobDB/linclust/slurm-4625318/cluster_more_than1_IDs.txt"
    
    # Run the splitting
    extract_seqs(fasta_file, fai_file, IDs_to_extract, output_file)

if __name__ == '__main__':
    main()
