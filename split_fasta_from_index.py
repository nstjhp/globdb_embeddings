from itertools import islice
import sys
from typing import List

def split_fasta(fasta_file: str, fai_file: str, part_size: int, output_string: str) -> None:
    """
    Splits a large FASTA file into smaller parts using its .fai index file.

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

    Key Steps:
    0. Make the .fai file e.g. `seqkit faidx file.fasta`.
    1. Use `islice` to read every `part_size`th line from the .fai file which avoids the need to
       read the entire file into memory. See https://docs.python.org/3/library/itertools.html#itertools.islice
       and https://stackoverflow.com/a/27108718
    2. For each part:
       - Extract the sequence ID and calculate its length in bytes.
       - Add 2 bytes (for '>' and newline) to account for the FASTA header format.
       - Adjust the byte offset to include the ID line.
    3. Write each part to a separate output file, ensuring valid FASTA format.

    Args:
        fasta_file (str): Path to the input FASTA file.
        fai_file (str): Path to the corresponding .fai index file.
        part_size (int): Number of sequences per output part.
        output_string (str): Template for output filenames (e.g., "path/to/file/prefix_"). `{part_num:03d}.fasta` will be appended.

    Output:
        Generates multiple FASTA files, each containing `part_size` sequences, with filenames
        formatted using `output_string`. For example, "part_001.fasta", "part_002.fasta", etc.

    Notes:
        - The first part always starts at byte 0 to include the header of the first sequence.
        - The script assumes the FASTA file is well-formatted, with each sequence starting
          with a '>' followed by the sequence ID and a newline.
        - The zero-padding is 3 digits in the output filename.
        - Having non-ASCII characters in the ID string might cause problems.
    """
    # Input validation
    assert part_size > 0, "part_size must be a positive integer."
    assert isinstance(fasta_file, str) and fasta_file.endswith(('.fasta', '.faa')), "fasta_file must be a valid .fasta or .faa file path."
    assert isinstance(fai_file, str) and fai_file.endswith('.fai'), "fai_file must be a valid .fai file path."
    assert isinstance(output_string, str) and "{part_num" in output_string, "output_string must be a valid format string with '{part_num}'."

    # Read the .fai file and calculate offsets
    offsets: List[int] = []
    try:
        with open(fai_file, 'r') as f:
            for line in islice(f, 0, None, part_size):
                cols = line.strip().split('\t')
                if len(cols) != 5:
                    raise ValueError(f"Invalid .fai file format in line: {line}")

                seq_id = cols[0]  # Sequence ID
                seq_offset = int(cols[2])  # Byte offset of the sequence

                # Calculate the byte position of the ID line
                id_length = len(seq_id.encode('utf-8'))  # Length of the ID in bytes
                id_line_bytes = id_length + 2  # Add 2 bytes for '>' and newline
                id_offset = seq_offset - id_line_bytes  # Start of the ID line

                offsets.append(id_offset)
    except FileNotFoundError:
        raise FileNotFoundError(f"The .fai file '{fai_file}' does not exist.")
    except Exception as e:
        raise RuntimeError(f"Error reading .fai file: {e}")

    # Ensure the first part starts at byte 0 i.e. start of first sequences's ID
    offsets[0] = 0

#    with open("/lisc/scratch/dome/pullen/GlobDB/fastas/offsets_for_splitting.txt", 'w') as offsets_file:
#        for offset in offsets:
#            offsets_file.write(f"{offset}\n")
#    print(f"len offset: {len(offsets)}")

    # Open the FASTA file and extract parts
    try:
        # Open the FASTA file for reading
        with open(fasta_file, 'rb') as fasta:
            part_num = 1
            total_parts = len(offsets)
    
            for i in range(total_parts):
                # Calculate the start and end of the part
                start_offset = offsets[i]
                if i + 1 < total_parts:
                    end_offset = offsets[i + 1]
                else:
                    # For the last part, read until the end of the file
                    end_offset = None
    
                    # Generate the output filename
                    output_file = output_string.format(part_num=part_num)
    
                    # Write the part to the output file
                    try:
                        with open(output_file, 'wb') as out_f:
                            # Seek the start of the part
                            fasta.seek(start_offset)
                            # With file.read(size) at most size characters (in text mode) or size bytes (in binary mode) are read and returned.
                            # So we read all that we need in 1 go
                            if end_offset is not None:
                                part_bytes = fasta.read(end_offset - start_offset)
                            else:
                                # Now it will read until the end of the file
                                part_bytes = fasta.read()
                            # Write the whole part to the output file
                            out_f.write(part_bytes)
                        print(f"Written {output_file}")
                        sys.stdout.flush()
                    except IOError as e:
                        raise IOError(f"Error writing to output file '{output_file}': {e}")
                    part_num += 1
    except FileNotFoundError:
        raise FileNotFoundError(f"The FASTA file '{fasta_file}' does not exist.")
    except Exception as e:
        raise RuntimeError(f"Error processing FASTA file: {e}")

def main():
    # Parameters
    fasta_file = "/lisc/scratch/dome/pullen/GlobDB/fastas/globdb_r226_all_prot.faa"
    fai_file = "/lisc/scratch/dome/pullen/GlobDB/fastas/globdb_r226_all_prot.faa.seqkit.fai"
    output_string = "/lisc/scratch/dome/pullen/GlobDB/fastas/globdb_r226_all_prot_part_{part_num:03d}.fasta"  # Zero-padded filenames
    part_size = 3354462  #sequences per split file, based on total sequences/some round number: 838615274/250=3354461.096
    
    # Run the splitting
    split_fasta(fasta_file, fai_file, part_size, output_string)

if __name__ == '__main__':
    main()
