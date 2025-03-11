import pandas as pd
import argparse
from collections import Counter
import subprocess

def generate_histogram_table(input_fasta, bin_size, threshold = 4000, save_table=False):
    """
    Generates a histogram table of protein lengths from a FASTA file.
    It replaces the following shell command:
        seqkit fx2tab -n -l -i uniprotkb_taxonomy_id_9606_AND_existence_2025_02_25.fasta | cut -f 2 | sort -n | uniq -c > human_length_dist.txt
    because the sort will take too long on big files.

    This function uses `seqkit` to extract sequence lengths and generates a histogram table with bin ranges, counts, and percentages of the total.
    
    Args:
        input_fasta (str): Path to the input FASTA file.
        bin_size (int): Size of each bin for the histogram.
        threshold (int): Threshold for counting proteins above a certain length.
        save_table (bool): Save the histogram table as a CSV?
    Returns:
        None: The histogram table is printed to the console and optionally saved to a CSV file.
    """
    # Step 1: Use seqkit to extract lengths
    print("Getting sequence lengths...")
    seqkit_command = f"seqkit fx2tab -n -l -i {input_fasta}"
    process = subprocess.Popen(seqkit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error running seqkit: {stderr.decode('utf-8')}")
        return    

    # Step 2: Parse seqkit output and count lengths
    print("Counting sequence lengths...")
    lengths = [int(line.split('\t')[1]) for line in stdout.decode('utf-8').splitlines()]
    length_counter = Counter(lengths)

    # Step 3: Create a DataFrame from the length counts
    df = pd.DataFrame(length_counter.items(), columns=['Length', 'Count'])    

    # Bin with integer division
    df['Bin'] = (df['Length'] // bin_size) * bin_size
    
    # Group by bins and sum the counts
    histogram_table = df.groupby('Bin')['Count'].sum().reset_index()
    
    # Better column names
    histogram_table.columns = ['Bin Start', 'Count']
    
    # Add column for the bin range
    histogram_table['Bin Range'] = histogram_table['Bin Start'].apply(
        lambda x: f'{x}-{x + bin_size - 1}'
    )

    # Calculate the percentage of total for each bin
    total_count = histogram_table['Count'].sum()
    histogram_table['% of Total'] = (histogram_table['Count'] / total_count * 100).round(2)
    
    histogram_table = histogram_table.sort_values('Bin Start')
    print(histogram_table)

    # Save the table to a file
    if save_table:
        output_file = f"{input_fasta}_histogram_bin{bin_size}.csv"
        histogram_table.to_csv(output_file, index=False)
        print(f"Histogram table saved to {output_file}")
    
    above_threshold = histogram_table[histogram_table['Bin Start'] >= threshold]['Count'].sum()
    print(f'Number of proteins above {threshold}: {above_threshold}')

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Generate a histogram table from a protein length distribution file.')
    parser.add_argument('input_file', type=str, help='Path to the input text file (e.g. human_length_dist.txt)')
    parser.add_argument('bin_size', type=int, help='Bin size for the histogram (e.g., 2000)')
    parser.add_argument('--threshold', type=int, default=4000, help='Threshold for counting proteins above a certain length (default: 4000)')
    parser.add_argument('--save_table', action="store_true", help='Save the table as a CSV file? Include this flag for yes, leave off for no')

    # Parse arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    generate_histogram_table(args.input_file, args.bin_size, args.threshold, args.save_table)
