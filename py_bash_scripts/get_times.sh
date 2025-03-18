#!/bin/bash
# get_times.sh
# Usage: ./get_times.sh '<input_glob>' <output_file>
#
# Example:
#   ./get_times.sh '4484222_*.out' timings_200_299.txt

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 '<input_glob>' <output_file>"
    echo "The input glob must be quoted!"
    exit 1
fi

INPUT_GLOB=$1
OUTPUT_FILE=$2

grep -hE '^(Job Array ID:|Array index:|Total )' $INPUT_GLOB > $OUTPUT_FILE
echo "Extracted timings saved to $OUTPUT_FILE"
