#!/usr/bin/env python3

import argparse
import os
import sys
import re
from collections import Counter



def read_fasta(fasta_file:str) -> dict:
    """
    Reads a FASTA file and stores its contents in a dictionary.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        dict: A dictionary where keys are FASTA headers and values are sequences.
    """
    fasta_dict = {}
    current_id = None
    current_seq = []
    with open(fasta_file,'r') as f_in:
        for line in f_in:
            line = line.strip()
            # The FASTA file is parsed. An empty dictionary (fasta_dict) is created that takes in a id and sequence i.e, key-value pair respectively.
            # The file is read line by line. An empty string called current_seq is created everytime a line starts with '>'. 
            if line.startswith('>'):
                if current_id:
                    # Save the previous sequence before starting the new one
                    fasta_dict[current_id]=''.join(current_seq) # This strategy can handle multiline fasta as well as single line fasta files.
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save the last sequence 
        if current_id:
            fasta_dict[current_id] = ''.join(current_seq) 
    return fasta_dict


def filtering_the_reads(fasta_dict, gc_content, scaffold_length):
    filtered_scaffolds = {}

    for head, seq in fasta_dict.items():
        # Extract the length from the header using regex
        threshold_length = re.search(r'length=(\d+)', head)
        if threshold_length:
            scaffold_len = int(threshold_length.group(1))
        else:
            scaffold_len = len(seq)  # If no length is provided in the header, use sequence length

        # Calculate the GC content of the sequence
        nucleotide_counts = Counter(seq.upper())  # Convert to uppercase to avoid case sensitivity
        gc_count = nucleotide_counts['G'] + nucleotide_counts['C']
        gc_percentage = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0  # Avoid division by zero

        # Apply the filtering conditions
        if scaffold_len >= scaffold_length and gc_percentage <= gc_content:
            filtered_scaffolds[head] = seq

    return filtered_scaffolds

       
def check_output_path(output_file_path: str) -> str:
    """
    Validates the output file path and handles the following scenarios:
    
    1. If the provided path is a directory, the output will be saved as 'default_output.txt' in that directory.
    2. If the path includes a valid directory and filename, use the specified filename.
    3. If the user provides only a filename without a directory, use the current working directory.
    4. If the directory is invalid, raise an error and exit.

    Args:
        output_file_path (str): The path where the output file should be saved.

    Returns:
        str: The final output file path.

    Raises:
        FileNotFoundError: If the provided directory does not exist.
    """

    # Normalize the path and handle relative paths
    output_file_path = os.path.normpath(output_file_path)
    print(f"Normalized path: {output_file_path}")

    # Get the current working directory
    wd = os.getcwd()

    # Case 1: If the provided path is a directory, use 'default_output.txt'
    if os.path.isdir(output_file_path):
        print(f"'{output_file_path}' is a directory. Using 'default_output.txt' as the filename.")
        return os.path.join(output_file_path, 'default_output.txt')

    # Extract the directory and filename from the provided path
    directory = os.path.dirname(output_file_path)
    filename = os.path.basename(output_file_path)

    # Case 2: If only a filename is provided without a directory, use the current working directory
    if not directory:
        print(f"No directory provided. Using the current directory with filename '{filename}'.")
        return os.path.join(wd, filename)

    # Case 3: If the directory exists, use the provided filename
    if os.path.isdir(directory):
        print(f"Valid directory found. Using '{filename}' as the output file.")
        return output_file_path

    # Case 4: If the directory is not valid, raise an error
    try:
        raise FileNotFoundError(f"The directory '{directory}' is invalid.")
    except FileNotFoundError as e:
        print(f"{e} Please try again...")
        sys.exit(1)


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run filtering on the scaffolds.")
    parser.add_argument("--fasta_file", dest="fasta_file", required=True, help="Enter the fasta file")
    parser.add_argument("--output_file", dest="output_file", required=True, help="Store the output file", nargs="?")
    parser.add_argument("--gc_content", dest="gc_content", type=int, required=True, help="Enter the threshold values of GC content. This will remove all scaffolds that are greater than the threshold value.")
    parser.add_argument("--scaffold_len", dest="scaffold_len", type=int, required=True, help="Enter the threshold value for scaffold length. This will remove all the scaffolds that are shorter than the threshold value.")
    
    # Store parsed arguments
    args = parser.parse_args()

    # Now you can access the parsed arguments
    fasta_file = args.fasta_file
    output_file = args.output_file
    gc_content_threshold = args.gc_content
    scaffold_length_threshold = args.scaffold_len

    # Check if the output file path is valid
    validated_output_file_path = check_output_path(args.output_file)
    print(validated_output_file_path)


    fasta_dict = read_fasta(fasta_file=fasta_file)
    filtered_scaffolds = filtering_the_reads(fasta_dict,28,3000)
    
    # Write the result into the output file.
    with open(validated_output_file_path,'w') as f_out:
        for header, seq in filtered_scaffolds.items():
            f_out.write(header)
            f_out.write('\n')
            f_out.write(seq.upper())
            f_out.write('\n')
         