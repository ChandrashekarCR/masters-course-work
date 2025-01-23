# FASTA Functional Annotation Script

## Author
Chandrashekar CR
Version: 1.00 | Date: 12th October 2024

## Description
This Python script processes a FASTA file and a BLASTx annotation file to generate a new FASTA file with annotated sequence headers. If no matching keys are found between the two files, the original FASTA file is output unchanged.

## Usage
python malaria.py <fasta_file> <blastx_file> <output_file_or_directory>
Eg - python malaria.py malaria.fna malaria.blastx.tab ../results/output.txt

- <fasta_file>: Input FASTA file path
- <blastx_file>: Input BLASTx file path
- <output_file_or_directory>: Path to output file or directory. If a directory is specified, the output will be named default_output.txt.

## Workflow
- 1. Input Validation:

    Verifies argument count and file paths.
    Checks file format (FASTA and BLASTx).

- 2. Reading Input Files:

    Stores FASTA sequences and BLASTx annotations in dictionaries.

- 3. Annotation Process:

    If matching keys are found, headers are annotated with protein names from BLASTx.
    If no matches, the original FASTA is copied as-is to the output.

## Error Handling
- Invalid paths: Script exits if paths are incorrect.
- Incorrect format: Prompts for correct input files if mismatched.
- No matching keys: Writes the original FASTA if files are unrelated.

## Dependencies
- Python Standard Libraries:

    os, sys

- No external dependencies.

