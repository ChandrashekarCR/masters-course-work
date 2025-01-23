# FastaAligner.py

## Description

FastaAligner is a Python script designed to calculate alignment scores between DNA sequences provided in a FASTA file. The script performs pairwise comparisons between sequences to determine matches, gaps, transitions, and transversions. Users can customize scoring parameters or use default values.

The script is highly customizable, includes comprehensive input validation, and is written to handle edge cases, ensuring robust and accurate results.

## Features

- Validates input FASTA files for correct format and valid nucleotide sequences.
- Performs all-against-all pairwise comparisons of sequences.
- Calculates alignment scores, number of matches, gaps, and mutations.
- Supports customizable scoring parameters via a parameter file.
- Provides detailed output, including percentages of gaps and identity matches.
- Saves results to a specified output file or defaults to `default_output.txt`.

## Requirements

This script uses only Python's standard library and does not require any non-standard modules.

## Usage

### Command-Line Arguments

```bash
python FastaAligner.py <FASTA file> <parameter file (optional)> <output file (optional)>
