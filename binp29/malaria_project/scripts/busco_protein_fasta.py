import os
import pandas as pd
import numpy as np
import argparse
import subprocess


sequence_list = {"species": [],
                 "busco_id": [],
                 "status": [],
                 "sequence": []}

def read_tsv(data_in, species_name):
    with open(data_in, "r") as f_in:
        for line in f_in:
            if not line.startswith("#"):
                status = line.strip().split()[1]
                if status == "Complete" or status == "Duplicated":
                    sequence = line.strip().split()[2]
                    sequence_list["sequence"].append(sequence)
                    sequence_list["busco_id"].append(line.strip().split()[0])
                    sequence_list["status"].append(line.strip().split()[1])
                    sequence_list["species"].append(species_name)


import os

def extract_gene_sequences(fasta_file, gene_id):
    """
    Extract sequences corresponding to the given gene ID from the FASTA file.
    
    :param fasta_file: Path to the FASTA file.
    :param gene_id: The gene ID to extract sequences for.
    :return: Dictionary with headers and sequences.
    """
    sequences = {"header": [], "seq": []}
    with open(fasta_file, 'r') as file:
        flag = False
        current_seq = ""
        for line in file:
            if line.startswith('>'):  # Header line
                if flag:  # Save the previous sequence
                    sequences["seq"].append(current_seq)
                    current_seq = ""  # Reset for the next sequence
                # Check if the header starts with the given gene ID
                if line[1:].split()[0] == gene_id:
                    sequences["header"].append(line.strip())  # Save the header
                    flag = True
                else:
                    flag = False
            elif flag:  # Sequence line and flag is True
                current_seq += line.strip()  # Concatenate sequence lines
        if flag:  # Save the last sequence
            sequences["seq"].append(current_seq)
    return sequences


def create_fasta_for_busco(filtered_df, fasta_dir, output_dir):
    """
    Create a fasta file for each BUSCO ID with sequences from different species.
    
    :param filtered_df: DataFrame containing BUSCO IDs, species, and gene IDs.
    :param fasta_dir: Directory containing species-specific FASTA files.
    :param output_dir: Directory to save the output FASTA files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create output directory if it doesn't exist

    for index, row in filtered_df.iterrows():
        busco_id = row['busco_id']
        species_list = row['species']  
        gene_ids = row['sequence']  

        # Ensure species_list and gene_ids are the same length (they should be matched pairs)
        if len(species_list) != len(gene_ids):
            print(f"Skipping BUSCO {busco_id} due to mismatched species-gene pairs.")
            continue

        output_file = os.path.join(output_dir, f"{busco_id}.faa")
        with open(output_file, 'w') as outfile:
            # Iterate over corresponding species-gene pairs
            for species, gene in zip(species_list, gene_ids):
                fasta_file = os.path.join(fasta_dir, f"{species}.faa")

                if os.path.exists(fasta_file):  # Ensure the FASTA file exists
                    sequence_data = extract_gene_sequences(fasta_file, gene)
                    # Write headers and sequences to the output file
                    for header, seq in zip(sequence_data["header"], sequence_data["seq"]):
                        outfile.write(f">{species}\n{seq}\n")

        print(f"Created FASTA file for BUSCO ID: {busco_id}")


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run filtering of busco proteins")
    parser.add_argument("--full_table", dest="fulltable", required=True, help="Enter the full table path")
    parser.add_argument("--fasta_directory", dest="fasta_dir", required=True, help="Enter the fasta directory")
    parser.add_argument("--output_dir",dest="output_dir",required=True,help="Enter the output directory")

    # Store parsed arguments
    args = parser.parse_args()

    species_list = ["Ht", "Pb", "Pc", "Pf", "Pk", "Pv", "Py", "Tg"]

    fulltable = args.fulltable
    fasta_dir = args.fasta_dir
    output_dir = args.output_dir

    # Read in data for each species
    for species in species_list:
        path_to_table = f"{fulltable}{species}/run_apicomplexa_odb12/full_table.tsv"
        read_tsv(path_to_table, species)

    # Create DataFrame from the sequence list
    seqeunce_df = pd.DataFrame.from_dict(sequence_list)
    df_exploded = seqeunce_df.explode(["species","sequence"])
    df_filtered = df_exploded.drop_duplicates(subset=["busco_id","species"],keep='first')
    df_final = df_filtered.groupby("busco_id").agg({
        "species":lambda x: list(x),
        "sequence": lambda x:list(x)
    }).reset_index()

    df_final = df_final[df_final["species"].apply(lambda x: set(x) == set(species_list))]

    create_fasta_for_busco(df_final,fasta_dir,output_dir)

    print('Done!')