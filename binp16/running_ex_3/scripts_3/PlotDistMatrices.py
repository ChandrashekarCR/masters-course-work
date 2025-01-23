#! ~/miniconda3/envs/ai_env/bin/python

'''
Script name - PlotDistMatrices.py

Description - This script plots a dendrogram using hireachiral clustering (UPGMA method) and a heatmap for the given matrices.


User defined function - 
1. read_dataframe: Read the dataframe which contains a matrix that is tab separated.
2. validate_file_paths: Validates if the file paths exists or not.
3. check_output_path: Checks if the output file path is valid or not.
4. get_valid_input: Gets an input from the user. Either a 'y' or 'n', and continues to ask if the user does not enter the valid characters.
5. check_the_file_name: Checks if the file name contains either Y or mt and the words alignment or identity. There is no way to know which matrix is given, when plotting.
                        Therefore, the user is prompted to enter the file names correctly for a valid output.
6. plot_dendrogram_comparison: Plot two dendrograms side-by-side for better comparison.
7. romanov_color: To color only the romonav label names with red and the rest with blue.
8. plot_heatmap_comparison: Plots two heatmaps side-by-side for better comparison.
9. plot_sngle_dendrogram: Plots a single dendrogram for the given matrix.
10. plot_single_heatmap: Plots a single heatmap for the given matrix.

Non-standard Modules:
This script requires the following non-standard Python modules to function properly:
    1. numpy: For numerical operations and matrix manipulations.
    2. pandas: For reading and processing TSV files.
    3. matplotlib: For plotting the dendrogram and heatmap visualizations.
    4. scipy: For hierarchical clustering and generating the dendrogram.
    5. seaborn: For enhanced heatmap visualizations.
    
To install these modules, you can use pip:
    pip install numpy pandas matplotlib scipy seaborn 

    Ensure that the python version you have is greater than Python 3.0.

Procedure:
    1. Input matrices (alignment and identity) and output filenames are accepted as arguments.
    2. Optional arguments allow the user to specify custom output filenames for both dendrograms and heatmaps.
    3. Verify that input files exist.
    4. Check the output paths and set defaults if paths are invalid or missing.
    5. For each input matrix, the script examines the filename to determine if is a mitochondrial (mt) or Y-chromosomal (Y) matrix.
    6. Identifies whether each matrix is an alignment or identity matrix for appropriate labeling, by just searching the words alignment or identity.
    7. Load the input matrices as DataFrames, extracting the labels and values for dendrogram and heatmap plotting.
    8. Perform hierarchical clustering on both matrices using UPGMA (average linkage) clustering.
    9. Label Romanov-associated names in red for easier identification in dendrograms.
    10. Plot heatmaps for alignment and identity matrices, displaying pairwise similarity values.
    11. Save each dendrogram and heatmap as a PNG file at a high DPI for publication-quality visuals.

Input - 
    1. Alignment score tsv file - Enter the file path for alignment matrix.
    2. Identity score tsv file - Enter the file path for identity matrix
    3. Output dendrogram 1 (optional) - Enter the file path to save the dendrogram plot 1 (alignment matrix)
    4. Output heatmap 1 (optional) - Enter the file path to save the heatmap plot 1 (alignment matrix)
    5. Output dendrogram 2 (optional) - Enter the file path to save the dendrogram plot 2 (identity matrix)
    6. Output heatmap 2 (optional) - Enter the file path to save the heatmap plot 2 (identity matrix)


Output - 
    1. 6 plots are saved as output. One individual dendrogram and heatmap for the two matrices and two comparitive plot for two matrices, i.e dendrogram and heatmap.
    2. 2 plots are saved as output. If only one matrix is provided by the user.

Usage - python plot_dendogram.py -i1 <alignment_score_file.tsv> -i2 <identity_score_file.tsv> -o1 <dendrogram_output1.png> -o2 <heatmap_output1.png> -o3 <dendrogram_ouput2> -o4 <heatmap_output2>

Version - 2.00
Date - 27/10/2024 (27th October 2024)
Name - Chandrashekar CR
'''



# Importing Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import seaborn as sns
import argparse
import time
import os
import sys
import re


# Default path for output file
DEFAULT_OUTPUT_FILE = os.getcwd()

def read_dataframe(in_tsv: str) -> tuple:

    """
    Reads a tab-separated values (TSV) file into a DataFrame and extracts the similarity matrix and labels.

    Args:
        in_tsv (str): The path to the input TSV file containing the similarity data.

    Returns:
        Tuple[numpy.ndarray, List[str]]: A tuple containing:
            - similarity_matrix_output (numpy.ndarray): The similarity matrix extracted from the DataFrame.
            - similarity_labels (List[str]): A list of labels corresponding to the rows of the similarity matrix.

    Raises:
        ValueError: If the input file cannot be read or if the expected format is not met.
    """
    try:
        in_df = pd.read_csv(in_tsv, sep='\t')
        in_df = in_df.rename(columns={in_df.columns[0]: 'names'})
        similarity_matrix_output = np.array(in_df.iloc[:, 1:])
        similarity_labels = list(in_df['names'])
    except FileNotFoundError:
        print(f'Error: The file "{in_tsv}" was not found. Please check the file path.')
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print('Error: The input file is empty. Please provide a valid TSV file.')
        sys.exit(1)
    except pd.errors.ParserError:
        print('Error: The input file is not in the correct format. Please ensure it is a valid TSV file.')
        sys.exit(1)
    except Exception as e:
        print(f'An unexpected error occurred: {e}. Please try again.')
        sys.exit(1)
    return similarity_matrix_output, similarity_labels

def check_shape_of_matrix_labels(matrix1, labels1, matrix2=None, labels2=None):
    """
    Checks if the shape of the provided matrix or matrices is compatible with the length of the provided label list(s).
    
    Args:
        matrix1 (np.ndarray): The first matrix to check.
        labels1 (list): The label list corresponding to the rows/columns of the first matrix.
        matrix2 (np.ndarray, optional): The second matrix to check. Defaults to None.
        labels2 (list, optional): The label list corresponding to the rows/columns of the second matrix. Defaults to None.
    
    Raises:
        ValueError: If the matrix shape and label length are incompatible.
    """
    try:
        # Check shape compatibility for matrix1 and labels1
        if matrix1.shape[0] != len(labels1) or matrix1.shape[1] != len(labels1):
            raise ValueError(f"Shape of matrix1 {matrix1.shape} does not match length of labels1 {len(labels1)}.")

        print("Matrix1 and labels1 are compatible.")
        print(matrix1.shape,len(labels1))
        
        # If matrix2 and labels2 are provided, check them as well
        if matrix2 is not None and labels2 is not None:
            if matrix2.shape[0] != len(labels2) or matrix2.shape[1] != len(labels2):
                raise ValueError(f"Shape of matrix2 {matrix2.shape} does not match length of labels2 {len(labels2)}.")
            
            print("Matrix2 and labels2 are compatible.")
            print(matrix2.shape,len(labels2))


    except ValueError as e:
        print(f"ValueError: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)



def validate_file_paths(*file_paths:str) -> None:
    """
    Checks if the provided file paths exist.

    Args:
        *file_paths (str): Variable number of file paths to validate.

    Raises:
        FileNotFoundError: If any of the provided paths do not exist.
    """
    for file_path in file_paths:
        try:
            if os.path.isfile(file_path):
                print(f'The file path {file_path} exists!')
            else:
                raise FileNotFoundError('Invalid file path')
        except FileNotFoundError as e:
            print(f'Expected a file but found {file_path}. {e}.')
            print(f'python plot_dendogram.py -i1 <alignment_score_file.tsv> -i2 <identity_score_file.tsv<optional>> -o1 <dendrogram_output1.png<optional>> -o2 <heatmap_output1.png<optional>> -o3 <dendrogram_ouput2<optional>> -o4 <heatmap_output2<optional>>')
            sys.exit(1)
        except TypeError as e:
            print('The file type is not valid. Please try again.')
            print('python plot_dendogram.py -i1 <alignment_score_file.tsv> -i2 <identity_score_file.tsv<optional>> -o1 <dendrogram_output1.png<optional>> -o2 <heatmap_output1.png<optional>> -o3 <dendrogram_ouput2<optional>> -o4 <heatmap_output2<optional>>')
            sys.exit(1)


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

    # Case 1: If the provided path is a directory, use 'default_output'
    if os.path.isdir(output_file_path):
        print(f"'{output_file_path}' is a directory. Using 'default_output' as the filename.")
        answer = get_valid_input('Warning: This script has multiple outputs. Consider giving your ouputs a unique name, else this will result in overwriting. You can also stick to default outputs by not giving any name. Do you still want to continue [Y/n]?: ')
        if answer == 'y':
            return os.path.join(output_file_path, 'default_output')
        else:
            print('Exiting the script.')
            sys.exit(1)

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

def get_valid_input(prompt: str) -> str:
    """
    Continuously prompts the user until they enter 'y' or 'n'.
    """
    while True:
        answer = input(prompt).strip().lower()
        if answer in ['y', 'n']:
            return answer
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

def check_the_file_name(in_file: str) -> tuple:
    in_file = os.path.basename(in_file)
    pattern_mt = r'^mt.*'
    pattern_identity = r'.*identity.*'
    pattern_alignment = r'.*alignment.*'
    pattern_Y = r'^Y.*'

    match_mt = bool(re.search(pattern_mt, in_file,re.IGNORECASE))
    match_identity = bool(re.search(pattern_identity, in_file,re.IGNORECASE))
    match_alignment = bool(re.search(pattern_alignment, in_file,re.IGNORECASE))
    match_Y = bool(re.search(pattern_Y, in_file,re.IGNORECASE))

    try:
        if match_mt and match_alignment:
            return 'mt', 'alignment'
        elif match_mt and match_identity:
            return 'mt', 'identity'
        elif match_Y and match_alignment:
            return 'Y', 'alignment'
        elif match_Y and match_identity:
            return 'Y', 'identity'
        else:
            raise ValueError("The file is not named properly")
    except ValueError as e:
        print(f'{e}. Rename the file so that it contains "mt" or "Y" at the start and either "alignment" or "identity" in the file name.')
        sys.exit(1)

def plot_the_dendogram_comparison(similarity_matrix: np.ndarray, identity_matrix: np.ndarray, y_plot_labels_alignment: list, y_plot_labels_identity: list, subplot1: str, subplot2: str, plot_title: str):
    """
    Generates and saves a side-by-side dendrogram comparison for similarity and identity matrices.

    Args:
        similarity_matrix (np.ndarray): 2D array representing the similarity matrix.
        identity_matrix (np.ndarray): 2D array representing the identity matrix.
        y_plot_labels_alignment (list): List of labels for the similarity matrix, used for the alignment dendrogram.
        y_plot_labels_identity (list): List of labels for the identity matrix, used for the identity dendrogram.
        subplot1 (str): Title descriptor for the first subplot (e.g., "alignment").
        subplot2 (str): Title descriptor for the second subplot (e.g., "identity").
        plot_title (str): Overall title for the plot.

    Raises:
        ValueError: If matrices are not 2D or labels are incompatible with matrix dimensions.
        FileNotFoundError: If the file paths are wrong.
        Exception" In case any error slips through. The exception case catches it.
    """
    try:
        # Convert matrices to condensed form for clustering
        alignment_matrix = squareform(similarity_matrix)
        identity_matrix = squareform(identity_matrix)

        # Perform hierarchical clustering on both matrices
        Z_alignment = linkage(alignment_matrix, method='average')
        Z_identity = linkage(identity_matrix, method='average')

        # Romanov names for color-coding
        romanov_names = [
            'Nicolas II Romanov', 'Alexandra Romanov', 'Olga Romanov', 
            'Tatiana Romanov', 'Maria Romanov', 'Alexei Romanov', 
            'Suspected body of Anastasia Romanov'
        ]

        # Function to assign colors
        def romanov_color(name):
            # Check if the name contains "Romanov" (case insensitive)
            if 'romanov' in ' '.join(name.split(' ')[:]).lower():
                return 'red'
            else:
                return 'blue'

        # Create subplots for side-by-side dendrograms
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

        # Dendrogram for alignment score
        dendro1 = dendrogram(
            Z_alignment, labels=y_plot_labels_alignment, orientation='left', leaf_font_size=10,
            ax=ax1
        )
        ax1.set_title(f'Dendrogram based on {subplot1.capitalize()} Score', fontsize=16, fontweight='bold')
        ax1.set_xlabel('Distance', fontsize=12)

        # Adjust y-axis tick parameters for alignment score
        ax1.tick_params(axis='y', labelsize=10)
        ax1.tick_params(axis='x', labelsize=10)

        # Color the labels manually for alignment score dendrogram
        for lbl in ax1.get_ymajorticklabels():
            lbl.set_color(romanov_color(lbl.get_text()))

        # Dendrogram for identity score
        dendro2 = dendrogram(
            Z_identity, labels=y_plot_labels_identity, orientation='right', leaf_font_size=10,
            ax=ax2
        )
        ax2.set_title(f'Dendrogram based on {subplot2.capitalize()} Score', fontsize=16, fontweight='bold')
        ax2.set_xlabel('Distance', fontsize=12)
        ax2.tick_params(axis='y', labelsize=10)
        ax2.tick_params(axis='x', labelsize=10)

        # Color the labels manually for identity score dendrogram
        for lbl in ax2.get_ymajorticklabels():
            lbl.set_color(romanov_color(lbl.get_text()))

        fig.suptitle(f'Genetic Distance Dendrogram of {plot_title.upper()}', fontsize=20, fontweight='bold')
        romanov_patch = mpatches.Patch(color='red', label='Romanov Family Members')
        other_patch = mpatches.Patch(color='blue', label='Other Samples')
        fig.legend(handles=[romanov_patch, other_patch], loc='center',bbox_to_anchor = (0.5,0.65), fontsize=12, ncol=1)
        plt.tight_layout()

        # Save the figure at high DPI for publication quality
        plt.savefig(f'comparison_dendrogram_{plot_title.lower()}.png', dpi=600)
        plt.show()

    except ValueError as e1:
        print(f"ValueError: {e1}. Ensure matrices are square and labels are compatible with matrix dimensions.")
        sys.exit(1)
    except FileNotFoundError as e2:
        print(f"FileNotFoundError: {e2}. The output path is not valid.")
        sys.exit(1)
    except Exception as e3:
        print(f"An error occurred: {e3}. Unable to generate and save the dendrogram comparison.")
        sys.exit(1)



def plot_the_heatmap_comparison(similarity_matrix: np.ndarray, identity_matrix: np.ndarray, plot_labels_alignment: list, plot_labels_identity: list, subplot1: str, subplot2: str, plot_title: str):
    """
    Generates and saves a comparison heatmap for two matrices, representing alignment and identity scores.

    Args:
        similarity_matrix (np.ndarray): 2D array representing the similarity matrix.
        identity_matrix (np.ndarray): 2D array representing the identity matrix.
        plot_labels_alignment (list): List of labels for the alignment matrix, used for x and y ticks.
        plot_labels_identity (list): List of labels for the identity matrix, used for x and y ticks.
        subplot1 (str): Title descriptor for the first subplot (e.g., "alignment").
        subplot2 (str): Title descriptor for the second subplot (e.g., "identity").
        plot_title (str): Overall title for the plot.

    Raises:
        ValueError: If matrices are not 2D or labels are incompatible with matrix dimensions.
        FileNotFoundError: If the file paths are wrong.
        Exception" In case any error slips through. The exception case catches it.
    """
    try:
        labels_alignment = plot_labels_alignment
        labels_identity = plot_labels_identity

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

        # Heatmap for alignment score
        sns.heatmap(
            similarity_matrix, ax=ax1, cmap="viridis", xticklabels=labels_alignment, yticklabels=labels_alignment,
            cbar=True, cbar_kws={'shrink': 0.5}, square=True
        )
        ax1.set_title(f'Heatmap based on {subplot1.capitalize()} Score', fontsize=16, fontweight='bold')
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=10)

        # Heatmap for identity score
        sns.heatmap(
            identity_matrix, ax=ax2, cmap="viridis", cbar=True, cbar_kws={'shrink': 0.5}, square=True,
            xticklabels=labels_identity, yticklabels=labels_identity
        )
        ax2.set_title(f'Heatmap based on {subplot2.capitalize()} Score', fontsize=16, fontweight='bold')
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax2.set_yticklabels(ax2.get_yticklabels(), fontsize=10)

        # Overall title for both heatmaps
        plt.suptitle(f'Genetic Distance Heatmap of {plot_title.upper()}', fontsize=20, fontweight='bold')
        plt.tight_layout()  # Adjust layout to fit the title rect=[0, 0, 1, 0.96]

        # Save the figure at high DPI for publication quality
        plt.savefig(f'comparison_heatmap_{plot_title.lower()}.png', dpi=600)
        plt.show()

    except ValueError as e1:
        print(f"ValueError: {e1}. Please ensure 'similarity_matrix' and 'identity_matrix' are 2D with compatible label dimensions.")
        sys.exit(1)
    except FileNotFoundError as e2:
        print(f"FileNotFoundError: {e2}. The output path is not valid.")
        sys.exit(1)
    except Exception as e3:
        print(f"An error occurred: {e3}. Unable to generate and save the heatmap comparison.")
        sys.exit(1)


    

def plot_single_dendrogram(matrix: np.ndarray, y_plot_labels: list, title: str, score_type: str, output_filename: str):
    """
    Generates and saves a dendrogram based on the provided matrix and labels.

    Args:
        matrix (np.ndarray): 2D array representing the similarity or distance matrix to be clustered.
        y_plot_labels (list): List of labels for each row in the matrix, corresponding to leaves in the dendrogram.
        title (str): Title for the dendrogram plot, used in the figure's title.
        score_type (str): A description of the type of score (e.g., "identity" or "alignment") to include in the plot's title.
        output_filename (str): The filename (including path) where the generated dendrogram image will be saved.

    Raises:
        ValueError: If matrices are not 2D or labels are incompatible with matrix dimensions.
        FileNotFoundError: If the file paths are wrong.
        Exception" In case any error slips through. The exception case catches it.
    """
    try:
        # Convert matrix to a condensed form for hierarchical clustering
        condensed_matrix = squareform(matrix)
        Z = linkage(condensed_matrix, method='average')

        # Romanov names for color-coding
        romanov_names = [
            'Nicolas II Romanov', 'Alexandra Romanov', 'Olga Romanov', 
            'Tatiana Romanov', 'Maria Romanov', 'Alexei Romanov', 
            'Suspected body of Anastasia Romanov'
        ]

        # Function to assign colors
        def romanov_color(name):
            if 'romanov' in ' '.join(name.split()).lower():
                return 'red'
            else:
                return 'blue'

        # Create the dendrogram plot
        fig, ax = plt.subplots(figsize=(16, 8))
        dendro = dendrogram(
            Z, labels=y_plot_labels, orientation='right', leaf_font_size=10, ax=ax
        )
        ax.set_xlabel('Distance', fontsize=12)
        ax.tick_params(axis='y', labelsize=10)
        ax.tick_params(axis='x', labelsize=10)

        # Color the labels manually
        for lbl in ax.get_ymajorticklabels():
            lbl.set_color(romanov_color(lbl.get_text()))

        fig.suptitle(f'Genetic Distance Dendrogram of {title.upper()} based on {score_type.capitalize()} Score',fontsize=18,fontweight='bold')
        romanov_patch = mpatches.Patch(color='red', label='Romanov Family Members')
        other_patch = mpatches.Patch(color='blue', label='Other People')
        fig.legend(handles=[romanov_patch, other_patch], loc='upper right', bbox_to_anchor=(0.98,0.94), fontsize=12, ncol=1)
        plt.tight_layout()

        # Save the figure at high DPI for publication quality
        plt.savefig(f'{output_filename}', dpi=600)
        plt.show()

    except ValueError as e1:
        print(f"ValueError: {e1}. Please check that 'matrix' is 2D and matches 'y_plot_labels' dimensions.")
        sys.exit(1)
    except FileNotFoundError as e2:
        print(f"FileNotFoundError: {e2}. The output path '{output_filename}' is not valid.")
        sys.exit(1)
    except Exception as e3:
        print(f"An error occurred: {e3}. Unable to generate and save the dendrogram.")
        sys.exit(1)


def plot_single_heatmap(matrix: np.ndarray, plot_labels: list, title: str, score_type: str,output_filename) -> None:
    """
    Generates and saves a heatmap based on the provided matrix and labels.

    Args:
        matrix (np.ndarray): 2D array representing the similarity or distance matrix to be visualized.
        plot_labels (list): List of labels for the x and y axes, corresponding to the rows/columns of the matrix.
        title (str): Title for the heatmap plot, used in the figure's title.
        score_type (str): A description of the type of score (e.g., "identity" or "alignment") to include in the plot's title.
        output_filename (str): The filename (including path) where the generated heatmap image will be saved.

    Raises:
        ValueError: If matrices are not 2D or labels are incompatible with matrix dimensions.
        FileNotFoundError: If the file paths are wrong.
        Exception" In case any error slips through. The exception case catches it.
    """
    try:
        fig, ax = plt.subplots(figsize=(16, 8))

        # Plot the heatmap
        sns.heatmap(
            matrix, ax=ax, cmap="viridis", xticklabels=plot_labels, yticklabels=plot_labels,
            cbar=True, square=False, annot=False, fmt=".2f", annot_kws={"size": 7}
        )
        #ax.set_title(f'Heatmap Based on {score_type} Score', fontsize=16, fontweight='bold')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)

        # Overall title
        fig.suptitle(f'Genetic Distance Heatmap of {title.upper()} based on {score_type.capitalize()} Score',fontsize=18,fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust to fit title without overlap

        # Save the figure at high DPI for publication quality
        plt.savefig(f'{output_filename}', dpi=600)
        plt.show()
    except ValueError as e1:
        print(f"ValueError: {e1}. Please check that 'matrix' is 2D and matches 'plot_labels' dimensions.")
        sys.exit(1)
    except FileNotFoundError as e2:
        print(f"FileNotFoundError: {e2}. The output path '{output_filename}' is not valid.")
        sys.exit(1)
    except Exception as e3:
        print(f"An error occurred: {e3}. Unable to generate and save the heatmap.")
        sys.exit(1)



if __name__ == '__main__':


    # Define command-line arguments
    parser = argparse.ArgumentParser(
        description='A program that reads a matrix and outputs a dendrogram and heatmap.',
        add_help=True,
        usage='python plot_dendogram.py -i1 <alignment_score_file.tsv> -i2 <identity_score_file.tsv<optional>> -o1 <dendrogram_output1.png<optional>> -o2 <heatmap_output1.png<optional>> -o3 <dendrogram_ouput2<optional>> -o4 <heatmap_output2<optional>>'
    )
    parser.add_argument('-i1', '--input1', dest='alignment_score', help='Enter the alignment similarity matrix file.', required=True)
    parser.add_argument('-i2', '--input2', dest='identity_score', help='Enter the identity similarity matrix file.', nargs='?')
    parser.add_argument('-o1', '--output_dendrogram1', dest='dendrogram_plot_input1', help='Enter the file name to save the dendrogram plot for the first matrix.', nargs='?',default=DEFAULT_OUTPUT_FILE+'/default_output_dendogram_1')
    parser.add_argument('-o2', '--output_heatmap1', dest='heatmap_plot_input1', help='Enter the file name to save the heatmap plot for the first matrix.', nargs='?',default=DEFAULT_OUTPUT_FILE+'/default_output_heatmap_1')
    parser.add_argument('-o3', '--output_dendrogram2', dest='dendrogram_plot_input2', help='Enter the file name to save the dendrogram plot for the second matrix.', nargs='?',default=DEFAULT_OUTPUT_FILE+'/default_output_dendogram_2')
    parser.add_argument('-o4', '--output_heatmap2', dest='heatmap_plot_input2', help='Enter the file name to save the heatmap plot for the second matrix.', nargs='?',default=DEFAULT_OUTPUT_FILE+'/default_output_heatmap_2')
    args = parser.parse_args()

    
    # Check if the output files are valid filepaths and names
    dendrogram_plot1 = check_output_path(args.dendrogram_plot_input1) + '.png'
    heatmap_plot1 = check_output_path(args.heatmap_plot_input1) + '.png'
    dendrogram_plot2 = check_output_path(args.dendrogram_plot_input2) + '.png'
    heatmap_plot2 = check_output_path(args.heatmap_plot_input2) + '.png'


    # Check if the second input exists
    validate_file_paths(args.alignment_score)
    if args.identity_score is not None:
        # Check if the file path exists
        validate_file_paths(args.identity_score)

        # Parse the file names for title and subplot titles.
        title1, subplot1 = check_the_file_name(args.alignment_score)
        title2, subplot2 = check_the_file_name(args.identity_score)
        if title1 == title2:
            title = title1
        else:
            title = '_'.join([title1,title2])
        print(f'Title: "{title.upper()}", Subplot Titles: "{subplot1.capitalize()}" and "{subplot2.capitalize()}"')
        
        # Read the dataframe
        alignment_matrix, alignment_labels = read_dataframe(args.alignment_score)
        identity_matrix, identity_labels = read_dataframe(args.identity_score)

        # Check the shape of the matrix
        check_shape_of_matrix_labels(matrix1=alignment_matrix,labels1=alignment_labels,matrix2=identity_matrix,labels2=identity_labels)
        
        # Plot the dendogram comparison
        plot_the_dendogram_comparison(alignment_matrix, identity_matrix, alignment_labels,identity_labels ,subplot1, subplot2, title)
        
        # Plot the heatmap for comparison
        plot_the_heatmap_comparison(alignment_matrix, identity_matrix, alignment_labels,identity_labels, subplot1, subplot2, title)
        
        # Plot a single dendrogram for alignment matrix
        plot_single_dendrogram(alignment_matrix, alignment_labels, title=title1, score_type=subplot1,output_filename=dendrogram_plot1)
        
        # Plot a single heatmap for alignment matrix
        plot_single_heatmap(alignment_matrix, alignment_labels, title1, subplot1,output_filename=heatmap_plot1)

        # Plot a single dendogram for identity matrix
        plot_single_dendrogram(identity_matrix, identity_labels, title=title2, score_type=subplot2,output_filename=dendrogram_plot2)

        # Plot a single heatmap for identity matrix
        plot_single_heatmap(identity_matrix,identity_labels,title2,subplot2,output_filename=heatmap_plot2)


        print('Done!')
        sys.exit(1)
    else:
        
        # Parse the file names for the title and subplot title
        title,subplot1 = check_the_file_name(args.alignment_score)
        print(f'Title: "{title.upper()}", Subplot Title: "{subplot1.capitalize()}"')

        # Read the dataframe
        alignment_matrix,alignment_labels = read_dataframe(args.alignment_score)

        # Check the shape of the matrix
        check_shape_of_matrix_labels(matrix1=alignment_matrix,labels1=alignment_labels)

        # Plot a single dendogram and single heatmap
        plot_single_dendrogram(alignment_matrix,alignment_labels,title,subplot1,output_filename=dendrogram_plot1)
        plot_single_heatmap(alignment_matrix,alignment_labels,title,subplot1,output_filename=heatmap_plot1)
        print('Done')
        sys.exit(1)
    