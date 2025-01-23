#! ~/miniconda3/bin/python

'''
Script name - FastaAligner.py

Description - This script calculates alignment scores between sequences from a FASTA file (Only DNA sequences). 
It performs pairwise comparisons to determine matches, gaps, transitions, and transversions between sequences. 
The scoring can be customized via a parameter file, or default values are used if not provided.

User defined function - 
1. validate_file_paths: Validates if the provided file paths exist.
2. read_fasta: Reads a FASTA file and stores its sequences in a dictionary.
3. calculate_score: Computes alignment scores between two sequences.
4. helper_calculate_scores: Facilitates all-against-all sequence comparisons.
5. calculate_percentage: Calculates the percentage of gaps or identity matches.
6. validate_nucleotides_in_sequence: Validates if sequences contain only valid nucleotide characters.
7. valid_parameters: Reads and validates the scoring parameters from a file.
8. get_valid_input: Repeatedly prompts the user for valid input.
9. check_output_path: Ensures the provided output path is valid and adjusts if necessary.

Non standard module - None (All used modules are part of the Python Standard Library.)

Procedure:
1. Parse command-line arguments for the input FASTA file, parameter file (optional), and output file path (optional).
2. Validate the input files and ensure sequences contain only valid nucleotides.
3. Read the sequences from the FASTA file into a dictionary.
4. If a parameter file is provided, validate and parse the custom scoring parameters; otherwise, use default parameters.
5. Perform pairwise comparisons of sequences, calculating alignment scores and tracking matches, gaps, and mutations.
6. Calculate percentages of gaps and identity matches for each comparison.
7. Print the results of the alignment comparisons. ?
8. Save the results to the specified output file path.


Input - 
1. FASTA file: Contains nucleotide sequences for alignment.
2. Parameter file (optional): Custom scoring parameters for match, transition, transversion, and gap penalties.
3. Output file (optional): Path to save the output; defaults to 'default_output.txt'.

Output - 
Printed and saved alignment scores for each sequence pair, including:
- Total score
- Number of matches, gaps, and mutations
- Percentage of gaps and identity matches

Usage - python <python script> <FASTA file> <parameter file (optional)> <output file (optional)>
Eg. python FastaAligner.py ../data/score.exampl.fna ../data/parameters.txt ../results/alignmet_score.txt

Version - 2.00
Date - 22/10/2024 (22nd October 2024)
Name - Chandrashekar CR
'''



# Importing libraries
import os
import argparse
import sys
import re
import textwrap


# Global variables
# Expected parameter keys in the correct order
EXPECTED_KEYS = ["match", "transition", "transversion", "gap"]
# Valid nucleotides
VALID_NUCLEOTIDES = {'A','T','G','C','N','a','t','g','c','n','-'}
# Types muatations in AGTC nucleotides
MUTATION_DICT = {
    'transition': ['AG', 'GA', 'CT', 'TC'],
    'transversion': ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']
                }
# Default scores
DEFAULT_SCORING_PARAMETERS = {'match':1,'transition':-1,'transversion':-2,'gap':-1}
# Default path for output file
DEFAULT_OUTPUT_FILE = os.getcwd()


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
            print(f'Usage: python <script_name_here>.py -i <fasta_file> -p <parameter file <optional>> -o <output file<optional>>')
            sys.exit(1)

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


def calculate_score(seq1: str, seq2: str,scoring_dict:dict) -> float:
    """
    Reads two sequences and the scoring matrix used. The score is calculated between the two sequences by considering
    a single nucleotide from each sequence for a given index position. The score is calculated based on a match, transistion,
    transversion or gap and the result is a float value.

    Args:
        seq1 (str): Sequence 1 for comparison
        seq2 (str): Sequence 2 for comparison
        scoring_dict (dict): Paramters for scoring the matches, mutations and gaps.

    Returns:
        float: The alignment score.

    """
    score, gap_count,identity_count,curated_length = 0,0,0,0  # Initialize score to 0
    

    for i in range(len(seq1)):

        try:

            # If both the sequences have a gap at the same postiion.
            if seq1[i] == '-' and seq2[i] == '-':
                continue

            curated_length+=1
        
            # Check if both positions are a gap, then add the gap-score to the score variable. Also increase the gap count by 1.
            if seq1[i] == '-' or seq2[i] == '-':
                #print('gap')
                score += scoring_dict['gap']
                gap_count+=1

            # Check if both nucleotides match, then add the match-score to the score variable. Also increase the identity count by 1.
            elif seq1[i] == seq2[i]:
                #print('match')
                score += scoring_dict['match']
                identity_count+=1

            # If both nucleotides are a mismatch, then if they are a transition, add the transition-score to the score varible.
            # Else add the transversion-score tot the score variable.
            elif f'{seq1[i]}{seq2[i]}' in MUTATION_DICT['transition']:
                #print('transition')
                score += scoring_dict['transition']

            elif f'{seq1[i]}{seq2[i]}' in MUTATION_DICT['transversion']:
                #print('transversion')
                score += scoring_dict['transversion']

            else:
                # Optional: Handle unknown cases (if needed)
                raise ValueError('unknown mutation type')
        # This is just in case the validation of the sequences is not done properly and something might have slipped through the quality check.
        except ValueError as e:
            print(f'{e}')
    try:
        if curated_length ==0:
            raise ZeroDivisionError('The sequence is either full of gaps or the wrong file is being parsed.')
    except ZeroDivisionError as e:
        print(f'{e}. Please try again...')
        sys.exit(1)
    return score, gap_count, identity_count, curated_length

def helper_calculate_scores(fasta_dict:dict) -> dict:
    """
    Read the sequences using two loops, which allows us to perfrom an all against all comparison.

    Args:
        fasta_dict (dict): Values stored in a dictionary are unpacked for comparions
    
    Returns:
        dict: The function returns a nested dictionary in which the key is a tuple, which comprises of two sequences that
        are being compared and a value associated with it is three dictionaries score, gap and indentity. In these other 
        three dictionaries the keys are the name of the parameters 'score', 'gap' and 'identity' and the values associated with it
        are the values of these parameters.

    """
    results = {}
    seq_ids = list(fasta_dict.keys())
    curated_length_list = []

    # Two for loops are used which helps in pefroming pairwise comparision between sequences.
    for i in range(len(seq_ids)):
        for j in range(i+1,len(seq_ids)):
            id1,id2 = seq_ids[i],seq_ids[j]
            seq1, seq2 = fasta_dict[id1].upper(), fasta_dict[id2].upper() # Sequences that are in lower cases are converted to upper case.

            # The function 'calculate_score' is called for a pair of sequences.
            score, gap, identity, curated_length = calculate_score(seq1,seq2,parameters_file)
            curated_length_list.append(curated_length) 

            results[(id1,id2)] = {
                'score':score,
                'gap':gap,
                'identity':identity
            }

    return results, curated_length_list

def calculate_percentage(val:float,total:float) -> float:
    """
    Calculates the percentage of the values given.

    Args:
        val (flaot): The value corresponding to the percentage to be calculated on.
        total (float): The value of the total to calculate the percentage on.
    
    Returns:
        float: Returns a float of calculated percentage.

    """
    # Calculates the percentage and returns a float value.
    return (100 * val)/total


def validate_nucleotides_in_sequence(fasta_file: str) -> None:
    """
    Validates the content of a FASTA file to ensure it is correctly formatted and
    contains only valid nucleotide sequences.

    Args:
        fasta_file (str): Path to the FASTA file.

    Raises:
        FileNotFoundError: If the file is not found or if the content does not follow the FASTA format.
        ValueError: If the sequence contains invalid characters outside the valid nucleotide set.
    """
    try:
        with open(fasta_file, 'r') as f_in:
            # Ensure the file contains at least one valid header
            valid_fasta = False  

            for line in f_in:
                line = line.strip()
                
                # Check if the line is a header
                if line.startswith('>'):
                    # Found a valid header, then set the 'valid_fasta' variable to True.
                    valid_fasta = True  
                else:
                    # Check if the sequence line contains only valid nucleotides
                    if not set(line).issubset(VALID_NUCLEOTIDES):
                        raise ValueError(f"Invalid character(s) in sequence: {line}")

            # If no valid headers were found, raise an error
            if not valid_fasta:
                raise FileNotFoundError("No valid FASTA headers found. Please provide a correct FASTA file.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except UnicodeDecodeError as e:
        print('This file format is not readable. Please try again...')
        sys.exit(1)

    print("The FASTA file is valid.")


def valid_parameters(parameters_file_in: str) -> dict:
    """
    Reads parameters from a file and ensures there are exactly four key-value pairs.
    If not, prompts the user to provide default values or raises an error.

    Args:
        parameters_file_in (str): Path to the parameters_file specified by the user.

    Returns: 
        dict: Extracts the parameters set by the user using regex and stores the inputs in a dicitonary.
        The key is the type of match and the value associated with it is a score.

    """
    # Read the file and apply regex
    try:
        with open(parameters_file_in, 'r') as f_in:
            read_params = f_in.read().lower()

        # Regex to capture key-value pairs
        # Regex Explanantion - (\w+) -> A group of letter characters that occur once or more.
        #                    - \s* -> Spaces zero times or more.
        #                    - [=:>\-]? -> Any of these characters.
        #                    - \s* -> Spaces zero times or more.
        #                    - ([+-]?\d+[\.]?[\d+]?) -> A plus or - either once or none of the times, followed by a digit, once or more, followed by a decimal point once or none, then digits once or none.
        pattern = re.compile(r'(\w+)\s*[=:>\-]?\s*([+-]?\d+[\.]?[\d+]?)')  # These were my initial regex tries -([A-Za-z]+)?\s*[;=:-\>+\s]?\s*([+-]?\d+) #(^\w*)[\s=\-:>]*([+-]?\d*$) #(^\w+)\s*([+-]?\d+$)

        matches = re.findall(pattern, read_params)
        # Check if we have exactly 4 parameters
        try:
            # If the number of matches by regex findall methods is less than or equal to 4, then create two lists -:
            # 'found'  which stores the values/names that were also in the 'EXPECTED_KEYS'
            # 'missing' which stores the values/names that were not in the 'EXPECTED_KEYS'
            if len(matches) <= 4:
                found = [key for key, _ in matches if key in EXPECTED_KEYS]
                missing = [key for key in EXPECTED_KEYS if key not in found]

                # If the found variable does not contain anything matching the EXPECTED_KEYS, then use the default scoring values.
                if len(found)==0:
                    print(f'Found no user input.')
                    print(f'Using default parameters {DEFAULT_SCORING_PARAMETERS}')
                    user_params_score = DEFAULT_SCORING_PARAMETERS
                    return user_params_score
                
                # If the found variable contains all the EXPECTED_KEYS, then we use the user defined values.
                if len(found)==4:
                    print('All the parameters are valid.')
                    user_params_score = dict(matches)
                    user_params_score = {param_type:float(param_val) for param_type, param_val in user_params_score.items()}
                    return user_params_score

                # If all the cases above fail. There some parameters that were missing, then ask the user to use the default values.
                # If yes, then use the default scoring parameters, else handle the case and exit the program.
                print(f"Found: {', '.join(found)} but missing: {', '.join(missing)}")
                answer = get_valid_input(f'4 scores needed but got {len(found)}. Do you want to use default values for the missing ones? [Y/n] ').strip().lower()
                
                if answer == 'y':
                    # Fill missing parameters with default values
                    for i in missing:
                        matches.append((f'{i}',f'{DEFAULT_SCORING_PARAMETERS[i]}'))
                elif answer == 'n':
                    raise SyntaxError('Need exactly 4 values to parse.')
            else:
                raise ValueError

        except SyntaxError as e:
            print(f'{e}')
            sys.exit(1)
        except ValueError:
            print(f'Too many values to unpack. Please follow the specified format.')
            sys.exit(1)

        
        user_params_score = dict(matches)
        user_params_score = {param_type:int(param_val) for param_type, param_val in user_params_score.items()}
        return user_params_score
    # If the file is not in readable format. This handles that error.
    except UnicodeDecodeError as e:
        print('The file format is not readable. Please try again...')
        sys.exit(1)
    # In case, any errors have slipped through that are syntax based. This handles that case as well.
    except SyntaxError as e:
        print('The parameters file contains an invalid format. Please try again...')
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


if __name__ == '__main__':

    # Parse command line arguements
    parser = argparse.ArgumentParser(
    prog='FastaAligner.py',
    formatter_class=argparse.RawTextHelpFormatter,
    usage='python FastaAligner.py -i <input FASTA file> -p <parameter file <optional>> -o <output file<optional>>',
    description='This is a program to calculate the alignment scores for aligned sequences.'
    )

    # FASTA file as input.
    parser.add_argument('-i','--input_fasta_file',dest='fasta_file', help='Enter the path of the fasta file.')

    # Parameters file.
    parser.add_argument('-p','--parameters_file',
        dest='parameters_file',
        help=textwrap.dedent('''\
        Enter the path for the parameters to be used.
        Use the specified format for writing the output file -:
                             
        match=(your value here)
        transition=(your value here)
        transversion=(your value here)
        gap=(your value here)
                             
        '''),
        default=DEFAULT_SCORING_PARAMETERS,
        nargs='?'
    )

    # Output file.
    parser.add_argument('-o','--output_file',dest='output_file',
        help='Enter the path of the output file and the name.',
        default=DEFAULT_OUTPUT_FILE,
        nargs='?'
    )

    args = parser.parse_args()
    
    # Check if the parameters file is not a dictionary.
    if type(args.parameters_file) != type(DEFAULT_SCORING_PARAMETERS):
        parameters_file = args.parameters_file
        # Check if the file paths are correct.
        validate_file_paths(parameters_file)   
        # Call the function 'valid_parameters' which uses regular expressions to parse the parameters file that is defined by the user. 
        parameters_file = valid_parameters(parameters_file)
        # Print the parameters being used for the user's convenience.
        print(f'The parameters being used are {parameters_file}')
    else:
        print('Using default parameters for scoring.')
        print(f'The default parameters are {DEFAULT_SCORING_PARAMETERS}')
        parameters_file = DEFAULT_SCORING_PARAMETERS

    
    # Check if the FASTA_file has a valid file path.
    validate_file_paths(args.fasta_file)
    # Check if the fasta file is valid
    validate_nucleotides_in_sequence(args.fasta_file)

    # Check if the output file path is valid
    validated_output_file_path = check_output_path(args.output_file)
    print(validated_output_file_path)

    # Call the function 'read_fasta'. It returns the fasta file in the form of a dictionary.
    fasta_dict = read_fasta(args.fasta_file)
    # Call the function 'helper_calculate_score'. It returns a dictionary that stores all the information, like alignment scores, gaps, identities etc. 
    results, new_length = helper_calculate_scores(fasta_dict)
    count = 0



    # Write the result into the output file.
    with open(validated_output_file_path,'w') as f_out:
        for (id1,id2),scores in results.items():
            f_out.write(f"{id1.replace('>','')}-{id2.replace('>','')}: "
            f"Identity: {scores['identity']}/{new_length[count]} "
            f"({calculate_percentage(scores['identity'], new_length[count]):.1f}%), "
            f"Gaps: {scores['gap']}/{new_length[count]} "
            f"({calculate_percentage(scores['gap'], new_length[count]):.1f}%), "
            f"Score={scores['score']:.5f}")
            f_out.write('\n')
            count+=1

print('Done!')