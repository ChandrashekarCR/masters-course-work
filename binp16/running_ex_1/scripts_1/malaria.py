#! ~/miniconda3/bin/python

'''
Script name - malaria.py

Description - This program takes three inputs from the user: a FASTA file, a BLASTx file, and output file name or directory. 
The goal is to generate a new FASTA file where each sequence header is joined with protein name that are not 'null' (if available) from the BLASTx file.
A key or identifier is defined as the part of the sequence header in the FASTA file, which is also present in the blastx file. This identifer is the
most crucial part of the program that helps in creating the output file by parsing the FASTA file and the blastx file. 
If there are no matching keys (identifiers) between the FASTA and blastx files, the script outputs the original FASTA file as it is. 

User defined function - None
Non standard module - None

Procedure:
1) Input validation and error handling: Ensures correct number of arguments, valid file paths, and correct file formats.
2) Reads and stores data from the FASTA file and blastx file into dictionaries.
3) Compares the keys (identifiers) from both files.
4) Writes a new FASTA file with protein names in the headers if matching keys are found. If not even a single match exists, then it writes the original FASTA file.
 

Input - FASTA file, Blastx file and output file or output folder.
Output - Fuctional annotated FASTA file.

Usage - python <python program file> <fasta file> <blastx file> <output file> # Note the ouput file can be a directory or a file.
Eg. python ../scripts/malaria.py malaria.fna malaria.blastx.tab ../results/output.txt

Version - 2.00
Date - 12/10/2024 (12th October 2024)
Name - Chandrashekar CR

'''


# Importing libraries.
import os
import sys

# Get the current working directory.
wd = os.getcwd()

# 1) Handling potential errors.

# a) Too many arguements error. This handles cases like, if the user provides arguements more than expected or less than expected.
if len(sys.argv)==4:
        print(f"Expected 4 arguments, and got {len(sys.argv)} arguements.")
else:
    print(f'Expected 4 arguements but got {len(sys.argv)}. Please try again....')
    sys.exit(1)

# Store the all the arguments in the list.
no_args = [sys.argv[0],sys.argv[1],sys.argv[2],sys.argv[3]]
# sys.argv[0] -> Python script
# sys.argv[1] -> FASTA file
# sys.argv[2] -> blastx file
# sys.argv[3] -> output file 

# b) File not found error. This handles cases, that check if the file path is a valid file path.

# Loop over all the arguments and check if the file paths are valid.
for i in no_args:
    # Special case. The fourth arguement (output file) is handled differently. The logic is mentioned below.
    '''Case 1: Check if the path corresponds to a directory. Then use the 'default_ouput.txt' as the file name for the ouput file to be saved.
    Case 2: First check if the directory is valid. If the user has specified output name for the file then use that file name for the ouput.
    Case 3: If both the cases above are not satisfied, then handle the error.''' 
    if i == sys.argv[3]:
        try:
            directory = os.path.join(wd,i)
            # Case 1
            if os.path.isdir(directory):
                print(f'{os.path.join(wd,directory)} is a valid directory for output file.')
                sys.argv[3] = os.path.join(directory, 'default_output.txt')
            # Case 2
            elif os.path.isdir(os.path.dirname(i)):
                print(f'Provided path is valid for a file: {os.path.join(wd, i)}')
                print(f'Keeping the file name: {os.path.basename(i)} as output file.')
            # Case 3
            else:
                raise FileNotFoundError('The output file path is invalid.')
        except FileNotFoundError as e:
            print(f'{e} Please try again...')
            sys.exit(1)

    # The other arguements (Python script, fasta file and blastx file) are checked to see if they are valid file paths.
    else:    
        while True:
            file_path = os.path.join(wd,i)
            try:
                # Check if the file path is valid.
                if os.path.isfile(file_path):
                    print(f'The file path {file_path} exists!')
                    break
                else:
                    raise FileNotFoundError('Invalid file path')
            except FileNotFoundError as e:
                print(f'Expected a file but found {i}. {e}. Please try again....')
                sys.exit(1)
            

# c) Check if it is a Fasta file.

# The FASTA file format always starts with header, wherein the first character of the line is '>'.
# The code below uses this property of the FASTA file and checks if the line starts with '>'. 
# If the line does not start with '>', then the error is handled by prompting the user to enter the correct FASTA file.

try:
    with open(f'{sys.argv[1]}','r') as f_in:
        valid_seq = ''
        for line in f_in:
            # Check if the line starts with '>'.
            if line.startswith('>'):
                print('This is a valid FASTA file.')
                break
            else:
                raise FileNotFoundError('Wrong file is being parsed.')
    
     
except FileNotFoundError as e:
    print(f'{e} Please enter the correct FASTA file and try again....')
    sys.exit(1)


# d) Check if it is a blastx file

# The blastx format contains values that are separated by tab. The first line in the blastx format contains column names.
# The first column is supposed to contain the queryName, that starts with '#' and hitDescription.
'''Case 1: If the line starts with '#' then we split the line based on tab separated values and store them in the list. 
        We then remove the '#' from the first element of the list.
Case 2: If the line does NOT start with '#' then we split the line based on the tab separated values and store them in the list.
Case 3: If both cases are not satisfied then we handle the error. This is because either the file is in the wrong format or it is the wrong file.'''


cols_to_check = ['queryName','hitDescription'] # Check if these two columns are present in the blastx file.

try:
    with open(sys.argv[2],'r') as f_blast:
        for line in f_blast:
            # Case 1
            if line.startswith('#'):
                cols_in_file = line.strip().split('\t')
                cols_in_file[0] = cols_in_file[0][1:]
                break
            # Case 2
            elif not line.startswith('#'):
                cols_in_file = line.strip().split('\t')
                break
            # Case 3
            else:
                raise NameError('Wrong file is being parsed.')
except NameError as e:
    print(f'{e} A blastx file is required. Please try again....')
    sys.exit(1)

# This is to check if the list, cols_to_check, is a subset of the list, cols_in_file.
# If it is a subset, i.e all the elements in the cols_to_check are present in the cols_in_file list, then the blastx file is valid. 
# Otherwise the error is handled by prompting the user to enter the correct format of the blastx file and exiting the program.
try:
    if set(cols_to_check).intersection(set(cols_in_file)) == set(cols_to_check): 
        print('This is a valid blastx file.')
    else:
        raise FileNotFoundError('Wrong file is being parsed. The columns required for analysis are not present. A proper blastx format file is required.')

except FileNotFoundError as e:
    print(f'{e}. Please try again ....')
    sys.exit(1)

# Main code begins here.

# 2) Reading in the FASTA file.
'''
After checking the FASTA files thoroughly. The FASTA file is parsed. An empty dictionary (fasta_dict) is created that takes in a id and sequence i.e, key-value pair respectively.
The file is read line by line. An empty string called fasta_seq is created everytime a line starts with '>'. This strategy can handle multiline fasta as well as single line fasta files.
'''
fasta_file =  os.path.join(wd,sys.argv[1]) 
fasta_dict = {}
with open(f'{fasta_file}','r') as f1:
    for line in f1:
        if line.startswith('>'):
            fasta_id = line.strip().split('\t')[0].replace('>','') # 1st index (zeroth in python language) of the fasta header.
            fasta_seq = '' # Creating and empty string to store the sequence for the corresponding fasta file.
        else:
            fasta_seq += line.strip() # Removes both leading and trailing white spaces.
        
        if fasta_id and fasta_seq:
            fasta_dict[fasta_id] = fasta_seq # Storing the fasta_id (identifier) and the sequence as a key value pair.
            
# 3) Reading in the blastx annotated file.
'''
After checking blastx files thoroughly. The blastx file is parsed. An empty dictionary (blastx_dict) is created that takes in an id and the function i.e, key-value pair respectively.
Two cases are handled,which are -
Case 1: The line starts with '#'. The line is split according to tab separators and the elements are stored in a list (blastx_cols). The '#' from the first element is removed. 
        The index for the queryName (identifier or key) and hitDescription (function) is stored in blastx_id and blastx_function variables.
Case 2: If the line does not start with '#'. Then line is split based on tab separator and the elements are stored in a list (blastx_cols).
        The index for the queryName (indentifier or key) and hitDescription (function) is stored in the blastx_id and blastx_function variables.
The id and the function are stored in the blastx_dict.
'''
blastx_file = os.path.join(wd,sys.argv[2]) 
blastx_dict = {}
with open (f'{blastx_file}','r') as f2:
    for line in f2:
        # Case 1
        if line.startswith('#'):
            blastx_cols = line.strip().split('\t')
            blastx_cols[0] = blastx_cols[0][1:]
            blastx_id = blastx_cols.index('queryName')
            blastx_function = blastx_cols.index('hitDescription')
        # Case 2
        else:
            blastx_cols = line.strip().split('\t')
            blastx_id = blastx_cols.index('queryName')
            blastx_function = blastx_cols.index('hitDescription')
        break
    for line in f2:
        # Additional information. 1st index in the list contains the id that is also present in the FASTA file. 9th index in the list contains the function. May or may not be true for every blastx file.
        blastx_dict[(line.strip().split('\t')[blastx_id])]=(line.strip().split('\t')[blastx_function])


# This is mainly done to handle a case where the user provides a file with a correct fasta file and a correct blastx file, but they are not related. 
# This is mainly used as a flag later on in the script to check if there are common keys are not.
common_keys = set(fasta_dict.keys()).intersection(set(blastx_dict.keys())) # To find if there are any common keys between the two dictionaries.

# 4) Writing results into a new file.
'''
The original FASTA file is opened in read-mode and the ouput file is opened in write-mode simultaneously. 
We check for two different cases -
Case 1: If the common_keys consists of common ids. Read the FASTA file header line and store the id/identifier (key_id) and remaining header (rest_header). 
Case 2: If the common keys does not even has a single common id. Write the FASTA file as it is to the new ouput file.
'''
results_file = os.path.join(wd,sys.argv[3]) 
with open(f'{fasta_file}','r') as fi, open(f'{results_file}','w') as fo:
    for line in fi:
        # Case 1
        if len(common_keys) != 0:
            if line.startswith('>'):
                key_id = line.strip().split('\t')[0].replace('>','') # This gets the key ID from the fasta file.
                rest_header = '\t'.join(line.strip().split('\t')[1:]) # This stores everything else apart from the key ID.

                '''Check if the key_id (id or identifier) is present in the keys of the blastx_dict and the value associated with it is not 'null'.
                Then we check if the key_id is also present in the fasta_dict.
                Finally, we format the output into a FASTA format by combining all the outputs.'''

                if key_id in blastx_dict.keys() and blastx_dict[key_id]!='null': 
                    if key_id in fasta_dict.keys():
                        seq_header = f'>{key_id}\t{rest_header}\tprotein={blastx_dict[key_id]}'
                        seq_seq = fasta_dict[key_id]
                        fo.write(seq_header)
                        fo.write('\n')
                        fo.write(seq_seq)
                        fo.write('\n')
        # Case 2
        else:
            print('The FASTA file and the blastx file are not compatible. Writing the FASTA file as the output.')
            fo.write(line.strip())
            fo.write('\n')
            fo.writelines(fi.readlines())
                

print('Done!')          
                       
            
            