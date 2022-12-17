#!/home/syntax_surgeon/anaconda3/envs/joseml/bin/python

# Import the necessary libraries
import pubchempy as pcp
import re
import pathlib
import signal
import readchar
import sys
import time

#Define functions to be used in main()

def welcome(): # Prints a simple welcome statement
    print("\n\033[1;33m============================================\n\
++++++DATABASE ID & PROPERTY EXTRACTOR++++++\n\
====================\033[3;34mSiddharth Yadav\033[1;33m=========\033[0m\n")

    
def available_properties(): # Returns a dictionary of available properties
    # All properties in a docstring form
    properties = '''atom_stereo_count
    bond_stereo_count
    canonical_smiles
    charge
    complexity
    conformer_id_3d
    conformer_rmsd_3d
    coordinate_type
    covalent_unit_count
    defined_atom_stereo_count
    defined_bond_stereo_count
    effective_rotor_count_3d
    exact_mass
    h_bond_acceptor_count
    h_bond_donor_count
    heavy_atom_count
    inchi
    inchikey
    isomeric_smiles
    isotope_atom_count
    iupac_name
    mmff94_energy_3d
    mmff94_partial_charges_3d
    molecular_formula
    molecular_weight
    monoisotopic_mass
    multipoles_3d
    pharmacophore_features_3d
    rotatable_bond_count
    shape_fingerprint_3d
    shape_selfoverlap_3d
    tpsa
    undefined_atom_stereo_count
    undefined_bond_stereo_count
    volume_3d
    xlogp'''

    properties = [prop.strip() for prop in properties.split('\n')] # Extract properties from the docstring using a list comprehension
    properties = dict(zip(range(1,len(properties)+1), properties)) # Create a dictionary using a zip function and dict constructor
    return properties


def get_mol_file(): # Asks user for the locations of the input file in smiles format
    mol_filepath = input("Please provide the full path to the file containing the molecules in smiles format:\nEnter path: ")
    while True:
        if mol_filepath == 'q': # Quit if user types 'q'
            return False
        elif pathlib.Path(mol_filepath).is_file(): # Returns the path if it exists and is a file
            return mol_filepath
        else: # User is asked to re-enter a valid file path
            mol_filepath = input("\n\033[31m***NO SUCH FILE***\033[0m\nPlease provide the correct path to the smiles file (or type 'q' to quit):\nEnter path: ")

            
def ask_properties(property_dict): # Prompts the user to select from a list of available properties
    print("\nList of available properties:")
    for index, mol_property in property_dict.items(): # Prints the available properties with their corresponding numbers
        print(f"\t{index}\t{mol_property}")

    selected_properties = input("\nChoose the properties to include in the 'molecular_properties' file (separated by commas):\nExample: 7,18,21  (or type 'a' for all)\n\nEnter properties: ").lower()

    if selected_properties == 'a': # Returns all properties id user types 'a'
        selected_properties = list(property_dict.keys())
        return selected_properties
    else:
        while True:
            try: # Checks if the user types a non-convertible string (such as alphabets, puncuations etc.)
                selected_properties = [int(prop.strip()) for prop in selected_properties.split(",") if prop.strip()]
                
                # Checks if the numbers entered by the use are in the correct range of available property numbers
                if all([True if int(prop) in list(range(1, len(property_dict)+1)) else False for prop in selected_properties]):
                    return selected_properties
                else:
                    selected_properties = input("\n\033[31m***INVALID INPUT***\033[0m\nPlease select valid property numbers separated by commas (example: 2,11,29):\nEnter properties: ")
                
            except: # Catches the exception caused by an invalid input and resumes the loop
                selected_properties = input("\n\033[31m***INVALID INPUT***\033[0m\nPlease select valid property numbers separated by commas (example: 2,11,29):\nEnter properties: ")

                
def get_ids_from_synonyms(synonyms): # Searches, extracts and returns database ids from synonyms
    # Dictionary of pre-compiled regular expression patterns for extracting database ids
    database_patterns = {
    "chebi": re.compile(r'\bCHEBI:[\d]+\b'),
    "chembl": re.compile(r'\bCHEMBL[\d]+\b'),
    "zinc": re.compile(r'\bZINC[\d]+\b'),
    "bdbm": re.compile(r'\bBDBM[\d]+\b')
    }
    
    long_string = ' '.join(synonyms) # Concatenate synonyms to a single 'space-separated' string
    
    # Creates a dictiory of the format ('database': 'database id') using a dictionary comprehenstion
    database_ids_dict = {database: pattern.search(long_string).group() for database, pattern in database_patterns.items() if pattern.search(long_string)}
    
    return database_ids_dict


def write_db_ids(molecule, write_file, db_id_dict): # FileI/O and console printing wrapper for database ids
    write_file.write(f"PUBCHEM ID: {molecule.cid}\n") # Writes the directly obtainable PubChem ID to the file
    print(f"PUBCHEM ID: {molecule.cid}")

    for database, db_ib in db_id_dict.items(): # Iterates over the provided database id dictionary and writes them to file
        write_file.write(f"{database.upper()} ID: {db_ib}\n")
        print(f"{database.upper()} ID: {db_ib}")

    
def write_properties(molecule, write_file, selected_properties, available_properties): # FileI/O and console printing wrapper for properties
    for prop in selected_properties: # Iterates over the list of selected propeties and writes them to file
        write_file.write(f"{available_properties[prop]}: {molecule.to_dict()[available_properties[prop]]}\n")
        print(f"{available_properties[prop]}: {molecule.to_dict()[available_properties[prop]]}")
    write_file.write("\n")
    

def not_found(write_file): # Writes the error log in the file and notifies the user
    write_file.write("***NO-COMPOUND-FOUND***\n\n")
    print("\033[31m***NO-COMPOUND-FOUND***\033[0m")


def print_stats(total, errors): # Displays brief statistics for the processing
    spacer = f"\033[34m+{' '*48}+\033[0m"
    _plus = "\033[34m+  \033[0m"
    plus_ = "\033[34m  +\033[0m"
    bottom_border = f"\033[34m+ {'='*46} +\033[0m"
    
    print(f'\n\n\033[34m+ ===============\033[0mBRIEF STATISTICS\033[34m=============== +\n\
{spacer}\n\
{_plus}{"Total molecules processed":-<40}{total:->4}{plus_}\n\
{_plus}{"Total number of errors":-<40}{errors:->4}{plus_}\n\
{spacer}\n\
{_plus}\033[32m{"Success rate":-<30}{((total-errors)*100)/total:->13.2f}%{plus_}\n\
{_plus}\033[31m{"Error rate":-<30}{(errors*100)/total:->13.2f}%{plus_}\n\
{spacer}\n\
{bottom_border}\n\n')    


def handler(signum, frame): # Handles the interruption of the program via Ctrl-C
    msg = '\n\033[1;33mCtrl-C was pressed. Type "q" to quit or press any other key to continue: ' # Print message to the user
    print(msg, end="", flush=True)
    res = readchar.readchar().lower() # Prompting a single character input from the user
    if res == 'q': # Aborting program if user typed 'q'
        if count >= 2:
            print_stats(count-1, error_count) # Checking if the count has been initialized
        print('\033[1;33mAborting program...')
        time.sleep(0.5)
        sys.exit(1)
    else:
        print('\n\033[1;34mResuming process...\033[0m\n') # Resuming the script if the user typed anything other than 'q'


def main(): # Primary script
    # Defining global count variables
    global count
    global error_count
    
    # Welcome statement
    welcome() 
    
    # Get the smiles file
    mol_file = get_mol_file() 
    
    # Quit the script if user types 'q'
    if not mol_file:
        return
    
    # Display and prompt for properties
    all_properties = available_properties() 
    selected_properties = ask_properties(all_properties)
    
    # Open the smiles file for reading the molecules and open the 'molecular_properties' file to writing the output
    with open(mol_file) as read_file, open('molecular_properties.txt', 'w') as prop_file:

        # Iterating over the smiles in the smiles file
        for mol in read_file:
            current_smile = mol.strip().split()[0] # Extracting the correct smile by excluding any appended identifiers
            identifier = ' '.join(mol.strip().split()[1:]) # Extracting any identifiers if present

            try:
                # Increment the counter, print it on the console and write to the file
                count += 1 
                print(f"\n\033[32mCompound #{count}:\033[0m {identifier}")
                prop_file.write(f"Compound #{count}: {identifier}\n")
                
                # Get the compound object from PubChem using the current smile
                current_mol = pcp.get_compounds(current_smile, 'smiles')[0]
                
                try:
                    # Get database ids from the synonyms and write them to file
                    synonyms = set(current_mol.synonyms)
                    other_database_ids = get_ids_from_synonyms(synonyms)
                    write_db_ids(current_mol, prop_file, other_database_ids)
                    
                    # Get the properties for the current compound and write them to file
                    write_properties(current_mol, prop_file, selected_properties, all_properties)

                except SystemExit: # Catching the SystemExit error bound to Ctrl-C
                   break 
                
                except: # Catch any exception thrown by an empty compound object retrieved from PubChem
                    error_count += 1 # Incrementing the error count
                    not_found(prop_file) # Display error message

            except SystemExit: # Catching the SystemExit error bound to Ctrl-C
                break 
                
            except: # Catch any exception thrown if no compound object was retrieved from PubChem
                error_count += 1 # Incrementing the error count
                not_found(prop_file) # Display error message
        else:
            print_stats(count, error_count) # Display final statistics


# Defining global count variables to record progress and enable display of statistics after deliberate abortion
count = 0
error_count = 0


# Run the script if it is executed directly and NOT via import
if __name__ == "__main__":
    signal.signal(signal.SIGINT, handler) # Loading the signal handler
    main() # Running the main script
