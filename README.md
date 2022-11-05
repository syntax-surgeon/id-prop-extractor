# id-prop-extractor
> A simiple database IDs and biochemical properties extractor for compounds in SMILES format using the pubchempy API.
> #### Author: *Siddharth Yadav (syntax-surgeon)*

### Description of files:
### main.py
* Contains the main program which prompts the user to provide a path to a smiles files and extracts PubChem, ZINC, CHEBI, CHEMBL & BDBM (BindingDB) ids
* REGEX is utilized to extract database ids from the synonyms column in the associated PubChem profile
* Several properties included in the PubChem profile of the compound can also be extracted
* The data is written to a file named 'molecular_properties.txt' in the directory from where the script was run
* Compounds not found are written as '\*\*\*NO-COMPOUND-FOUND\*\*\*' in the 'molecular_properties.txt' file
* Ctrl+C can be utilized to quit/pause the script

### test.py
* Tests for the appropriate and intended functionality of the pubchempy API
* Uses three test cases based on the name, molecular formula and Inchi-Key of the drug 'Atorvastatin'
