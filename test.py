# Importing unittest & pubchempy modules
import unittest
import pubchempy as pcp


# Defining a TestCase class
class TestCaseForPropertyExtractor(unittest.TestCase):
    def setUp(self): # Setup code to avoid reinitialization of expensive connection objects 
        self.smile = "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4" # Sample smile for atorvastatin
        self.current_mol = pcp.get_compounds(self.smile, 'smiles')[0] # Initilizing the compound object

    def test_for_name(self): # Testing for the string 'atorvastatin' in the list of retrieved synonyms
        self.assertTrue("".join(self.current_mol.synonyms).lower().count("atorvastatin") >= 1)
        print("\n---> Testing for name")

    def test_for_mol_formula(self): # Testing for molecular formula
        self.assertEqual(self.current_mol.molecular_formula, "C33H35FN2O5")
        print("\n---> Testing for molecular weight")

    def test_for_inchi_key(self): # Testing for InChi-key
        self.assertEqual(self.current_mol.inchikey, "XUKUURHRXDUEBC-UHFFFAOYSA-N")
        print("\n---> Testing for InChi-key")

        
# Run the script if it is executed directly and NOT via import
if __name__ == '__main__':
    unittest.main()
