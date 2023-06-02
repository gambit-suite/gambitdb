import unittest
import os
import filecmp
from gambitdb.GtdbSpreadsheetParser  import GtdbSpreadsheetParser

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','gtdb_spreadsheet_parser')

class TestGtdbSpreadsheetParser(unittest.TestCase):

    def comp_files(self, file1, file2):
        result = filecmp.cmp(file1, file2)

        if result:
            return True
        else:
            with open(file1, 'r') as file1_fh, open(file2, 'r') as file2_fh:
                diff = set(file1_fh).difference(file2_fh)
                print("The files are different. The following lines are different:")
                for line in diff:
                    print(line.strip())
            return False    
    
    def test_generating_signatures(self):
        self.cleanup()
        g = GtdbSpreadsheetParser(os.path.join(data_dir, 'bac120_metadata_r214'),
                                    95,
                                    5,
                                    400,
                                    False,
                                    os.path.join(data_dir, 'output_species.csv'),
                                    os.path.join(data_dir, 'output_genome_metadata.csv'),
                                    os.path.join(data_dir, 'output_accessions_to_download.csv'),
                                    False,
                                    False)
        g.generate_spreadsheets()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_genome_metadata.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_species.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_accessions_to_download.csv')))

        self.comp_files(os.path.join(data_dir, 'output_genome_metadata.csv'), os.path.join(data_dir, 'expected_genome_assembly_metadata.csv'))
        self.comp_files(os.path.join(data_dir, 'output_species.csv'), os.path.join(data_dir, 'expected_species_taxon.csv'))
        self.comp_files(os.path.join(data_dir, 'output_accessions_to_download.csv'), os.path.join(data_dir, 'expected_accession'))  

        self.cleanup()

    def cleanup(self):
        for f in ['output_genome_metadata.csv', 'output_species.csv', 'output_accessions_to_download.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
