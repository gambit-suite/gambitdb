import unittest
import os
import filecmp
from gambitdb.GtdbSpreadsheetParser  import GtdbSpreadsheetParser

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','gtdb_spreadsheet_parser')

class TestGtdbSpreadsheetParser(unittest.TestCase):
    """
    Test class for GtdbSpreadsheetParser.
    """

    def comp_files(self, file1, file2):
        """
    Compares two files and prints out the differences.
    Args:
      file1 (str): Path to the first file.
      file2 (str): Path to the second file.
    Returns:
      bool: True if the files are the same, False otherwise.
    Side Effects:
      Prints out the differences between the two files.
    Examples:
      >>> comp_files('file1.txt', 'file2.txt')
      The files are different. The following lines are different:
      <differences>
    """
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
        """
    Tests the generation of spreadsheets from a GTDB metadata file.
    Args:
      None
    Returns:
      None
    Side Effects:
      Generates 3 spreadsheets in the data directory.
    Examples:
      >>> test_generating_signatures()
      <No output>
    """
        self.cleanup()
        g = GtdbSpreadsheetParser(os.path.join(data_dir, 'bac120_metadata_r214'),
                                    95,
                                    5,
                                    400,
                                    False,
                                    False,
                                    1,
                                    '',
                                    os.path.join(data_dir, 'output_species.csv'),
                                    os.path.join(data_dir, 'output_genome_metadata.csv'),
                                    os.path.join(data_dir, 'output_accessions_to_download.csv'),
                                    None,
                                    False,
                                    False)
        g.generate_spreadsheets()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_genome_metadata.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_species.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_accessions_to_download.csv')))

        self.assertTrue(self.comp_files(os.path.join(data_dir, 'output_genome_metadata.csv'), os.path.join(data_dir, 'expected_genome_assembly_metadata.csv')))
        self.assertTrue(self.comp_files(os.path.join(data_dir, 'output_species.csv'), os.path.join(data_dir, 'expected_species_taxon.csv')))
        self.assertTrue(self.comp_files(os.path.join(data_dir, 'output_accessions_to_download.csv'), os.path.join(data_dir, 'expected_accession')))

        self.cleanup()

    def test_representative_genomes(self):
        self.cleanup()
        g = GtdbSpreadsheetParser(os.path.join(data_dir, 'representative_genomes.csv'),
                                    95,
                                    5,
                                    400,
                                    False,
                                    False,
                                    1,
                                    '',
                                    os.path.join(data_dir, 'output_species.csv'),
                                    os.path.join(data_dir, 'output_genome_metadata.csv'),
                                    os.path.join(data_dir, 'output_accessions_to_download.csv'),
                                    None,
                                    False,
                                    False)
        g.generate_spreadsheets()
        output_filename = os.path.join(data_dir, 'output_representative_genomes.csv')
        g.save_representative_genome_accessions_to_file( output_filename)

        self.assertTrue(os.path.exists(output_filename))
        self.assertTrue(self.comp_files(output_filename, os.path.join(data_dir, 'expected_representative_genomes.csv')))

        #self.cleanup()
        

    def cleanup(self):
        """
    Removes the generated spreadsheets from the data directory.
    Args:
      None
    Returns:
      None
    Side Effects:
      Removes 3 spreadsheets from the data directory.
    Examples:
      >>> cleanup()
      <No output>
    """
        for f in ['output_genome_metadata.csv', 'output_species.csv', 'output_accessions_to_download.csv', 'output_representative_genomes.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
