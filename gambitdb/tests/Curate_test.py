# tests for the Curate class
import unittest
import os
import pandas
from gambitdb.Curate import Curate

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','curate')

class TestCurate(unittest.TestCase):
    """
    Tests for the Curate class.
    """
    def test_remove_species_with_zero_diameter(self):
        """
    Tests the remove_species_with_zero_diameter method of the Curate class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Sets c.species_taxonids_removed and c.species_removed.
    Examples:
      >>> test_remove_species_with_zero_diameter()
      None
    """
        c = Curate( None, None, None, None, None, None, None, None, None, None, None, 1, 0.99, 0.7, 1, False, False)
        species = pandas.read_csv(os.path.join(data_dir, 'species.csv'), index_col=False)
        species = species.set_index('species_taxid')

        # get the number of rows in species
        n_rows_before = species.shape[0]
        self.assertEqual(n_rows_before, 5)
        filtered_species = c.remove_species_with_zero_diameter(species)
        n_rows_after = filtered_species.shape[0]
        self.assertEqual(n_rows_after, 3)
        self.assertEqual(c.species_taxonids_removed, [2, 6])
        self.assertEqual(c.species_removed,['Red white','G3'])


    def test_remove_species_with_fewer_than_n_genomes(self):
        """
    Tests the remove_species_with_fewer_than_n_genomes method of the Curate class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Sets c.species_taxonids_removed and c.species_removed.
    Examples:
      >>> test_remove_species_with_fewer_than_n_genomes()
      None
    """
        c = Curate( None, None, None, None, None, None, None, None, None, None, 1, 1, 0.99, 0.7, 1, False, False)
        species = pandas.read_csv(os.path.join(data_dir, 'species.csv'), index_col=False)
        species = species.set_index('species_taxid')

        # get the number of rows in species
        n_rows_before = species.shape[0]
        self.assertEqual(n_rows_before, 5)
        filtered_species = c.remove_species_with_fewer_than_n_genomes(species)
        n_rows_after = filtered_species.shape[0]
        self.assertEqual(n_rows_after, 3)
        self.assertEqual(c.species_taxonids_removed, [5, 6])
        self.assertEqual(c.species_removed,['G2','G3'])

    def test_remove_species_using_input_file(self):
        """
    Tests the remove_species_using_input_file method of the Curate class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Sets c.species_taxonids_removed and c.species_removed.
    Examples:
      >>> test_remove_species_using_input_file()
      None
    """
        c = Curate( None, None, None, None, os.path.join(data_dir, 'species_to_remove'), None, None, None, None, None, None, 1, 0.99, 0.7, 1, False, False)
        species = pandas.read_csv(os.path.join(data_dir, 'species.csv'), index_col=False)
        species = species.set_index('species_taxid')

        filtered_species = c.remove_species_using_input_file(species)
        n_rows_after = filtered_species.shape[0]
        self.assertEqual(n_rows_after, 3)
        self.assertEqual(c.species_taxonids_removed, [3, 5])
        self.assertEqual(c.species_removed,['Orange teal','G2'])

    def test_remove_accessions_using_input_file(self):
        """
    Tests the remove_accessions_using_input_file method of the Curate class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Sets c.accessions_removed.
    Examples:
      >>> test_remove_accessions_using_input_file()
      None
    """
        c = Curate( None, None, None, None, None, os.path.join(data_dir, 'accessions_to_remove'), None, None, None, None, None, 1, 0.99, 0.7, 1, False, False)
        accessions = pandas.read_csv(os.path.join(data_dir, 'assembly_metadata.csv'), index_col=False)
        accessions = accessions.set_index('assembly_accession')

        # get the number of accessions before
        n_rows_before = accessions.shape[0]
        self.assertEqual(n_rows_before, 3)
        filtered_accessions = c.remove_accessions_using_input_file(accessions)
        n_rows_after = filtered_accessions.shape[0]
        self.assertEqual(n_rows_after, 2)
        self.assertEqual(c.accessions_removed, ['C789'])

    def test_real_dataset(self):
        """
    Tests the Curate class with a real dataset.
    Args:
      None
    Returns:
      None
    Side Effects:
      Creates and removes files in the data_dir.
    Examples:
      >>> test_real_dataset()
    """
        self.cleanup()
        c = Curate( os.path.join(data_dir, 'ec_species'),
                    os.path.join(data_dir, 'ec_genome_assembly'),
                    None, 
                    os.path.join(data_dir, 'ec_pairwise'),
                    None,
                    None,                     
                    os.path.join(data_dir, 'species_taxon_output.csv'),
                    os.path.join(data_dir, 'assembly_metadata_output.csv'),
                    os.path.join(data_dir, 'accessions_removed'), 
                    os.path.join(data_dir, 'species_removed'), 
                    1, 1, 0.99, 0.7, 1, False, False)
        c.filter_spreadsheets_and_output_new_files()

        # The input assembly metadata file shouldnt be equal to the output assembly metadata file because we have removed some rows
        self.assertNotEqual(self.file_len(os.path.join(data_dir, 'ec_genome_assembly')), self.file_len(os.path.join(data_dir, 'assembly_metadata_output.csv')))
        self.assertEqual(self.file_len(os.path.join(data_dir, 'assembly_metadata_output.csv')), 87)

        # The input species file shouldnt be equal to the output species file because we have removed some rows
        self.assertNotEqual(self.file_len(os.path.join(data_dir, 'ec_species')), self.file_len(os.path.join(data_dir, 'species_taxon_output.csv')))
        self.assertEqual(self.file_len(os.path.join(data_dir, 'species_taxon_output.csv')), 18)
        self.cleanup()

    def test_filter_files(self):
        """
    Tests the Curate class with a filter files.
    Args:
      None
    Returns:
      None
    Side Effects:
      Creates and removes files in the data_dir.
    Examples:
      >>> test_filter_files()
    """
        self.cleanup()
        c = Curate( os.path.join(data_dir, 'species.csv'), 
                    os.path.join(data_dir, 'assembly_metadata.csv'), 
                    None,
                    os.path.join(data_dir, 'pairwise.csv'),
                    os.path.join(data_dir, 'species_to_remove'), 
                    os.path.join(data_dir, 'accessions_to_remove'), 
                    os.path.join(data_dir, 'species_taxon_output.csv'),
                    os.path.join(data_dir, 'assembly_metadata_output.csv'),
                    os.path.join(data_dir, 'accessions_removed'), 
                    os.path.join(data_dir, 'species_removed'), 
                    1, 1, 0.99, 0.7, 1, False, False)
        c.filter_spreadsheets_and_output_new_files()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'species_taxon_output.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'assembly_metadata_output.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'accessions_removed')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'species_removed')))
        self.cleanup()

    # method to return the number of lines in a file
    def file_len(self, fname):
        """
    Returns the number of lines in a file.
    Args:
      fname (str): The path to the file.
    Returns:
      int: The number of lines in the file.
    Examples:
      >>> file_len('data/curate/species.csv')
      18
    """
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def cleanup(self):
        """
    Removes files in the data_dir.
    Args:
      None
    Returns:
      None
    Side Effects:
      Removes files in the data_dir.
    Examples:
      >>> cleanup()
    """
        for f in ['species_taxon_output.csv', 'assembly_metadata_output.csv', 'accessions_removed', 'species_removed']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
