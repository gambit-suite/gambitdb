# tests for the Diameters class
import unittest
import os
import pandas
from gambitdb.Curate import Curate

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','curate')

# Curate( options.species_taxon_filename, 
#         options.genome_assembly_metadata, 
#         options.species_to_remove, 
#         options.accessions_to_remove, 
#         options.species_taxon_output_filename,
#         options.genome_assembly_metadata_output_filename,
#         options.accession_removed_output_filename,
#         options.species_removed_output_filename,
#         options.minimum_ngenomes,
#         options.debug,
#         options.verbose)

class TestCurate(unittest.TestCase):
    def test_remove_species_with_zero_diameter(self):
        c = Curate( None, None, None, None, None, None, None, None, None, False, False)
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
        c = Curate( None, None, None, None, None, None, None, None, 1, False, False)
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
        c = Curate( None, None, os.path.join(data_dir, 'species_to_remove'), None, None, None, None, None, None, False, False)
        species = pandas.read_csv(os.path.join(data_dir, 'species.csv'), index_col=False)
        species = species.set_index('species_taxid')

        filtered_species = c.remove_species_using_input_file(species)
        n_rows_after = filtered_species.shape[0]
        self.assertEqual(n_rows_after, 3)
        self.assertEqual(c.species_taxonids_removed, [3, 5])
        self.assertEqual(c.species_removed,['Orange teal','G2'])

    def test_remove_accessions_using_input_file(self):
        c = Curate( None, None, None, os.path.join(data_dir, 'accessions_to_remove'), None, None, None, None, None, False, False)
        accessions = pandas.read_csv(os.path.join(data_dir, 'assembly_metadata.csv'), index_col=False)
        accessions = accessions.set_index('assembly_accession')

        # get the number of accessions before
        n_rows_before = accessions.shape[0]
        self.assertEqual(n_rows_before, 3)
        filtered_accessions = c.remove_accessions_using_input_file(accessions)
        n_rows_after = filtered_accessions.shape[0]
        self.assertEqual(n_rows_after, 2)
        self.assertEqual(c.accessions_removed, ['C789'])

    def test_filter_files(self):
        self.cleanup()
        c = Curate( os.path.join(data_dir, 'species.csv'), 
                    os.path.join(data_dir, 'assembly_metadata.csv'), 
                    os.path.join(data_dir, 'species_to_remove'), 
                    os.path.join(data_dir, 'accessions_to_remove'), 
                    os.path.join(data_dir, 'species_taxon_output.csv'),
                    os.path.join(data_dir, 'assembly_metadata_output.csv'),
                    os.path.join(data_dir, 'accessions_removed'), 
                    os.path.join(data_dir, 'species_removed'), 
                    1, False, False)
        c.filter_spreadsheets_and_output_new_files()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'species_taxon_output.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'assembly_metadata_output.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'accessions_removed')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'species_removed')))
        self.cleanup()

    def cleanup(self):
        for f in ['species_taxon_output.csv', 'assembly_metadata_output.csv', 'accessions_removed', 'species_removed']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
        
