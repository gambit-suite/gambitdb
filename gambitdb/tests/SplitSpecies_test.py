# tests for the SplitSpecies class
import unittest
import os
import pandas
from gambitdb.SplitSpecies import SplitSpecies

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','split_species')

class TestSplitSpecies(unittest.TestCase):
    def test_read_files(self):
        ss = SplitSpecies(os.path.join(data_dir, 'species_taxon.csv'),
                          os.path.join(data_dir, 'assembly_metadata.csv'), 
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          0.7,
                          1, False)
        species, assemblies, pairwise = ss.read_input_files()
        self.assertEqual(species.shape[0], 3)
        self.assertEqual(assemblies.shape[0], 13)
        self.assertEqual(pairwise.shape[0], 13)

    def test_split_high_diameter_species(self):
        ss = SplitSpecies(os.path.join(data_dir, 'species_taxon.csv'),
                          os.path.join(data_dir, 'assembly_metadata.csv'), 
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          0.7,
                          1, False)
        species, a, p = ss.read_input_files()
        self.assertEqual(species.shape[0], 3)
        high_diameter_species  = ss.filter_high_diameter_species(species)
        self.assertEqual(high_diameter_species.shape[0], 1)

    def test_split_species(self):
        ss = SplitSpecies(os.path.join(data_dir, 'species_taxon.csv'),
                          os.path.join(data_dir, 'assembly_metadata.csv'), 
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          0.7,
                          1, False)
        species, a, p = ss.split_high_diameter_species()
