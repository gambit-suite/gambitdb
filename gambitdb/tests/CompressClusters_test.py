# tests for the CompressClusters class
import unittest
import os
import pandas as pd
from gambitdb.CompressClusters import CompressClusters

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','compress_clusters')

class TestCompressClusters(unittest.TestCase):
    """
    Test class for CompressClusters class.
    """

    def test_compress_clusters(self):
        pairwise_filename = os.path.join(data_dir, 'pairwise_distances.csv')
        representative_genomes_filename = os.path.join(data_dir, 'representative_genomes.csv')

        c = CompressClusters(pairwise_filename, 'output', 0.2, representative_genomes_filename, False)
        samples_to_remove, num_samples = c.compress()
        self.assertEqual(num_samples, 13)
        self.assertEqual(len(samples_to_remove), 7)
        self.assertNotIn('GCA_10', samples_to_remove)
