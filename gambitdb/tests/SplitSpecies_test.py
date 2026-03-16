# tests for the SplitSpecies class
import unittest
import os
import pandas as pd
from gambitdb.SplitSpecies import SplitSpecies

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','split_species')

class TestSplitSpecies(unittest.TestCase):
    """
    Test class for SplitSpecies class.
    """
    def test_read_files(self):
        """
    Tests the read_files method of the SplitSpecies class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Asserts that the species, assemblies, and pairwise dataframes have the expected shapes.
    Examples:
      >>> test_read_files()
      None
    """
        species = pd.read_csv(os.path.join(data_dir, 'species_taxon.csv'), index_col=False).set_index('species_taxid')
        genome_metadata = pd.read_csv(os.path.join(data_dir, 'assembly_metadata.csv')).set_index('assembly_accession')

        ss = SplitSpecies(species,
                          genome_metadata,
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          [],
                          0.7,
                          1, 'average', False)
        species, assemblies, pairwise = ss.read_input_files()
        self.assertEqual(species.shape[0], 3)
        self.assertEqual(assemblies.shape[0], 13)
        self.assertEqual(pairwise.shape[0], 13)

    def test_split_high_diameter_species(self):
        """
    Tests the filter_high_diameter_species method of the SplitSpecies class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Asserts that the high_diameter_species dataframe has the expected shape.
    Examples:
      >>> test_split_high_diameter_species()
      None
    """
        species = pd.read_csv(os.path.join(data_dir, 'species_taxon.csv'), index_col=False).set_index('species_taxid')
        genome_metadata = pd.read_csv(os.path.join(data_dir, 'assembly_metadata.csv')).set_index('assembly_accession')

        ss = SplitSpecies(species,
                          genome_metadata,
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          [],
                          0.7,
                          1, 'average', False)
        species, a, p = ss.read_input_files()
        self.assertEqual(species.shape[0], 3)
        high_diameter_species  = ss.filter_high_diameter_species(species)
        self.assertEqual(high_diameter_species.shape[0], 1)

    def test_split_species(self):
        """
    Tests the split_high_diameter_species method of the SplitSpecies class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Asserts that the species, genome_metadata, and accessions_removed dataframes have the expected shapes.
    Examples:
      >>> test_split_species()
      None
    """
        species = pd.read_csv(os.path.join(data_dir, 'species_taxon.csv'), index_col=False).set_index('species_taxid')
        genome_metadata = pd.read_csv(os.path.join(data_dir, 'assembly_metadata.csv')).set_index('assembly_accession')

        ss = SplitSpecies(species,
                          genome_metadata,
                          os.path.join(data_dir, 'pairwise_distances.csv'),
                          [],
                          0.7,
                          1, 'average', False)
        s, g, accessions_removed = ss.split_high_diameter_species()
        self.assertEqual(s.shape[0], 5)
        self.assertEqual(g.shape[0], 13)
        self.assertEqual(len(accessions_removed), 0)

    def test_two_genome_high_diameter_species_removed(self):
        """
        Tests that a species with only 2 genomes and diameter > 0.7 is removed entirely.
        """
        species = pd.read_csv(os.path.join(data_dir, 'species_taxon_two_genome.csv'), index_col=False).set_index('species_taxid')
        genome_metadata = pd.read_csv(os.path.join(data_dir, 'assembly_metadata_two_genome.csv')).set_index('assembly_accession')

        ss = SplitSpecies(species,
                          genome_metadata,
                          os.path.join(data_dir, 'pairwise_distances_two_genome.csv'),
                          [],
                          0.7,
                          1, 'average', False)
        s, g, accessions_removed = ss.split_high_diameter_species()
        # Yellow black (2 genomes, diameter 0.9) should be removed entirely
        # Remaining: Yellow white + Yellow genus = 2 species rows, diameter set to 0.0 for Yellow black
        self.assertNotIn('Yellow black subspecies', ' '.join(s['name'].tolist()))
        # GCA_1 and GCA_2 should be in accessions_removed
        self.assertIn('GCA_1', accessions_removed)
        self.assertIn('GCA_2', accessions_removed)
        self.assertEqual(len(accessions_removed), 2)

    def test_singleton_outliers_removed_species_kept(self):
        """
        Tests that singleton outliers are removed and the species is kept with recalculated diameter.
        Red black has 5 genomes: GCA_1,2,3 are close (0.1), GCA_4 is an outlier (0.9 from all),
        GCA_5 is another outlier (0.8 from cluster, 0.9 from GCA_4).
        Splitting produces: cluster of 3 (GCA_1,2,3), singleton GCA_4, singleton GCA_5.
        Singletons are removed, species kept with diameter recalculated from GCA_1,2,3 only.
        """
        species = pd.read_csv(os.path.join(data_dir, 'species_taxon_singleton.csv'), index_col=False).set_index('species_taxid')
        genome_metadata = pd.read_csv(os.path.join(data_dir, 'assembly_metadata_singleton.csv')).set_index('assembly_accession')

        ss = SplitSpecies(species,
                          genome_metadata,
                          os.path.join(data_dir, 'pairwise_distances_singleton.csv'),
                          [],
                          0.7,
                          1, 'average', False)
        s, g, accessions_removed = ss.split_high_diameter_species()
        # GCA_4 and GCA_5 should be removed as singletons
        self.assertIn('GCA_4', accessions_removed)
        self.assertIn('GCA_5', accessions_removed)
        self.assertEqual(len(accessions_removed), 2)
        # Species should still exist (not subspeciated) with recalculated diameter
        red_black = s[s['name'] == 'Red black']
        self.assertEqual(len(red_black), 1)
        # Diameter should be 0.1 (max distance among GCA_1,2,3)
        self.assertAlmostEqual(float(red_black['diameter'].iloc[0]), 0.1, places=2)
        # ngenomes should be 3 (outliers removed)
        self.assertEqual(int(red_black['ngenomes'].iloc[0]), 3)
        # No subspecies should have been created
        self.assertNotIn('subspecies', ' '.join(s['name'].tolist()))

    def test_all_singletons_species_removed(self):
        """
        Tests that a species where all genomes are singletons (all >0.7 apart) is removed entirely.
        """
        species = pd.DataFrame({
            'name': ['All apart', 'Good species', 'Parent'],
            'rank': ['species', 'species', 'genus'],
            'parent_taxid': [5, 5, None],
            'ncbi_taxid': [10, 11, 2],
            'gambit_taxid': [10, 11, 2],
            'diameter': [0.9, 0.3, 0.1],
            'ngenomes': [3, 3, 0],
            'report': [1, 1, 0]
        }, index=pd.Index([1, 2, 5], name='species_taxid'))

        genome_metadata = pd.DataFrame({
            'uuid': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3'],
            'species_taxid': [1, 1, 1, 2, 2, 2],
            'species': ['All apart', 'All apart', 'All apart', 'Good species', 'Good species', 'Good species']
        }, index=pd.Index(['GCA_1', 'GCA_2', 'GCA_3', 'GCA_11', 'GCA_12', 'GCA_13'], name='assembly_accession'))

        # All 3 genomes of "All apart" are >0.7 from each other
        pw_data = {
            'GCA_1':  [0.0, 0.9, 0.8, 0.9, 0.9, 0.9],
            'GCA_2':  [0.9, 0.0, 0.9, 0.9, 0.9, 0.9],
            'GCA_3':  [0.8, 0.9, 0.0, 0.9, 0.9, 0.9],
            'GCA_11': [0.9, 0.9, 0.9, 0.0, 0.1, 0.1],
            'GCA_12': [0.9, 0.9, 0.9, 0.1, 0.0, 0.1],
            'GCA_13': [0.9, 0.9, 0.9, 0.1, 0.1, 0.0],
        }
        idx = ['GCA_1', 'GCA_2', 'GCA_3', 'GCA_11', 'GCA_12', 'GCA_13']
        pairwise = pd.DataFrame(pw_data, index=idx)

        # Write temp CSV for pairwise distances
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            pairwise.to_csv(f)
            pw_path = f.name

        try:
            ss = SplitSpecies(species, genome_metadata, pw_path, [],
                              0.7, 1, 'average', False)
            s, g, accessions_removed = ss.split_high_diameter_species()
            # All 3 genomes of "All apart" should be removed
            self.assertIn('GCA_1', accessions_removed)
            self.assertIn('GCA_2', accessions_removed)
            self.assertIn('GCA_3', accessions_removed)
        finally:
            os.unlink(pw_path)
