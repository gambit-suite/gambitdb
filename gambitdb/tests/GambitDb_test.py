# tests for the GambitDb class
import unittest
import os
from gambitdb.GambitDb  import GambitDb

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','gambitdb')

class TestGambitDb(unittest.TestCase):
    """
    Test class for the GambitDb class.
    """
    
    def test_generate_pairwise_and_diameters(self):
        """
    Tests the generate_gambit_db method of the GambitDb class.
    Args:
      None
    Returns:
      None
    Side Effects:
      Creates files in the data_dir directory.
    Notes:
      Asserts that the files created by generate_gambit_db exist.
    Examples:
      >>> test_generate_pairwise_and_diameters()
    """
        self.cleanup()
        gdb = GambitDb(data_dir,
                        os.path.join(data_dir, 'assemblies'), 
                        os.path.join(data_dir, 'assembly_metadata.csv'),
                        os.path.join(data_dir, 'species.csv'),
                        None,
                        None,
                        os.path.join(data_dir, 'output_accessions_removed.csv'),
                        os.path.join(data_dir, 'output_species_removed.csv'),
                        os.path.join(data_dir, 'signatures.h5'), 
                        os.path.join(data_dir, 'pw-dists.csv'), 
                        os.path.join(data_dir, 'output_species_taxon.csv'),
                        os.path.join(data_dir, 'output_assembly_metadata.csv'),
                        11, 
                        'ATGAC', 
                        1,
                        1,
                        1,
                        0.99,
                        0.7,
                        1,
                        1,
                        None,
                        False)
        gdb.generate_gambit_db()
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'signatures.h5')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'pw-dists.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_species_taxon.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_assembly_metadata.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_species_removed.csv')))
        
        self.cleanup()

    def cleanup(self):
        """
    Deletes files created by the test_generate_pairwise_and_diameters method.
    Args:
      None
    Returns:
      None
    Side Effects:
      Deletes files in the data_dir directory.
    Notes:
      None
    Examples:
      >>> cleanup()
    """
        for f in ['signatures.h5', 'pw-dists.csv', 'output_species_taxon.csv', 'output_assembly_metadata.csv', 'output_accessions_removed.csv', 'output_species_removed.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))

