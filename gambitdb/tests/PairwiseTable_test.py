import unittest
import os
from gambitdb.PairwiseTable  import PairwiseTable

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','pairwise_table')

class TestPairwiseTable(unittest.TestCase):
    """
    Test class for PairwiseTable.
    """
    
    def test_generating_signatures(self):
        """
    Tests the generation of signatures and pairwise table.
    Args:
      self (TestPairwiseTable): The instance of the TestPairwiseTable class.
    Side Effects:
      Creates files in the data_dir directory.
    Notes:
      Asserts that the files were created.
    """
        self.cleanup()
        pw = PairwiseTable(os.path.join(data_dir, 'assemblies'), os.path.join(data_dir, 'signatures.h5'), os.path.join(data_dir, 'pw-dists.csv'), 11, 'ATGAC',  None, 1, False)
        pw.generate_sigs_and_pairwise_table()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'signatures.h5')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'pw-dists.csv')))
        self.cleanup()

    def cleanup(self):
        """
    Removes files from the data_dir directory.
    Args:
      self (TestPairwiseTable): The instance of the TestPairwiseTable class.
    Side Effects:
      Removes files from the data_dir directory.
    """
        for f in ['signatures.h5', 'pw-dists.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
