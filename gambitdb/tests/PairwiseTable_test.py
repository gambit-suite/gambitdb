import unittest
import os
from gambitdb.PairwiseTable  import PairwiseTable

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','pairwise_table')

class TestPairwiseTable(unittest.TestCase):
    
    def test_generating_signatures(self):
        self.cleanup()
        pw = PairwiseTable(os.path.join(data_dir, 'assemblies'), os.path.join(data_dir, 'signatures.h5'), os.path.join(data_dir, 'pw-dists.csv'), 11, 'ATGAC',  None, 1, False)
        pw.generate_sigs_and_pairwise_table()

        self.assertTrue(os.path.exists(os.path.join(data_dir, 'signatures.h5')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'pw-dists.csv')))
        self.cleanup()

    def cleanup(self):
        for f in ['signatures.h5', 'pw-dists.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
