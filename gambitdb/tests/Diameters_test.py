# tests for the Diameters class
import unittest
import os
from gambitdb.Diameters  import Diameters

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','diameters')

class TestDiameters(unittest.TestCase):
    
    def test_generating_diameters(self):
        self.cleanup()
        d = Diameters(os.path.join(data_dir, 'assembly_metadata.csv'), 
                      os.path.join(data_dir, 'pairwise_dist.csv'), 
                      os.path.join(data_dir, 'species.csv'), 
                      os.path.join(data_dir, 'output_species.csv'), 
                      os.path.join(data_dir, 'output_min_inter.csv'), 
                      False)
        d.calculate_diameters()
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_species.csv')))
        self.assertTrue(os.path.exists(os.path.join(data_dir, 'output_min_inter.csv')))
        
        self.cleanup()

    def cleanup(self):
        for f in ['output_species.csv', 'output_min_inter.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))