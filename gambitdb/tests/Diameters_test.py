# tests for the Diameters class
import unittest
import os
import filecmp
from gambitdb.Diameters  import Diameters

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','diameters')

class TestDiameters(unittest.TestCase):
    
    # Theres a standard method to do this but cant remember the name
    # Compare the contents of two files, if different, print the diff
    def comp_files(self, file1, file2):
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

        self.assertTrue(self.comp_files(os.path.join(data_dir, 'output_species.csv'), os.path.join(data_dir, 'expected_output_species.csv')))
        self.assertTrue(self.comp_files(os.path.join(data_dir, 'output_min_inter.csv'), os.path.join(data_dir, 'expected_output_min_inter.csv')))
        
        self.cleanup()

    def cleanup(self):
        for f in ['output_species.csv', 'output_min_inter.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
