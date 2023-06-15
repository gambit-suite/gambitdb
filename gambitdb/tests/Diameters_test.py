# tests for the Diameters class
import unittest
import os
import filecmp
from gambitdb.Diameters  import Diameters

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','diameters')

class TestDiameters(unittest.TestCase):
    """
    Test class for the Diameters class.
    """
    
    # Theres a standard method to do this but cant remember the name
    # Compare the contents of two files, if different, print the diff
    def comp_files(self, file1, file2):
        """
    Compares the contents of two files, if different, prints the diff.
    Args:
      file1 (str): Path to the first file.
      file2 (str): Path to the second file.
    Returns:
      bool: True if the files are the same, False if they are different.
    Side Effects:
      Prints the diff of the two files if they are different.
    Examples:
      >>> TestDiameters.comp_files('file1.txt', 'file2.txt')
      The files are different. The following lines are different:
      <diff lines>
    """
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
        """
    Tests the Diameters class's calculate_diameters method.
    Args:
      None
    Returns:
      None
    Side Effects:
      Creates output files in the data/diameters directory.
    Examples:
      >>> TestDiameters.test_generating_diameters()
      None
    """
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


    def test_calc_thresholds_no_genomes(self):
        """
    Tests the Diameters class's calculate_thresholds method when there are no genomes.
    Args:
      None
    Returns:
      None
    Examples:
      >>> TestDiameters.test_calc_thresholds_no_genomes()
      None
    """
        self.cleanup()
        d = Diameters(os.path.join(data_dir, 'assembly_metadata.csv'), 
                      os.path.join(data_dir, 'pairwise_dist.csv'), 
                      os.path.join(data_dir, 'species_no_genomes.csv'), 
                      os.path.join(data_dir, 'output_species.csv'), 
                      os.path.join(data_dir, 'output_min_inter.csv'), 
                      False)
        
        genome_metadata, species, pairwise_distances = d.read_files()
        genomes_grouped_by_species_name = genome_metadata.groupby('species')
        number_of_species = species.shape[0]
        species_genomes = {}
        for species_name in species['name'].tolist():
            if species_name in genomes_grouped_by_species_name.groups:                
                species_genomes[species_name] = genomes_grouped_by_species_name.get_group(species_name)['assembly_accession'].values
            else:
                species_genomes[species_name] = []

        diameters, min_inter, ngenomes = d.calculate_thresholds(number_of_species, species_genomes, pairwise_distances)
        self.assertEqual(diameters.shape[0], 4)
        
    def cleanup(self):
        """
    Removes output files from the data/diameters directory.
    Args:
      None
    Returns:
      None
    Side Effects:
      Removes output files from the data/diameters directory.
    Examples:
      >>> TestDiameters.cleanup()
      None
    """
        for f in ['output_species.csv', 'output_min_inter.csv']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))
