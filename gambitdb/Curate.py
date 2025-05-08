# this is a method which will allow you to curate the input genomes and species taxon files
# it will output new files which can be used for the core database
# it will also output a list of accessions removed and a list of species removed
# it will also output a list of species which have been removed because they have too few genomes
# it will also output a list of species which have been removed because they are in the species_to_remove file
# it will also output a list of accessions which have been removed because they are in the accessions_to_remove file
# it will also output a list of species which have been removed because they are in the accessions_to_remove file

import pandas
import logging
import os
import gzip
from gambitdb.SplitSpecies import SplitSpecies


class Curate:
    """
    A class for curating input genomes and species taxon files.
    """
    def __init__(self, species_taxon_filename, genome_assembly_metadata, assembly_directory, pairwise_distances_filename, 
                 species_to_remove, accessions_to_remove, species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, accession_removed_output_filename, 
                 species_removed_output_filename, minimum_ngenomes,small_cluster_ngenomes,
                 small_cluster_diameter, maximum_diameter, minimum_cluster_size, linkage_method,  debug, verbose):
        """
    Initializes the Curate class.
    Args:
      species_taxon_filename (str): The path to the species taxon file.
      genome_assembly_metadata (str): The path to the genome assembly metadata file.
      assembly_directory (str): The path to the assembly directory.
      pairwise_distances_filename (str): The path to the pairwise distances file.
      species_to_remove (str): The path to the species to remove file.
      accessions_to_remove (str): The path to the accessions to remove file.
      species_taxon_output_filename (str): The path to the species taxon output file.
      genome_assembly_metadata_output_filename (str): The path to the genome assembly metadata output file.
      accession_removed_output_filename (str): The path to the accession removed output file.
      species_removed_output_filename (str): The path to the species removed output file.
      minimum_ngenomes (int): The minimum number of genomes for a species to be included.
      small_cluster_ngenomes (int): The minimum number of genomes for a small cluster to be included.
      small_cluster_diameter (float): The maximum diameter for a small cluster to be included.
      maximum_diameter (float): The maximum diameter for a cluster to be included.
      minimum_cluster_size (int): The minimum size for a cluster to be included.
      linkage_method (str): The linkage method for the clustering.
      debug (bool): Whether to enable debug logging.
      verbose (bool): Whether to enable verbose logging.
    Attributes:
      species_removed (list): A list of species removed.
      species_taxonids_removed (list): A list of species taxonids removed.
      accessions_removed (list): A list of accessions removed.
      logger (logging.Logger): The logger for the class.
    Examples:
      >>> curate = Curate(species_taxon_filename, genome_assembly_metadata, assembly_directory, pairwise_distances_filename, 
                          species_to_remove, accessions_to_remove, species_taxon_output_filename, 
                          genome_assembly_metadata_output_filename, accession_removed_output_filename, 
                          species_removed_output_filename, minimum_ngenomes,small_cluster_ngenomes,
                          small_cluster_diameter, maximum_diameter, minimum_cluster_size,  debug, verbose)
    """
        self.logger = logging.getLogger(__name__)
        self.species_taxon_filename = species_taxon_filename
        self.genome_assembly_metadata = genome_assembly_metadata
        self.assembly_directory = assembly_directory
        self.pairwise_distances_filename = pairwise_distances_filename
        self.species_to_remove = species_to_remove
        self.accessions_to_remove = accessions_to_remove
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accession_removed_output_filename = accession_removed_output_filename
        self.species_removed_output_filename = species_removed_output_filename
        self.minimum_ngenomes = minimum_ngenomes
        self.small_cluster_ngenomes = small_cluster_ngenomes
        self.small_cluster_diameter = small_cluster_diameter
        self.maximum_diameter = maximum_diameter
        self.minimum_cluster_size = minimum_cluster_size
        self.linkage_method = linkage_method
        
        self.debug = debug
        self.verbose = verbose

        self.species_removed = []
        self.species_taxonids_removed = []
        self.accessions_removed = []
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def remove_species_using_input_file(self, species):
        """
    Removes species from the species dataframe using a file of species names.
    Args:
      species (pandas.DataFrame): The species dataframe.
    Returns:
      pandas.DataFrame: The species dataframe with the species removed.
    Side Effects:
      Adds the removed species to the list self.species_removed.
      Adds the removed species taxonids to the list self.species_taxonids_removed.
    Notes:
      The file of species names must be in the same directory as the code.
    Examples:
      >>> species = pandas.DataFrame({'name': ['species1', 'species2', 'species3'], 'ngenomes': [1, 2, 3], 'diameter': [1, 2, 3]}, index=[1, 2, 3])
      >>> species_to_remove = pandas.DataFrame({'name': ['species2']}, index=[2])
      >>> curate.remove_species_using_input_file(species)
      pandas.DataFrame({'name': ['species1', 'species3'], 'ngenomes': [1, 3], 'diameter': [1, 3]}, index=[1, 3])
    """
        self.logger.debug('remove_species')
        if self.species_to_remove and os.path.isfile(self.species_to_remove):
            # species to remove contains the names of species. need to translate to species_taxonids
            species_names_to_remove = pandas.read_csv(self.species_to_remove, header=None)
            species_names_to_remove = species_names_to_remove[0].tolist()

            species_taxonids_to_remove = species[species.isin(species_names_to_remove)['name']].index.tolist()
            species = species[~species.isin(species_names_to_remove)['name']]

            # append species_to_remove to list self.species_removed
            self.species_taxonids_removed = self.species_taxonids_removed + species_taxonids_to_remove
            self.species_removed = self.species_removed + species_names_to_remove 
        
        self.logger.debug('remove_species: species_removed: %s' % len(self.species_removed))
        self.logger.debug('remove_species: species_taxonids_removed: %s' % len(self.species_taxonids_removed))
        self.logger.debug('remove_species: species: %s' % species.shape[0])
        return species

    def remove_accessions_using_input_file(self, genome_metadata):
        """
    Removes accessions from the genome metadata dataframe using a file of accessions.
    Args:
      genome_metadata (pandas.DataFrame): The genome metadata dataframe.
    Returns:
      pandas.DataFrame: The genome metadata dataframe with the accessions removed.
    Side Effects:
      Adds the removed accessions to the list self.accessions_removed.
    Notes:
      The file of accessions must be in the same directory as the code.
    Examples:
      >>> genome_metadata = pandas.DataFrame({'species_taxid': [1, 2, 3], 'assembly_accession': ['accession1', 'accession2', 'accession3']}, index=['accession1', 'accession2', 'accession3'])
      >>> accessions_to_remove = pandas.DataFrame({'accession': ['accession2']}, index=['accession2'])
      >>> curate.remove_accessions_using_input_file(genome_metadata)
      pandas.DataFrame({'species_taxid': [1, 3], 'assembly_accession': ['accession1', 'accession3']}, index=['accession1', 'accession3'])
    """
        self.logger.debug('remove_accessions')
        if self.accessions_to_remove and os.path.isfile(self.accessions_to_remove) and os.stat(self.accessions_to_remove).st_size != 0:
            accessions_to_remove = pandas.read_csv(self.accessions_to_remove, header=None)
            accessions_to_remove = accessions_to_remove[0].tolist()
            genome_metadata = genome_metadata[~genome_metadata.index.isin(accessions_to_remove)]
            self.logger.debug('remove_accessions: accessions_to_remove using input file: %s' % len(accessions_to_remove))

            # append accessions_to_remove to list self.accessions_removed
            self.accessions_removed = self.accessions_removed + accessions_to_remove
        return genome_metadata

    # given the genome_metadata dataframe, remove assembly_accession which are not contained in the filenames of the assembly_directory
    def remove_accessions_for_truncated_assemblies(self, genome_metadata):
        self.logger.debug('remove_accessions_for_truncated_assemblies')
        if self.assembly_directory is None:
            return genome_metadata

        # get the list of filenames in the assembly_directory
        assembly_filenames = os.listdir(self.assembly_directory)
        # get the list of assembly_accession in the genome_metadata dataframe
        assembly_accessions = genome_metadata.index.tolist()

        accessions_with_no_matching_files = []
        for accession in assembly_accessions:
            accession_found = False
            for assembly_file in assembly_filenames:
                if accession in assembly_file:
                    accession_found = True
                    # open the accession file which is gzipped. Open the file to test if it is truncated
                    with gzip.open(self.assembly_directory + '/' + assembly_file, 'rb') as f:
                        file_content = f.read()
                        if file_content[-1:] != b'\x00':
                            accession_found = False
                            accessions_with_no_matching_files.append(accession)
                            self.logger.debug('remove_accessions_for_truncated_assemblies: accession: %s' % accession)
                    break
            if accession_found == False:
                accessions_with_no_matching_files.append(accession)
                self.logger.debug('remove_accessions_not_in_assembly_directory: accession: %s' % accession)

        # remove the assembly_accession_to_remove from the genome_metadata dataframe
        genome_metadata = genome_metadata[~genome_metadata.index.isin(accessions_with_no_matching_files)]

        # append assembly_accession_to_remove to list self.accessions_removed
        self.accessions_removed = self.accessions_removed + accessions_with_no_matching_files
        self.logger.debug('remove_accessions_not_in_assembly_directory: accessions_removed: %s' % len(self.accessions_removed))
        return genome_metadata

    # given the genome_metadata dataframe, remove assembly_accession which are not contained in the filenames of the assembly_directory
    def remove_accessions_not_in_assembly_directory(self, genome_metadata):
        """
    Removes accessions from the genome metadata dataframe if they are not contained in the filenames of the assembly directory.
    Args:
      genome_metadata (pandas.DataFrame): The genome metadata dataframe.
    Returns:
      pandas.DataFrame: The genome metadata dataframe with the accessions removed.
    Side Effects:
      Adds the removed accessions to the list self.accessions_removed.
    Notes:
      The assembly directory must be in the same directory as the code.
    Examples:
      >>> genome_metadata = pandas.DataFrame({'species_taxid': [1, 2, 3], 'assembly_accession': ['accession1', 'accession2', 'accession3']}, index=['accession1', 'accession2', 'accession3'])
      >>> assembly_filenames = ['accession1.fa', 'accession2.fa', 'accession4.fa']
      >>> curate.remove_accessions_not_in_assembly_directory(genome_metadata)
      pandas.DataFrame({'species_taxid': [1, 2], 'assembly_accession': ['accession1', 'accession2']}, index=['accession1', 'accession2'])
    """
        self.logger.debug('remove_accessions_not_in_assembly_directory')
        if self.assembly_directory is None:
            return genome_metadata

        # get the list of filenames in the assembly_directory
        assembly_filenames = os.listdir(self.assembly_directory)
        # get the list of assembly_accession in the genome_metadata dataframe
        assembly_accessions = genome_metadata.index.tolist()

        accessions_with_no_matching_files = []
        for accession in assembly_accessions:
            accession_found = False
            for assembly_file in assembly_filenames:
                if accession in assembly_file:
                    accession_found = True
                    break
            if accession_found == False:
                accessions_with_no_matching_files.append(accession)
                self.logger.debug('remove_accessions_not_in_assembly_directory: accession: %s' % accession)

        # remove the assembly_accession_to_remove from the genome_metadata dataframe
        genome_metadata = genome_metadata[~genome_metadata.index.isin(accessions_with_no_matching_files)]

        # append assembly_accession_to_remove to list self.accessions_removed
        self.accessions_removed = self.accessions_removed + accessions_with_no_matching_files
        self.logger.debug('remove_accessions_not_in_assembly_directory: accessions_removed: %s' % len(self.accessions_removed))
        return genome_metadata

    # remove species where ngenomes < an integer number provided
    def remove_species_with_fewer_than_n_genomes(self, species):
        """
    Removes species from the species dataframe if they have fewer than n genomes.
    Args:
      species (pandas.DataFrame): The species dataframe.
    Returns:
      pandas.DataFrame: The species dataframe with the species removed.
    Side Effects:
      Adds the removed species to the list self.species_removed.
      Adds the removed species taxonids to the list self.species_taxonids_removed.
    Notes:
      The number of genomes is determined by the parameter self.minimum_ngenomes.
    Examples:
      >>> species = pandas.DataFrame({'name': ['species1', 'species2', 'species3'], 'ngenomes': [1, 2, 3], 'diameter': [1, 2, 3]}, index=[1, 2, 3])
      >>> curate.minimum_ngenomes = 2
      >>> curate.remove_species_with_fewer_than_n_genomes(species)
      pandas.DataFrame({'name': ['species2', 'species3'], 'ngenomes': [2, 3], 'diameter': [2, 3]}, index=[2, 3])
    """
        self.logger.debug('remove_species_with_fewer_than_n_genomes')
        # create a list of species names where the species have fewer than n genomes
        species_taxonids_to_remove = species[species['ngenomes'] < self.minimum_ngenomes].index.tolist()
        # get the species column called names and add them to a list species_to_remove
        species_to_remove = species[species['ngenomes'] < self.minimum_ngenomes]['name'].tolist()
        species = species[species['ngenomes'] >= self.minimum_ngenomes]

        self.species_removed = self.species_removed + species_to_remove
        self.species_taxonids_removed = self.species_taxonids_removed + species_taxonids_to_remove

        self.logger.debug('remove_species_with_fewer_than_n_genomes: species_removed: %s' % len(self.species_removed))
        self.logger.debug('remove_species_with_fewer_than_n_genomes: species_taxonids_removed: %s' % len(self.species_taxonids_removed))
        self.logger.debug('remove_species_with_fewer_than_n_genomes: species: %s' % species.shape[0])
        return species

    def remove_species_with_zero_diameter(self, species):
        """
    Removes species from the species dataframe if they have a diameter of 0.
    Args:
      species (pandas.DataFrame): The species dataframe.
    Returns:
      pandas
    """
        self.logger.debug('remove_species_with_zero_diameter')
        # create a list of species names where the species have a diameter of 0
        species_taxonids_to_remove = species[species['diameter'] == 0].index.tolist()
        species_to_remove = species[species['diameter'] == 0]['name'].tolist()

        species = species[species['diameter'] != 0]
        self.species_removed = self.species_removed + species_to_remove
        self.species_taxonids_removed = self.species_taxonids_removed + species_taxonids_to_remove

        self.logger.debug('remove_species_with_zero_diameter: species_removed: %s' % len(self.species_removed))
        self.logger.debug('remove_species_with_zero_diameter: species_taxonids_removed: %s' % len(self.species_taxonids_removed))
        self.logger.debug('remove_species_with_zero_diameter: species: %s' %  species.shape[0])
        return species

    def remove_species_where_cluster_is_small_and_diameter_is_large(self, species):
        self.logger.debug('remove_species_where_cluster_is_small_and_diameter_is_large')
        species_taxonids_to_remove = species[(species['ngenomes'] < self.small_cluster_ngenomes) & (species['diameter'] > self.small_cluster_diameter)].index.tolist()
        species_to_remove = species[(species['ngenomes'] < self.small_cluster_ngenomes) & (species['diameter'] > self.small_cluster_diameter)]['name'].tolist()

        species = species[~((species['ngenomes'] < self.small_cluster_ngenomes) & (species['diameter'] > self.small_cluster_diameter))]
        self.species_removed = self.species_removed + species_to_remove
        self.species_taxonids_removed = self.species_taxonids_removed + species_taxonids_to_remove

        self.logger.debug('remove_species_where_cluster_is_small_and_diameter_is_large: species_removed: %s' % len(self.species_removed))
        self.logger.debug('remove_species_where_cluster_is_small_and_diameter_is_large: species_taxonids_removed: %s' % len(self.species_taxonids_removed))
        self.logger.debug('remove_species_where_cluster_is_small_and_diameter_is_large: species: %s' % species.shape[0])
        return species
    
    def remove_genomes_where_the_species_has_been_removed(self, genome_metadata):
        self.logger.debug('remove_genomes_where_the_species_has_been_removed')
        # remove genomes where the species has been removed

        genomes_to_remove = genome_metadata['species_taxid'].isin(self.species_taxonids_removed)
        genome_metadata = genome_metadata[~genomes_to_remove]

        self.logger.debug('remove_genomes_where_the_species_has_been_removed: genomes_to_remove: %s' % genome_metadata['species_taxid'].isin(self.species_taxonids_removed).tolist().count(True))
        self.logger.debug('remove_genomes_where_the_species_has_been_removed: genome_metadata: %s' % genome_metadata.shape[0])
        return genome_metadata

    def filter_spreadsheets_and_output_new_files(self):
        """
    Reads in the species taxon file and genome assembly filenames with path, filters them, and outputs new files.
    Args:
      self (Curate): The Curate object.
    Returns:
      None
    Side Effects:
      Outputs new files.
    Notes:
      This method calls other methods in the Curate class.
    Examples:
      >>> curate = Curate()
      >>> curate.filter_spreadsheets_and_output_new_files()
    """
        self.logger.debug('filter_spreadsheets_and_output_new_files')
        # Read in the species taxon file
        species = pandas.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')

        # Filter species
        species = self.remove_species_using_input_file(species)
        species = self.remove_species_with_fewer_than_n_genomes(species)
        species = self.remove_species_with_zero_diameter(species)
        species = self.remove_species_where_cluster_is_small_and_diameter_is_large(species)

        # Read in the genome assembly filenames with path
        genome_metadata = pandas.read_csv(self.genome_assembly_metadata)
        genome_metadata = genome_metadata.set_index('assembly_accession')

        # Filter genome metadata
        genome_metadata = self.remove_accessions_using_input_file(genome_metadata)
        genome_metadata = self.remove_genomes_where_the_species_has_been_removed(genome_metadata)
        genome_metadata = self.remove_accessions_not_in_assembly_directory(genome_metadata)
        #genome_metadata = self.remove_accessions_for_truncated_assemblies(genome_metadata)
        
        # Split species into subspecies
        split_species_obj = SplitSpecies(species, 
                                         genome_metadata, 
                                         self.pairwise_distances_filename,
                                         self.accessions_removed, 
                                         self.maximum_diameter, self.minimum_cluster_size, self.linkage_method, self.verbose)
        species, genome_metadata, self.accessions_removed = split_species_obj.split_high_diameter_species()
        self.write_output_files(species, genome_metadata)

    def write_output_files(self, species, genome_metadata):
        """
    Writes out the species taxon table, genome assembly metadata table, accessions removed, and species removed to new files.
    Args:
      self (Curate): The Curate object.
      species (pandas.DataFrame): The species taxon table.
      genome_metadata (pandas.DataFrame): The genome assembly metadata table.
    Returns:
      None
    Side Effects:
      Outputs new files.
    Examples:
      >>> curate = Curate()
      >>> curate.write_output_files(species, genome_metadata)
    """
        self.logger.debug('write_output_files')
        # set the species_taxonids to integers in the species dataframe
        species.index = species.index.astype(int)
        species['parent_taxid'] = species['parent_taxid'].astype(int)
        species['ncbi_taxid'] = species['ncbi_taxid'].astype(int)
        species['gambit_taxid'] = species['gambit_taxid'].astype(int)

        # Write out the species taxon table to a new file
        species.to_csv(self.species_taxon_output_filename)

        # set the species_taxonids to integers in the genome_metadata dataframe
        genome_metadata['species_taxid'] = genome_metadata['species_taxid'].astype(int)
        # Write out the genome assembly metadata table to a new file
        genome_metadata.to_csv(self.genome_assembly_metadata_output_filename)
        # Write out the accessions removed to a new file
        self.write_accessions_removed_to_file()
        # Write out the species removed to a new file
        self.write_species_removed_to_file()

        # Write out the number of species removed
        self.logger.info('write_output_files: Number of species removed: ' + str(len(self.species_removed)))
        self.logger.info('write_output_files: Number of accessions removed: ' + str(len(self.accessions_removed)))
        self.logger.debug('write_output_files: No. genomes in metadata' + str(genome_metadata.shape[0]))

    def write_accessions_removed_to_file(self):
        """
    Writes out the accessions removed to a new file.
    Args:
      self (Curate): The Curate object.
    Returns:
      None
    Side Effects:
      Outputs a new file.
    Examples:
      >>> curate = Curate()
      >>> curate.write_accessions_removed_to_file()
    """
        self.logger.debug('write_accessions_to_file')
        self.write_list_to_file(self.accessions_removed, self.accession_removed_output_filename)
     
    def write_species_removed_to_file(self):
        """
    Writes out the species removed to a new file.
    Args:
      self (Curate): The Curate object.
    Returns:
      None
    Side Effects:
      Outputs a new file.
    Examples:
      >>> curate = Curate()
      >>> curate.write_species_removed_to_file()
    """
        self.logger.debug('write_species_removed_to_file')
        self.write_list_to_file(self.species_removed, self.species_removed_output_filename)

    def write_list_to_file(self, list_to_write, output_filename):
        """
    Writes out a list to a new file.
    Args:
      self (Curate): The Curate object.
      list_to_write (list): The list to write.
      output_filename (str): The output filename.
    Returns:
      None
    Side Effects:
      Outputs a new file.
    Examples:
      >>> curate = Curate()
      >>> curate.write_list_to_file(list_to_write, output_filename)
    """
        if len(list_to_write) > 0:
            with open(output_filename, 'w') as f:
                for item in list_to_write:
                    f.write(str(item) + '\n')
