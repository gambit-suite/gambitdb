# A class which will take in the pairwise distance file, accessions, and species taxon data and calculate diameters and minimums
import logging
import numpy
import pandas

class Diameters:
    """
    A class which will take in the pairwise distance file, accessions, and species taxon data and calculate diameters and minimums.
    """
    def __init__(self, genome_assembly_metadata, pairwise_distances_filename, species_taxon_filename, species_taxon_output_filename, min_inter_output_filename, verbose):
        """
    Initializes the Diameters class.
    Args:
      genome_assembly_metadata (str): The path to the genome assembly metadata file.
      pairwise_distances_filename (str): The path to the pairwise distances file.
      species_taxon_filename (str): The path to the species taxon file.
      species_taxon_output_filename (str): The path to the species taxon output file.
      min_inter_output_filename (str): The path to the min inter output file.
      verbose (bool): Whether to print debug messages.
    Side Effects:
      Sets the logger level to DEBUG if verbose is True, otherwise sets it to ERROR.
    """
        self.logger = logging.getLogger(__name__)

        # input files
        self.genome_assembly_metadata = genome_assembly_metadata
        self.pairwise_distances_filename = pairwise_distances_filename
        self.species_taxon_filename = species_taxon_filename

        # output files
        self.species_taxon_output_filename = species_taxon_output_filename
        self.min_inter_output_filename = min_inter_output_filename

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def read_files(self):
        """
    Reads in the genome assembly filenames with path, species taxon file, and pairwise distances file.
    Returns:
      tuple: A tuple containing the genome metadata, species, and pairwise distances dataframes.
    """
        # Read in the genome assembly filenames with path
        genome_metadata = pandas.read_csv(self.genome_assembly_metadata)
        # sort the genome_metadata dataframe by assembly_accession
        genome_metadata = genome_metadata.sort_values(by='assembly_accession')
        self.logger.debug('calculate_diameters: genome_metadata size: %s' % genome_metadata.shape[0])
        
        # Read in the species taxon file
        species = pandas.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')

        # Read in the pairwise distances file
        pairwise_distances = pandas.read_csv(self.pairwise_distances_filename, index_col=0)
        # sort the pairwise_distances dataframe by the index column
        pairwise_distances = pairwise_distances.sort_index()
        self.logger.debug('calculate_diameters: pairwise_distances size: %s' % pairwise_distances.shape[0])

        return genome_metadata, species, pairwise_distances

    def calculate_diameters(self):
        """
    Calculates the diameters and minimums.
    Side Effects:
      Writes out the species taxon table and min inter output files.
    """
        self.logger.debug('calculate_diameters')
        genome_metadata, species, pairwise_distances = self.read_files()

        genomes_grouped_by_species_name = genome_metadata.groupby('species')
        self.logger.debug('calculate_diameters: genomes_grouped_by_species_name size: %s' % genomes_grouped_by_species_name.size())

        number_of_species = species.shape[0]
        self.logger.debug('calculate_diameters: species size: %s' % species.shape[0])

        # a dictionary of species_taxid: [assembly_accessions]
        species_genomes = {}
        # Loop over each species
        for species_name in species['name'].tolist():
            # get a list of the assembly_accessions for the species and add to the species_genomes dictionary
            #check to see if species_name is in genomes_grouped_by_species_name
            if species_name in genomes_grouped_by_species_name.groups:                
                species_genomes[species_name] = genomes_grouped_by_species_name.get_group(species_name)['assembly_accession'].values
            else:
                species_genomes[species_name] = []

        diameters, min_inter, ngenomes = self.calculate_thresholds(number_of_species, species_genomes, pairwise_distances)
        # Extend the species taxon table and add in the diameters and number of genomes
        species['diameter'] = diameters
        species['ngenomes'] = ngenomes
        species['ngenomes'] = species['ngenomes'].astype(int)

        if number_of_species >0:
        # Take the min-inter values, add to a dataframe and write out to a new file
            mininter_df = pandas.DataFrame(min_inter, index=species.index, columns=species.index)
            mininter_df.to_csv(self.min_inter_output_filename)

        # The parent TaxonIDs (the genus) need to exist for the species
        species = self.create_mock_genus_rows(species)
        # Write out the species taxon table to a new file
        species.to_csv(self.species_taxon_output_filename)

    def create_mock_genus_rows(self, species):
        """
    Creates mock genus rows for the species taxon table.
    Args:
      species (pandas.DataFrame): The species taxon dataframe.
    Returns:
      pandas.DataFrame: The species taxon dataframe with mock genus rows.
    """
        self.logger.debug('add_parent_ids')

        # Get all unique parent_taxids
        parent_taxids = species['parent_taxid'].unique()
        species_taxids = species.index.unique().tolist()

        # get all species names and parent_taxid.  Split the species name to return the first word (genus name)
        species_names = species[['name', 'parent_taxid']].assign(genus_name=species['name'].str.split(' ').str[0])

        # create a dictionary of parent_taxid: genus_name
        parent_taxid_genus_name = dict(zip(species_names.parent_taxid, species_names.genus_name))

        # Group the species by parent_taxid and get the maximum diameter for each parent_taxid
        parent_taxid_diameters = species.groupby('parent_taxid')['diameter'].max()
        # create a dictionary where the key is the parent_taxid and the value is the diameter
        parent_taxid_diameters_dict = parent_taxid_diameters.to_dict()
        
        genus_list = []
        # loop over the parent_taxids and add them to the parent_ids dataframe
        for parent_taxid in parent_taxids:
            # A subspecies will have the species as the parent taxid, so skip these
            if parent_taxid in species_taxids:
                continue
            genus_name = parent_taxid_genus_name[parent_taxid]
            genus_list.append([parent_taxid, genus_name, 'genus', '', parent_taxid, parent_taxid, parent_taxid_diameters_dict[parent_taxid], 0,1])
            
        df_extended = pandas.DataFrame(genus_list, columns=['species_taxid', 
                                                            'name', 
                                                            'rank', 
                                                            'parent_taxid', 
                                                            'ncbi_taxid', 
                                                            'gambit_taxid', 
                                                            'diameter', 
                                                            'ngenomes', 'report']).set_index('species_taxid')
        species = pandas.concat([species, df_extended])

        return species

    # Calculate the diameters and minimums
    def calculate_thresholds(self, number_of_species, species_genomes, pairwise_distances):
        """
    Calculates the diameters and minimums for a given set of species.
    Args:
      number_of_species (int): The number of species in the dataset.
      species_genomes (dict): A dictionary of species taxon and assembly accessions.
      pairwise_distances (pandas.DataFrame): A dataframe of pairwise distances.
    Returns:
      tuple: A tuple containing the diameters, minimums, and number of genomes for each species.
    Examples:
      >>> diameters, min_inter, ngenomes = Diameters.calculate_thresholds(number_of_species, species_genomes, pairwise_distances)
      >>> diameters
      array([0.1, 0.2, 0.3])
      >>> min_inter
      array([[0.1, 0.2, 0.3],
             [0.2, 0.1, 0.4],
             [0.3, 0.4, 0.1]])
      >>> ngenomes
      array([2, 3, 4])
    """
        self.logger.debug('calculate_thresholds')   
        # initalise the diameters and minimums, these get returned at the end
        diameters = numpy.zeros(number_of_species)
        ngenomes = numpy.zeros(number_of_species)
        min_inter = numpy.zeros((number_of_species, number_of_species))

        # loop over the species_genomes dictionary, looking up the index of the genome in the pairwise_distances dataframe
        for i, (species_taxid, assembly_accessions) in enumerate(species_genomes.items()):
            # if the accessbly accessions list is empty then continue
            if len(assembly_accessions) == 0:
                continue

            inds1 = pairwise_distances.index.get_indexer(assembly_accessions)
            # Find the maximum diameters for each species. Basically look at the pairwise distances
            # and find the maximum distance between any two genomes in the species
            diameters[i] = pairwise_distances.values[numpy.ix_(inds1, inds1)].max()
            ngenomes[i] = len(assembly_accessions)
            for j, (species_taxid2, assembly_accessions2) in enumerate(species_genomes.items()):
                if len(assembly_accessions2) == 0:
                    continue
                inds2 = pairwise_distances.index.get_indexer(assembly_accessions2)
                mi = pairwise_distances.values[numpy.ix_(inds1, inds2)].min()
                min_inter[i, j] = min_inter[j, i] = mi

        self.logger.debug('calculate_thresholds: diameters: %s' % diameters.shape[0])
        self.logger.debug('calculate_thresholds: min_inter: %s' % min_inter.shape[0])
        # return the diameters and minimums
        return diameters, min_inter, ngenomes
