# A class which will take in the pairwise distance file, accessions, and species taxon data and calculate diameters and minimums
import logging
import numpy
import pandas

class Diameters:
    def __init__(self, genome_assembly_metadata, pairwise_distances_filename, species_taxon_filename, species_taxon_output_filename, min_inter_output_filename, verbose):
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

    def calculate_diameters(self):
        self.logger.debug('calculate_diameters')

        # Read in the genome assembly filenames with path
        genome_metadata = pandas.read_csv(self.genome_assembly_metadata)
        # sort the genome_metadata dataframe by assembly_accession
        genome_metadata = genome_metadata.sort_values(by='assembly_accession')
        self.logger.debug('calculate_diameters: genome_metadata size: %s' % genome_metadata.shape[0])
        genomes_grouped_by_species_taxid = genome_metadata.groupby('species_taxid')
        self.logger.debug('calculate_diameters: genomes_grouped_by_species_taxid size: %s' % genomes_grouped_by_species_taxid.size())

        # Read in the species taxon file
        species = pandas.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')
        number_of_species = species.shape[0]
        self.logger.debug('calculate_diameters: species size: %s' % species.shape[0])

        # Read in the pairwise distances file
        pairwise_distances = pandas.read_csv(self.pairwise_distances_filename, index_col=0)
        # sort the pairwise_distances dataframe by the index column
        pairwise_distances = pairwise_distances.sort_index()
        self.logger.debug('calculate_diameters: pairwise_distances size: %s' % pairwise_distances.shape[0])

        # a dictionary of species_taxid: [assembly_accessions]
        species_genomes = {}
        # Loop over each species
        for species_taxid in species.index:
            # get a list of the assembly_accessions for the species and add to the species_genomes dictionary
            species_genomes[species_taxid] = genomes_grouped_by_species_taxid.get_group(species_taxid)['assembly_accession'].values

        diameters, min_inter = self.calculate_thresholds(number_of_species, species_genomes, pairwise_distances)

        # Extend the species taxon table and add in the diameters and number of genomes
        species['diameter'] = diameters
        species['ngenomes'] = genome_metadata.groupby('species_taxid').size()

        if number_of_species >0:
        # Take the min-inter values, add to a dataframe and write out to a new file
            mininter_df = pandas.DataFrame(min_inter, index=species.index, columns=species.index)
            mininter_df.to_csv(self.min_inter_output_filename)

        # The parent TaxonIDs (the genus) need to exist for the species
        species = self.create_mock_genus_rows(species)
        # Write out the species taxon table to a new file
        species.to_csv(self.species_taxon_output_filename)

    def create_mock_genus_rows(self, species):
        self.logger.debug('add_parent_ids')

        # Get all unique parent_taxids
        parent_taxids = species['parent_taxid'].unique()
        # get all species names and parent_taxid.  Split the species name to return the first word (genus name)
        species_names = species[['name', 'parent_taxid']].assign(genus_name=species['name'].str.split(' ').str[0])

        # create a dictionary of parent_taxid: genus_name
        parent_taxid_genus_name = dict(zip(species_names.parent_taxid, species_names.genus_name))

        genus_list = []
        # loop over the parent_taxids and add them to the parent_ids dataframe
        for parent_taxid in parent_taxids:
            genus_name = parent_taxid_genus_name[parent_taxid]
            genus_list.append([parent_taxid, genus_name, 'genus', '', parent_taxid, parent_taxid, 0, 0,0])
            
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
        self.logger.debug('calculate_thresholds')   
        # initalise the diameters and minimums, these get returned at the end
        diameters = numpy.zeros(number_of_species)
        min_inter = numpy.zeros((number_of_species, number_of_species))

        # loop over the species_genomes dictionary, looking up the index of the genome in the pairwise_distances dataframe
        for i, (species_taxid, assembly_accessions) in enumerate(species_genomes.items()):
            inds1 = pairwise_distances.index.get_indexer(assembly_accessions)
            # Find the maximum diameters for each species. Basically look at the pairwise distances
            # and find the maximum distance between any two genomes in the species
            diameters[i] = pairwise_distances.values[numpy.ix_(inds1, inds1)].max()
    
            for j, (species_taxid2, assembly_accessions2) in enumerate(species_genomes.items()):
                inds2 = pairwise_distances.index.get_indexer(assembly_accessions2)
                mi = pairwise_distances.values[numpy.ix_(inds1, inds2)].min()
                min_inter[i, j] = min_inter[j, i] = mi

        self.logger.debug('calculate_thresholds: diameters: %s' % diameters.shape[0])
        self.logger.debug('calculate_thresholds: min_inter: %s' % min_inter.shape[0])
        # return the diameters and minimums
        return diameters, min_inter
