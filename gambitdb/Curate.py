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

class Curate:
    def __init__(self, species_taxon_filename, genome_assembly_metadata, assembly_directory, species_to_remove, accessions_to_remove, species_taxon_output_filename, genome_assembly_metadata_output_filename, accession_removed_output_filename, species_removed_output_filename, minimum_ngenomes, debug, verbose):
        self.logger = logging.getLogger(__name__)
        self.species_taxon_filename = species_taxon_filename
        self.genome_assembly_metadata = genome_assembly_metadata
        self.assembly_directory = assembly_directory
        self.species_to_remove = species_to_remove
        self.accessions_to_remove = accessions_to_remove
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accession_removed_output_filename = accession_removed_output_filename
        self.species_removed_output_filename = species_removed_output_filename
        self.minimum_ngenomes = minimum_ngenomes
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
        self.logger.debug('remove_accessions')
        if self.accessions_to_remove and os.path.isfile(self.accessions_to_remove):
            accessions_to_remove = pandas.read_csv(self.accessions_to_remove, header=None)
            accessions_to_remove = accessions_to_remove[0].tolist()
            genome_metadata = genome_metadata[~genome_metadata.index.isin(accessions_to_remove)]
            self.logger.debug('remove_accessions: accessions_to_remove using input file: %s' % len(accessions_to_remove))

            # append accessions_to_remove to list self.accessions_removed
            self.accessions_removed = self.accessions_removed + accessions_to_remove
        return genome_metadata

    # given the genome_metadata dataframe, remove assembly_accession which are not contained in the filenames of the assembly_directory
    def remove_accessions_not_in_assembly_directory(self, genome_metadata):
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
    
    def remove_genomes_where_the_species_has_been_removed(self, genome_metadata):
        self.logger.debug('remove_genomes_where_the_species_has_been_removed')
        # remove genomes where the species has been removed

        genomes_to_remove = genome_metadata['species_taxid'].isin(self.species_taxonids_removed)
        genome_metadata = genome_metadata[~genomes_to_remove]

        self.logger.debug('remove_genomes_where_the_species_has_been_removed: genomes_to_remove: %s' % genome_metadata['species_taxid'].isin(self.species_taxonids_removed).tolist().count(True))
        self.logger.debug('remove_genomes_where_the_species_has_been_removed: genome_metadata: %s' % genome_metadata.shape[0])
        return genome_metadata

    def filter_spreadsheets_and_output_new_files(self):
        self.logger.debug('filter_spreadsheets_and_output_new_files')
        # Read in the species taxon file
        species = pandas.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')

        # Filter species
        species = self.remove_species_using_input_file(species)
        species = self.remove_species_with_fewer_than_n_genomes(species)
        species = self.remove_species_with_zero_diameter(species)

        # Read in the genome assembly filenames with path
        genome_metadata = pandas.read_csv(self.genome_assembly_metadata)
        genome_metadata = genome_metadata.set_index('assembly_accession')

        # Filter genome metadata
        genome_metadata = self.remove_accessions_using_input_file(genome_metadata)
        genome_metadata = self.remove_genomes_where_the_species_has_been_removed(genome_metadata)
        genome_metadata = self.remove_accessions_not_in_assembly_directory(genome_metadata)

        self.write_output_files(species, genome_metadata)

    def write_output_files(self, species, genome_metadata):
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
        self.logger.debug('write_accessions_to_file')
        self.write_list_to_file(self.accessions_removed, self.accession_removed_output_filename)
     
    def write_species_removed_to_file(self):
        self.logger.debug('write_species_removed_to_file')
        self.write_list_to_file(self.species_removed, self.species_removed_output_filename)

    def write_list_to_file(self, list_to_write, output_filename):
        if len(list_to_write) > 0:
            with open(output_filename, 'w') as f:
                for item in list_to_write:
                    f.write(str(item) + '\n')
