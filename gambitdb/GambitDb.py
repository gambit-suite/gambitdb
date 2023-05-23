# A class to generate a Gambit database

import os
import logging

from gambitdb.PairwiseTable import PairwiseTable
from gambitdb.Diameters import Diameters

class GambitDb:

    def __init__(self, output_directory, assembly_directory, genome_assembly_metadata, species_taxon_filename,  signatures_output_filename, database_output_filename, kmer, kmer_prefix, verbose):
        self.logger = logging.getLogger(__name__)
        self.output_directory = output_directory
        self.assembly_directory = assembly_directory
        self.genome_assembly_metadata = genome_assembly_metadata
        self.signatures_output_filename = signatures_output_filename
        self.database_output_filename = database_output_filename
        self.kmer = kmer
        self.kmer_prefix = kmer_prefix

        # This file may not exist at this point
        self.species_taxon_filename = species_taxon_filename 

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def check_species_taxonid_file_exists_or_create_one(self, species_taxon_filename, genome_assembly_metadata):
        # check if the species_taxon_filename exists, if it does return the filename
        # otherwise generate a file using the genome_assembly_metadata
        if species_taxon_filename is not None and os.path.isfile(species_taxon_filename):
            return species_taxon_filename
        else: 
            self.logger.debug('species_taxon_file: %s does not exist', species_taxon_filename)
            # generate the species_taxon_filename
            #self.generate_species_taxon_file(species_taxon_filename)
            self.logger.error('not implemented yet. this is to allow for GTDB in the future')
            os.system('exit')   
            return 1

    def generate_gambit_db(self):
        self.logger.debug('generate_gambit_db')
        pairwise = self.generate_pairwise_table()

        self.species_taxon_filename = self.check_species_taxonid_file_exists_or_create_one(self.species_taxon_filename, self.genome_assembly_metadata)
        diameters = self.generate_diameters(pairwise.distance_table_output_filename, self.species_taxon_filename)

    def generate_pairwise_table(self):
        self.logger.debug('generate_pairwise_table')
    
        pw = PairwiseTable(self.assembly_directory, 
                           os.path.join(self.output_directory, self.signatures_output_filename), 
                           os.path.join(self.output_directory, 'pw-dists.csv'), 
                           self.kmer, 
                           self.kmer_prefix, 
                           self.verbose)
        pw.generate_sigs_and_pairwise_table()
        return pw

    def generate_diameters(self, distance_table, species_taxon_filename):
        self.logger.debug('generate_diameters')
        d = Diameters(self.genome_assembly_metadata,
                      distance_table, 
                      species_taxon_filename,
                      os.path.join(self.output_directory, 'species_taxon_output.csv'),
                      os.path.join(self.output_directory, 'min_inter_output.csv'),
                      self.verbose)
        d.calculate_diameters()
        return d

    