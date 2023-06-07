# A class to generate a Gambit database

import os
import logging

from gambitdb.PairwiseTable import PairwiseTable
from gambitdb.Diameters import Diameters
from gambitdb.Curate import Curate

class GambitDb:

    def __init__(self, output_directory, assembly_directory, genome_assembly_metadata, 
                 species_taxon_filename,species_to_remove, accessions_to_remove, 
                 accession_removed_output_filename, species_removed_output_filename , 
                 signatures_output_filename, database_output_filename, species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, kmer, kmer_prefix, minimum_ngenomes, cpus, 
                 small_cluster_ngenomes, small_cluster_diameter, maximum_diameter, minimum_cluster_size, verbose):
        self.logger = logging.getLogger(__name__)
        self.output_directory = output_directory
        self.assembly_directory = assembly_directory
        self.genome_assembly_metadata = genome_assembly_metadata
        self.species_taxon_filename = species_taxon_filename
        self.species_to_remove = species_to_remove
        self.accessions_to_remove = accessions_to_remove
        self.accession_removed_output_filename = accession_removed_output_filename
        self.species_removed_output_filename = species_removed_output_filename
        self.signatures_output_filename = signatures_output_filename
        self.database_output_filename = database_output_filename
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.kmer = kmer
        self.kmer_prefix = kmer_prefix
        self.minimum_ngenomes = minimum_ngenomes

        # This file may not exist at this point
        self.species_taxon_filename = species_taxon_filename 

        self.check_output_directory_exists_or_create_one()
        self.cpus = cpus
        self.small_cluster_ngenomes = small_cluster_ngenomes
        self.small_cluster_diameter = small_cluster_diameter
        self.maximum_diameter = maximum_diameter
        self.minimum_cluster_size = minimum_cluster_size
        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    # Check if the output_directory exists, if not create it
    def check_output_directory_exists_or_create_one(self):
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

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

    def intermediate_species_taxon_filename(self):
        return os.path.join(self.output_directory, 'species_taxon_output.csv')

    def curated_species_taxon_filename(self):
        return os.path.join(self.output_directory, 'species_taxon_output_curated.csv')
    
    # generate the pairwise tables, then curate the input genomes and species, and do a second pass to generatethe pairwise tables and signatures database
    def generate_gambit_db(self):
        self.logger.debug('generate_gambit_db')
        pairwise = self.generate_pairwise_table(None)

        self.species_taxon_filename = self.check_species_taxonid_file_exists_or_create_one(self.species_taxon_filename, self.genome_assembly_metadata)
        diameters = self.generate_diameters(pairwise.distance_table_output_filename, 
                                            self.species_taxon_filename,  
                                            self.intermediate_species_taxon_filename())

        # Curate the outputs
        Curate( self.intermediate_species_taxon_filename(), 
                self.genome_assembly_metadata,
                self.assembly_directory,
                pairwise.distance_table_output_filename,
                self.species_to_remove,
                self.accessions_to_remove,
                self.curated_species_taxon_filename(),
                os.path.join(self.output_directory,self.genome_assembly_metadata_output_filename),
                os.path.join(self.output_directory,self.accession_removed_output_filename),
                os.path.join(self.output_directory,self.species_removed_output_filename),
                self.minimum_ngenomes,
                self.small_cluster_ngenomes,
                self.small_cluster_diameter,
                self.maximum_diameter,
                self.minimum_cluster_size,
                self.verbose,
                self.verbose).filter_spreadsheets_and_output_new_files()
        
        # Generate the pairwise table and signatures database again
        # need to actually filter the genomes passed in so that ignored genomes arent used.
        pairwise = self.generate_pairwise_table(os.path.join(self.output_directory, self.accession_removed_output_filename))
        diameters = self.generate_diameters(pairwise.distance_table_output_filename, 
                                            self.curated_species_taxon_filename(),
                                            os.path.join(self.output_directory, self.species_taxon_output_filename))

    def generate_pairwise_table(self, accessions_to_ignore_file):
        self.logger.debug('generate_pairwise_table')
    
        pw = PairwiseTable(self.assembly_directory, 
                           os.path.join(self.output_directory, self.signatures_output_filename), 
                           os.path.join(self.output_directory, 'pw-dists.csv'), 
                           self.kmer, 
                           self.kmer_prefix, 
                           accessions_to_ignore_file,
                           self.cpus,
                           self.verbose)
        pw.generate_sigs_and_pairwise_table()
        return pw

    def generate_diameters(self, distance_table, species_taxon_filename, species_taxon_output):
        self.logger.debug('generate_diameters')
        d = Diameters(self.genome_assembly_metadata,
                      distance_table, 
                      species_taxon_filename,
                      species_taxon_output,
                      os.path.join(self.output_directory, 'min_inter_output.csv'),
                      self.verbose)
        d.calculate_diameters()
        return d
