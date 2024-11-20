# A class to generate a Gambit database

import os
import logging

from gambitdb.PairwiseTable import PairwiseTable
from gambitdb.Diameters import Diameters
from gambitdb.Curate import Curate
from gambitdb.CompressClusters import CompressClusters
import tempfile

class GambitDb:
    """
    Generates a Gambit database.
    """

    def __init__(self, output_directory, assembly_directory, genome_assembly_metadata, 
                 species_taxon_filename,species_to_remove, accessions_to_remove, 
                 accession_removed_output_filename, species_removed_output_filename , 
                 signatures_output_filename, database_output_filename, species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, kmer, kmer_prefix, minimum_ngenomes, cpus, 
                 small_cluster_ngenomes, small_cluster_diameter, maximum_diameter, minimum_cluster_size, linkage_method,
                 compress_max_distance, representative_genomes, verbose):
        """
    Initializes a GambitDb object.
    Args:
      output_directory (str): The directory to store the output files.
      assembly_directory (str): The directory containing the assembly files.
      genome_assembly_metadata (str): The path to the genome assembly metadata file.
      species_taxon_filename (str): The path to the species taxon file.
      species_to_remove (list): A list of species to remove.
      accessions_to_remove (list): A file containing a list of accessions to remove.
      accession_removed_output_filename (str): The path to the accession removed output file.
      species_removed_output_filename (str): The path to the species removed output file.
      signatures_output_filename (str): The path to the signatures output file.
      database_output_filename (str): The path to the database output file.
      species_taxon_output_filename (str): The path to the species taxon output file.
      genome_assembly_metadata_output_filename (str): The path to the genome assembly metadata output file.
      kmer (int): The kmer size.
      kmer_prefix (str): The kmer prefix.
      minimum_ngenomes (int): The minimum number of genomes.
      cpus (int): The number of CPUs to use.
      small_cluster_ngenomes (int): The minimum number of genomes in a small cluster.
      small_cluster_diameter (int): The maximum diameter of a small cluster.
      maximum_diameter (int): The maximum diameter of a cluster.
      minimum_cluster_size (int): The minimum size of a cluster.
      linkage_method (str): The linkage method.
      verbose (bool): Whether to print verbose output.
    Side Effects:
      Creates the output directory if it does not exist.
    """
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
        self.linkage_method = linkage_method
        self.compress_max_distance = compress_max_distance
        self.representative_genomes = representative_genomes
        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    # Check if the output_directory exists, if not create it
    def check_output_directory_exists_or_create_one(self):
        """
    Checks if the output directory exists, and creates it if it does not.
    Side Effects:
      Creates the output directory if it does not exist.
    """
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def check_species_taxonid_file_exists_or_create_one(self, species_taxon_filename, genome_assembly_metadata):
        """
    Checks if the species taxon file exists, and creates it if it does not.
    Args:
      species_taxon_filename (str): The path to the species taxon file.
      genome_assembly_metadata (str): The path to the genome assembly metadata file.
    Side Effects:
      Creates the species taxon file if it does not exist.
    """
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
        """
    Returns the path to the intermediate species taxon file.
    Returns:
      str: The path to the intermediate species taxon file.
    """
        return os.path.join(self.output_directory, 'species_taxon_output.csv')

    def curated_species_taxon_filename(self):
        """
    Returns the path to the curated species taxon file.
    Returns:
      str: The path to the curated species taxon file.
    """
        return os.path.join(self.output_directory, 'species_taxon_output_curated.csv')
    
    # generate the pairwise tables, then curate the input genomes and species, and do a second pass to generatethe pairwise tables and signatures database
    def generate_gambit_db(self):
        """
    Generates a Gambit database.
    Args:
      self (GambitDb): The GambitDb instance.
    Returns:
      None
    Side Effects:
      Generates a pairwise table, curates the input genomes and species, and does a second pass to generate the pairwise tables and signatures database.
    Examples:
      >>> GambitDb.generate_gambit_db()
    """
        self.logger.debug('generate_gambit_db')
        pairwise = self.generate_pairwise_table(None)

        self.species_taxon_filename = self.check_species_taxonid_file_exists_or_create_one(self.species_taxon_filename, self.genome_assembly_metadata)
        diameters = self.generate_diameters(pairwise.distance_table_output_filename, 
                                            self.species_taxon_filename,  
                                            self.intermediate_species_taxon_filename(), 
                                            self.genome_assembly_metadata)

        # Many highly similar samples are in the public databases which offer little additional information for GAMBIT. Cluster similar samples together, keep the centroid and remove the rest.
        if self.compress_max_distance is not None:
            compress = CompressClusters(pairwise.distance_table_output_filename, 
                                        os.path.join(self.output_directory, 'pw-dists-compressed.csv'), 
                                        self.compress_max_distance, 
                                        self.representative_genomes,
                                        self.verbose)
            sample_accessions_highly_similar, num_samples = compress.compress()
            # dont remove too many samples, if that happens remove some from the sample_accessions_highly_similar list
            if len(sample_accessions_highly_similar) + self.minimum_ngenomes >= num_samples:
                no_samples_to_keep = (len(sample_accessions_highly_similar) + self.minimum_ngenomes) - num_samples
                sample_accessions_highly_similar = sample_accessions_highly_similar[:-no_samples_to_keep]

            # if a file path exists for accessions_to_remove then open it for appending and write the sample_accessions_highly_similar, one per line
            # otherwise create a new file and write the sample_accessions_highly_similar, one per line
            if self.accessions_to_remove is not None:
                with open(self.accessions_to_remove, 'a') as f:
                    for accession in sample_accessions_highly_similar:
                        f.write(accession + '\n')
            else:
                with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
                  self.accessions_to_remove = f.name
                  for accession in sample_accessions_highly_similar:
                    f.write(accession + '\n')

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
                self.linkage_method,
                self.verbose,
                self.verbose).filter_spreadsheets_and_output_new_files()
        
        # Generate the pairwise table and signatures database again
        # need to actually filter the genomes passed in so that ignored genomes arent used.
        pairwise = self.generate_pairwise_table(os.path.join(self.output_directory, self.accession_removed_output_filename))
        diameters = self.generate_diameters(pairwise.distance_table_output_filename, 
                                            self.curated_species_taxon_filename(),
                                            os.path.join(self.output_directory, self.species_taxon_output_filename),
                                            os.path.join(self.output_directory,self.genome_assembly_metadata_output_filename))

    def generate_pairwise_table(self, accessions_to_ignore_file):
        """
    Generates a pairwise table.
    Args:
      self (GambitDb): The GambitDb instance.
      accessions_to_ignore_file (str): The file containing accessions to ignore.
    Returns:
      PairwiseTable: The generated pairwise table.
    Examples:
      >>> GambitDb.generate_pairwise_table(accessions_to_ignore_file)
    """
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

    def generate_diameters(self, distance_table, species_taxon_filename, species_taxon_output, genome_assembly_metadata):
        """
    Generates diameters.
    Args:
      self (GambitDb): The GambitDb instance.
      distance_table (str): The distance table.
      species_taxon_filename (str): The species taxon filename.
      species_taxon_output (str): The species taxon output.
      genome_assembly_metadata (str): The genome assembly metadata.
    Returns:
      Diameters: The generated diameters.
    Examples:
      >>> GambitDb.generate_diameters(distance_table, species_taxon_filename, species_taxon_output, genome_assembly_metadata)
    """
        self.logger.debug('generate_diameters')
        d = Diameters(genome_assembly_metadata,
                      distance_table, 
                      species_taxon_filename,
                      species_taxon_output,
                      os.path.join(self.output_directory, 'min_inter_output.csv'),
                      self.verbose)
        d.calculate_diameters()
        return d
