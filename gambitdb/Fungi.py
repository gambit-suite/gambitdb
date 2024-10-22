# a class which will read in a GTDB spreadsheet, parse each row and output a modified spreadsheet.

import pandas
import logging
import sys

class FungiParser:
    """
    Reads in a Fungi spreadsheet, parses each row and outputs a modified spreadsheet.
    """
    def __init__(self, 
                 fungi_metadata_spreadsheet, 
                 max_contigs, 
                 minimum_genomes_per_species,
                 species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, 
                 accessions_output_filename, 
                 representative_genomes_filename):
        """
    Initializes the FungiParser class.
    Args:
      fungi_metadata_spreadsheet (str): The path to the Fungi metadata spreadsheet.
      max_contigs (int): The maximum number of contigs to include a genome.
      minimum_genomes_per_species (int): The minimum number of genomes per species to include a species.
      species_taxon_output_filename (str): The path to the output file for species taxon information.
      genome_assembly_metadata_output_filename (str): The path to the output file for genome assembly metadata.
      accessions_output_filename (str): The path to the output file for accessions.
    """
        self.fungi_metadata_spreadsheet = fungi_metadata_spreadsheet
        self.max_contigs = max_contigs
        self.minimum_genomes_per_species = minimum_genomes_per_species
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accessions_output_filename = accessions_output_filename
        self.representative_genomes_filename = representative_genomes_filename

        self.representative_genomes = []

        