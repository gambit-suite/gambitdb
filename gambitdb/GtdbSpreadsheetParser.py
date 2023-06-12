# a class which will read in a GTDB spreadsheet, parse each row and output a modified spreadsheet.

import pandas
import logging
import sys

class GtdbSpreadsheetParser:
    def __init__(self, 
                 gtdb_metadata_spreadsheet, 
                 checkm_completeness, 
                 checkm_contamination, 
                 max_contigs, 
                 include_derived_samples, 
                 include_novel_species,
                 minimum_genomes_per_species,
                 species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, 
                 accessions_output_filename, 
                 debug, 
                 verbose):
        self.gtdb_metadata_spreadsheet = gtdb_metadata_spreadsheet
        self.checkm_completeness = checkm_completeness
        self.checkm_contamination = checkm_contamination
        self.max_contigs = max_contigs
        self.include_derived_samples = include_derived_samples
        self.include_novel_species = include_novel_species
        self.minimum_genomes_per_species = minimum_genomes_per_species
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accessions_output_filename = accessions_output_filename
        self.debug = debug
        self.verbose = verbose
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler(sys.stdout))
        if self.debug:
            self.logger.setLevel(logging.DEBUG)
        if self.verbose:
            self.logger.setLevel(logging.INFO)

        self.logger.debug("Debugging mode enabled")
        self.logger.debug("gtdb_metadata_spreadsheet: " + self.gtdb_metadata_spreadsheet)
        self.logger.debug("species_taxon_output_filename: " + self.species_taxon_output_filename)
        self.logger.debug("genome_assembly_metadata_output_filename: " + self.genome_assembly_metadata_output_filename)
        self.logger.debug("accessions_output_filename: " + self.accessions_output_filename)
        self.logger.debug("debug: " + str(self.debug))
        self.logger.debug("verbose: " + str(self.verbose))
        self.logger.debug("logger: " + str(self.logger))

        self.stats_starting_assemblies = 0
        self.stats_checkm_completeness = 0
        self.stats_checkm_contamination = 0
        self.stats_contig_count = 0
        self.stats_include_novel_species = 0
        self.stats_include_derived_samples = 0
        self.stats_minimum_genomes_per_species = 0 
        self.stats_genus = 0
        self.stats_starting_species = 0
        self.stats_species = 0

    def filter_input_spreadsheet(self, input_spreadsheet_df):
        initial_samples = len(input_spreadsheet_df.index)
        # filter spreadsheet to only include genomes with a checkm completeness of x% or greater
        input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['checkm_completeness'] >= self.checkm_completeness]
        self.stats_checkm_completeness = len(input_spreadsheet_df.index)
        # filter spreadsheet to only include genomes with a checkm contamination of x% or less
        input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['checkm_contamination'] <= self.checkm_contamination]
        self.stats_checkm_contamination = len(input_spreadsheet_df.index)

        # filter spreadsheet to only include genomes below the maximum number of contigs on the contig_count column
        input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['contig_count'] <= self.max_contigs]
        self.stats_contig_count = len(input_spreadsheet_df.index)

        # filter spreadsheet so that if the gtdb_taxonomy column ends with ' sp' followed by digits, then remove the row
        # These are novel species that GTDB has made up that dont exist in NCBI.
        if not self.include_novel_species:
            input_spreadsheet_df = input_spreadsheet_df[~input_spreadsheet_df['gtdb_taxonomy'].str.contains(' sp\d+$')]
        self.stats_include_novel_species =  len(input_spreadsheet_df.index)

        # if include_derived_samples is False then only include rows with 'none' from ncbi_genome_category
        if not self.include_derived_samples:
            input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['ncbi_genome_category'] == 'none']
        self.stats_include_derived_samples = len(input_spreadsheet_df.index)
        return input_spreadsheet_df

    def generate_spreadsheets(self):
        # This is a massive spreadsheet with all the GTDB metadata
        input_spreadsheet_df = self.read_in_gtdb_spreadsheet()
        accessions_spreadsheet_df = self.generate_accessions_df()

        # Get statistics for what gets passed in before filtering.
        species_taxon_ids = self.create_mock_taxon_ids_for_species(input_spreadsheet_df['species'])
        self.stats_starting_species = len(species_taxon_ids)

        # filter spreadsheet to only include genomes with a checkm completeness of 95% or greater
        input_spreadsheet_df = self.filter_input_spreadsheet(input_spreadsheet_df)

        # copy the accession column from the input spreadsheet to the accessions spreadsheet
        accessions_spreadsheet_df['assembly_accession'] = input_spreadsheet_df['accession']
        accessions_spreadsheet_df['uuid'] = input_spreadsheet_df['accession']
        accessions_spreadsheet_df['species'] = input_spreadsheet_df['species']

        species_taxon_ids = self.create_mock_taxon_ids_for_species(input_spreadsheet_df['species'])
        self.stats_species = len(species_taxon_ids)
        accessions_spreadsheet_df['species_taxid'] = input_spreadsheet_df['species'].map(species_taxon_ids)
        # count the number of rows in accessions_spreadsheet_df with each species_taxid and remove rows where the species_taxid is less than 3
        species_taxid_counts = accessions_spreadsheet_df['species_taxid'].value_counts()
        accessions_spreadsheet_df = accessions_spreadsheet_df[accessions_spreadsheet_df['species_taxid'].isin(species_taxid_counts.index[species_taxid_counts >= self.minimum_genomes_per_species])]
        self.stats_minimum_genomes_per_species = len(accessions_spreadsheet_df.index)

        # output the accessions_spreadsheet_df to a csv file saving as self.genome_assembly_metadata_output_filename
        accessions_spreadsheet_df.to_csv(self.genome_assembly_metadata_output_filename, index=False)

        # now create the species taxon table
        species_taxon_spreadsheet_df = self.generate_species_taxon_df()
        # add the species_taxon_ids values to the species_taxon_spreadsheet_df
        species_taxon_spreadsheet_df['species_taxid'] = species_taxon_ids.values()
        species_taxon_spreadsheet_df['name'] = species_taxon_ids.keys()
        species_taxon_spreadsheet_df['rank'] = 'species'
        species_taxon_spreadsheet_df['report'] = 1

        # Create mock parent_ids based on the genus from the first word of the species name
        species_taxon_spreadsheet_df['genus'] = species_taxon_spreadsheet_df['name'].str.split(' ').str[0]
        genus_taxon_ids = self.create_mock_taxon_ids_for_genus(species_taxon_spreadsheet_df['genus'], len(species_taxon_ids) + 1)
        species_taxon_spreadsheet_df['parent_taxid'] = species_taxon_spreadsheet_df['genus'].map(genus_taxon_ids)
        # delete the genus column
        del species_taxon_spreadsheet_df['genus']
        self.stats_genus = len(genus_taxon_ids)

        species_taxon_spreadsheet_df['ncbi_taxid'] = species_taxon_spreadsheet_df['species_taxid']
        species_taxon_spreadsheet_df['gambit_taxid'] = species_taxon_spreadsheet_df['species_taxid']

        # output the species_taxon_spreadsheet_df to a csv file saving as self.species_taxon_output_filename
        species_taxon_spreadsheet_df.to_csv(self.species_taxon_output_filename, index=False)

        # create a file called self.accessions_output_filename and output the column from accessions_spreadsheet_df['assembly_accession'] to the file
        accessions_spreadsheet_df.to_csv(self.accessions_output_filename, index=False)
        accessions_spreadsheet_df['assembly_accession'].to_csv(self.accessions_output_filename, header=False, columns=['assembly_accession'], index=False)
   
        self.print_summary_statistics()

    # create a dataframe for storing the species taxons with the following columns: species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid
    def generate_species_taxon_df(self):
        self.logger.debug("generate_species_taxon_df")

        # Create a new pandas dataframe with the following columns: species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid
        species_taxon_spreadsheet = pandas.DataFrame(columns=['species_taxid', 'name', 'rank', 'parent_taxid', 'ncbi_taxid', 'gambit_taxid'])
        return species_taxon_spreadsheet

    def create_mock_taxon_ids_for_species(self, species_list):
        unique_species_names = self.get_unique_species(species_list) 
        # create a dictionary with each species name as the key and an incremented integer, starting at 1, as the value
        species_taxon_ids = {species: i for i, species in enumerate(unique_species_names, 1)}
        return species_taxon_ids
    
    def create_mock_taxon_ids_for_genus(self, genus_list, offset):
        unique_genus_names = self.get_unique_genus(genus_list) 
        # create a dictionary with each genus name as the key and an incremented integer, starting at 1, as the value
        genus_taxon_ids = {genus: i + offset for i, genus in enumerate(unique_genus_names, 1)}
        return genus_taxon_ids
    
    def get_unique_genus(self, genus_list):
        self.logger.debug("get_unique_genus")

        # remove any empty values from the list
        genus_list = [x for x in genus_list if str(x) != 'nan']

        # remove any duplicate values from the list
        genus_list = list(dict.fromkeys(genus_list))

        return genus_list

    # Given a list of species names in a dataframe column, return a list of unique species names
    def get_unique_species(self, species_list):
        self.logger.debug("get_unique_species")

        # remove any empty values from the list
        species_list = [x for x in species_list if str(x) != 'nan']

        # remove any duplicate values from the list
        species_list = list(dict.fromkeys(species_list))

        return species_list

    # create a new spreadsheet with the following columns: uuid,species_taxid,assembly_accession
    def generate_accessions_df(self):
        self.logger.debug("generate_accessions_df")

        # Create a new pandas dataframe with the following columns: uuid,species_taxid,assembly_accession
        accessions_spreadsheet = pandas.DataFrame(columns=['uuid', 'species_taxid', 'assembly_accession','species'])
        return accessions_spreadsheet

    # This will read in a TSV file and return the contents as a pandas DataFrame
    def read_in_gtdb_spreadsheet(self):
        self.logger.debug("read_in_gtdb_spreadsheet")

        # Read in the GTDB spreadsheet
        spreadsheet = pandas.read_csv(self.gtdb_metadata_spreadsheet, sep='\t', index_col=False)
        spreadsheet = spreadsheet.fillna('')
        spreadsheet = spreadsheet.replace('NA', '')
        spreadsheet = spreadsheet.replace('NaN', '')
        spreadsheet = spreadsheet.replace('nan', '')

        # update all values in the accession column to remove the database prefixes
        spreadsheet['accession'] = spreadsheet['accession'].str.replace('GB_', '')
        spreadsheet['accession'] = spreadsheet['accession'].str.replace('RS_', '')

        # create a new column called species which will be formed from the gtdb_taxonomy column where the text after s__ is the species
        spreadsheet['species'] = spreadsheet['gtdb_taxonomy'].str.extract(r's__([a-zA-Z0-9_\-\s]+)', expand=False)

        self.starting_assemblies = len(spreadsheet.index)
        return spreadsheet

    # move this functionality to a separate class
    def print_summary_statistics(self):
        print("Starting Assemblies:\t" + str(self.starting_assemblies))
        print("After checkm_completeness:\t" + str(self.stats_checkm_completeness))
        print("After checkm_contamination:\t" + str(self.stats_checkm_contamination))
        print("After contig_count:\t" + str(self.stats_contig_count))
        print("After include_novel_species:\t" + str(self.stats_include_novel_species))
        print("After include_derived_samples:\t" + str(self.stats_include_derived_samples))
        print("After removing species with low numbers of genomes:\t" + str(self.stats_minimum_genomes_per_species))

        print("No. Genus:\t" + str(self.stats_genus))
        print("Starting species:\t" + str(self.stats_starting_species))
        print("Final Species:\t" + str(self.stats_species))