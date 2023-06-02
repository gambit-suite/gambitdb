# a class which will read in a GTDB spreadsheet, parse each row and output a modified spreadsheet.

import pandas
import logging
import sys

class GtdbSpreadsheetParser:
    def __init__(self, gtdb_metadata_spreadsheet, checkm_completeness, max_contigs, species_taxon_output_filename, genome_assembly_metadata_output_filename, accessions_output_filename, debug, verbose):
        self.gtdb_metadata_spreadsheet = gtdb_metadata_spreadsheet
        self.checkm_completeness = checkm_completeness
        self.max_contigs = max_contigs
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
        self.logger.info("Verbose mode enabled")
        self.logger.debug("gtdb_metadata_spreadsheet: " + self.gtdb_metadata_spreadsheet)
        self.logger.debug("species_taxon_output_filename: " + self.species_taxon_output_filename)
        self.logger.debug("genome_assembly_metadata_output_filename: " + self.genome_assembly_metadata_output_filename)
        self.logger.debug("accessions_output_filename: " + self.accessions_output_filename)
        self.logger.debug("debug: " + str(self.debug))
        self.logger.debug("verbose: " + str(self.verbose))
        self.logger.debug("logger: " + str(self.logger))

    def filter_input_spreadsheet(self, input_spreadsheet_df):
        # filter spreadsheet to only include genomes with a checkm completeness of 95% or greater
        input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['checkm_completeness'] >= self.checkm_completeness]
        # filter spreadsheet to only include genomes below the maximum number of contigs on the contig_count column
        input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['contig_count'] <= self.max_contigs]
        # filter spreadsheet so that if the gtdb_taxonomy column ends with ' sp' followed by digits, then remove the row
        # These are novel species that GTDB has made up that dont exist in NCBI.
        input_spreadsheet_df = input_spreadsheet_df[~input_spreadsheet_df['gtdb_taxonomy'].str.contains(' sp\d+$')]
        return input_spreadsheet_df

    def generate_spreadsheets(self):
        # This is a massive spreadsheet with all the GTDB metadata
        input_spreadsheet_df = self.read_in_gtdb_spreadsheet()
        accessions_spreadsheet_df = self.generate_accessions_df()

        # filter spreadsheet to only include genomes with a checkm completeness of 95% or greater
        input_spreadsheet_df = self.filter_input_spreadsheet(input_spreadsheet_df)

        # copy the accession column from the input spreadsheet to the accessions spreadsheet
        accessions_spreadsheet_df['assembly_accession'] = input_spreadsheet_df['accession']
        accessions_spreadsheet_df['uuid'] = input_spreadsheet_df['accession']
        accessions_spreadsheet_df['species'] = input_spreadsheet_df['species']

        species_taxon_ids = self.create_mock_taxon_ids_for_species(input_spreadsheet_df['species'])
        accessions_spreadsheet_df['species_taxid'] = input_spreadsheet_df['species'].map(species_taxon_ids)

        # output the accessions_spreadsheet_df to a csv file saving as self.genome_assembly_metadata_output_filename
        accessions_spreadsheet_df.to_csv(self.genome_assembly_metadata_output_filename, index=False)

        # now create the species taxon table
        species_taxon_spreadsheet_df = self.generate_species_taxon_df()
        # add the species_taxon_ids values to the species_taxon_spreadsheet_df
        species_taxon_spreadsheet_df['species_taxid'] = species_taxon_ids.values()
        species_taxon_spreadsheet_df['name'] = species_taxon_ids.keys()
        species_taxon_spreadsheet_df['rank'] = 'species'

        # Create mock parent_ids based on the genus from the first word of the species name
        species_taxon_spreadsheet_df['genus'] = species_taxon_spreadsheet_df['name'].str.split(' ').str[0]
        genus_taxon_ids = self.create_mock_taxon_ids_for_genus(species_taxon_spreadsheet_df['genus'], len(species_taxon_ids) + 1)
        species_taxon_spreadsheet_df['parent_taxid'] = species_taxon_spreadsheet_df['genus'].map(genus_taxon_ids)
        # delete the genus column
        del species_taxon_spreadsheet_df['genus']

        species_taxon_spreadsheet_df['ncbi_taxid'] = species_taxon_spreadsheet_df['species_taxid']
        species_taxon_spreadsheet_df['gambit_taxid'] = species_taxon_spreadsheet_df['species_taxid']

        # output the species_taxon_spreadsheet_df to a csv file saving as self.species_taxon_output_filename
        species_taxon_spreadsheet_df.to_csv(self.species_taxon_output_filename, index=False)

        # create a file called self.accessions_output_filename and output the column from accessions_spreadsheet_df['assembly_accession'] to the file
        accessions_spreadsheet_df.to_csv(self.accessions_output_filename, index=False)
        accessions_spreadsheet_df['assembly_accession'].to_csv(self.accessions_output_filename, header=False, columns=['assembly_accession'], index=False)
   
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
        spreadsheet['species'] = spreadsheet['gtdb_taxonomy'].str.extract(r's__([a-zA-Z0-9_\s]+)', expand=False)

        return spreadsheet
