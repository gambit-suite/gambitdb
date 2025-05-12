# a class which will read in a GTDB spreadsheet, parse each row and output a modified spreadsheet.

import pandas
import logging
import sys

class GtdbSpreadsheetParser:
    """
    Reads in a GTDB spreadsheet, parses each row and outputs a modified spreadsheet.
    """
    def __init__(self, 
                 gtdb_metadata_spreadsheet, 
                 checkm_completeness, 
                 checkm_contamination, 
                 max_contigs, 
                 include_derived_samples, 
                 include_novel_species,
                 minimum_genomes_per_species,
                 species_filter,
                 species_taxon_output_filename, 
                 genome_assembly_metadata_output_filename, 
                 accessions_output_filename, 
                 representative_genomes_filename,
                 debug, 
                 verbose):
        """
    Initializes the GtdbSpreadsheetParser class.
    Args:
      gtdb_metadata_spreadsheet (str): The path to the GTDB metadata spreadsheet.
      checkm_completeness (float): The minimum CheckM completeness to include a genome.
      checkm_contamination (float): The maximum CheckM contamination to include a genome.
      max_contigs (int): The maximum number of contigs to include a genome.
      include_derived_samples (bool): Whether to include derived samples.
      include_novel_species (bool): Whether to include novel species.
      minimum_genomes_per_species (int): The minimum number of genomes per species to include a species.
      species_taxon_output_filename (str): The path to the output file for species taxon information.
      genome_assembly_metadata_output_filename (str): The path to the output file for genome assembly metadata.
      accessions_output_filename (str): The path to the output file for accessions.
      debug (bool): Whether to enable debugging mode.
      verbose (bool): Whether to enable verbose mode.
    Side Effects:
      Sets the following attributes:
        gtdb_metadata_spreadsheet
        checkm_completeness
        checkm_contamination
        max_contigs
        include_derived_samples
        include_novel_species
        minimum_genomes_per_species
        species_taxon_output_filename
        genome_assembly_metadata_output_filename
        accessions_output_filename
        debug
        verbose
        logger
        stats_starting_assemblies
        stats_checkm_completeness
        stats_checkm_contamination
        stats_contig_count
        stats_include_novel_species
        stats_include_derived_samples
        stats_minimum_genomes_per_species
        stats_genus
        stats_starting_species
        stats_species
    Examples:
      >>> parser = GtdbSpreadsheetParser(
          gtdb_metadata_spreadsheet="metadata.tsv",
          checkm_completeness=0.5,
          checkm_contamination=0.1,
          max_contigs=1000,
          include_derived_samples=True,
          include_novel_species=True,
          minimum_genomes_per_species=2,
          species_taxon_output_filename="species_taxon.tsv",
          genome_assembly_metadata_output_filename="genome_assembly_metadata.tsv",
          accessions_output_filename="accessions.tsv",
          debug=False,
          verbose=True
      )
    """
        self.gtdb_metadata_spreadsheet = gtdb_metadata_spreadsheet
        self.checkm_completeness = checkm_completeness
        self.checkm_contamination = checkm_contamination
        self.max_contigs = max_contigs
        self.include_derived_samples = include_derived_samples
        self.include_novel_species = include_novel_species
        self.minimum_genomes_per_species = minimum_genomes_per_species
        self.species_filter = species_filter
        self.species_taxon_output_filename = species_taxon_output_filename
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accessions_output_filename = accessions_output_filename
        self.representative_genomes_filename = representative_genomes_filename
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

        self.representative_genomes = []

    def filter_input_spreadsheet(self, input_spreadsheet_df):
        """
    Filters a given spreadsheet to only include genomes with a checkm completeness of x% or greater, a checkm contamination of x% or less, and below the maximum number of contigs.
    Args:
      input_spreadsheet_df (DataFrame): The input spreadsheet to filter.
    Returns:
      DataFrame: The filtered spreadsheet.
    Examples:
      >>> filter_input_spreadsheet(input_spreadsheet_df)
      Filtered spreadsheet
    """

        if self.species_filter != '':
            input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['species'] == self.species_filter]
            
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
            input_spreadsheet_df = input_spreadsheet_df[~input_spreadsheet_df['gtdb_taxonomy'].str.contains(r' sp\d+$')]
        self.stats_include_novel_species =  len(input_spreadsheet_df.index)

        # if include_derived_samples is False then only include rows with 'none' from ncbi_genome_category
        if not self.include_derived_samples:
            input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['ncbi_genome_category'] == 'none']
        self.stats_include_derived_samples = len(input_spreadsheet_df.index)
        return input_spreadsheet_df

    def generate_spreadsheets(self):
        """
    Generates spreadsheets for storing the species taxons and accessions.
    Args:
      None
    Returns:
      None
    Examples:
      >>> generate_spreadsheets()
      Generates spreadsheets
    """
        # This is a massive spreadsheet with all the GTDB metadata
        input_spreadsheet_df = self.read_in_gtdb_spreadsheet()
        accessions_spreadsheet_df = self.generate_accessions_df()

        self.representative_genomes = self.get_representative_genomes(input_spreadsheet_df)
        if self.representative_genomes_filename:
          self.save_representative_genome_accessions_to_file(self.representative_genomes_filename)

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

        # create a data hash where the key is the 'species_taxid' and the value is the number of genomes with that species_taxid
        species_taxid_counts = accessions_spreadsheet_df['species_taxid'].value_counts()
        # add a column called ngenomes to species_taxon_spreadsheet_df and add the corresponding species_taxid_counts value for each species_taxid
        species_taxon_spreadsheet_df['ngenomes'] = species_taxon_spreadsheet_df['species_taxid'].map(species_taxid_counts)
        

        # output the species_taxon_spreadsheet_df to a csv file saving as self.species_taxon_output_filename
        species_taxon_spreadsheet_df.to_csv(self.species_taxon_output_filename, index=False)

        # create a file called self.accessions_output_filename and output the column from accessions_spreadsheet_df['assembly_accession'] to the file
        accessions_spreadsheet_df.to_csv(self.accessions_output_filename, index=False)
        accessions_spreadsheet_df['assembly_accession'].to_csv(self.accessions_output_filename, header=False, columns=['assembly_accession'], index=False)
   
        if self.verbose:
          self.print_summary_statistics()

    # create a dataframe for storing the species taxons with the following columns: species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid
    def generate_species_taxon_df(self):
        """
    Creates a dataframe for storing the species taxons with the following columns: species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid.
    Args:
      None
    Returns:
      DataFrame: A dataframe for storing the species taxons.
    Examples:
      >>> generate_species_taxon_df()
      Species taxon dataframe
    """
        self.logger.debug("generate_species_taxon_df")

        # Create a new pandas dataframe with the following columns: species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid
        species_taxon_spreadsheet = pandas.DataFrame(columns=['species_taxid', 'name', 'rank', 'parent_taxid', 'ncbi_taxid', 'gambit_taxid'])
        return species_taxon_spreadsheet

    def create_mock_taxon_ids_for_species(self, species_list):
        """
    Creates a dictionary with each species name as the key and an incremented integer, starting at 1, as the value.
    Args:
      species_list (list): A list of species names.
    Returns:
      dict: A dictionary with each species name as the key and an incremented integer, starting at 1, as the value.
    Examples:
      >>> create_mock_taxon_ids_for_species(species_list)
      {'species1': 1, 'species2': 2, ...}
    """
        unique_species_names = self.get_unique_species(species_list) 
        # create a dictionary with each species name as the key and an incremented integer, starting at 1, as the value
        species_taxon_ids = {species: i for i, species in enumerate(unique_species_names, 1)}
        return species_taxon_ids
    
    def create_mock_taxon_ids_for_genus(self, genus_list, offset):
        """
    Creates a dictionary with each genus name as the key and an incremented integer, starting at 1, as the value.
    Args:
      genus_list (list): A list of genus names.
      offset (int): The offset to start the incremented integer from.
    Returns:
      dict: A dictionary with each genus name as the key and an incremented integer, starting at the given offset, as the value.
    Examples:
      >>> create_mock_taxon_ids_for_genus(genus_list, 10)
      {'genus1': 10, 'genus2': 11, ...}
    """
        unique_genus_names = self.get_unique_genus(genus_list) 
        # create a dictionary with each genus name as the key and an incremented integer, starting at 1, as the value
        genus_taxon_ids = {genus: i + offset for i, genus in enumerate(unique_genus_names, 1)}
        return genus_taxon_ids
    
    def get_unique_genus(self, genus_list):
        """
    Given a list of genus names in a dataframe column, return a list of unique genus names.
    Args:
      genus_list (list): A list of genus names.
    Returns:
      list: A list of unique genus names.
    Examples:
      >>> get_unique_genus(genus_list)
      ['genus1', 'genus2', ...]
    """
        self.logger.debug("get_unique_genus")

        # remove any empty values from the list
        genus_list = [x for x in genus_list if str(x) != 'nan']

        # remove any duplicate values from the list
        genus_list = list(dict.fromkeys(genus_list))

        return genus_list

    # Given a list of species names in a dataframe column, return a list of unique species names
    def get_unique_species(self, species_list):
        """
    Given a list of species names in a dataframe column, return a list of unique species names.
    Args:
      species_list (list): A list of species names.
    Returns:
      list: A list of unique species names.
    Examples:
      >>> get_unique_species(species_list)
      ['species1', 'species2', ...]
    """
        self.logger.debug("get_unique_species")

        # remove any empty values from the list
        species_list = [x for x in species_list if str(x) != 'nan']

        # remove any duplicate values from the list
        species_list = list(dict.fromkeys(species_list))

        return species_list

    # create a new spreadsheet with the following columns: uuid,species_taxid,assembly_accession
    def generate_accessions_df(self):
        """
    Creates a new spreadsheet with the following columns: uuid,species_taxid,assembly_accession.
    Args:
      None
    Returns:
      DataFrame: A new spreadsheet with the given columns.
    Examples:
      >>> generate_accessions_df()
      Accessions dataframe
    """
        self.logger.debug("generate_accessions_df")

        # Create a new pandas dataframe with the following columns: uuid,species_taxid,assembly_accession
        accessions_spreadsheet = pandas.DataFrame(columns=['uuid', 'species_taxid', 'assembly_accession','species'])
        return accessions_spreadsheet

    # This will read in a TSV file and return the contents as a pandas DataFrame
    def read_in_gtdb_spreadsheet(self):
        """
    Reads in a TSV file and returns the contents as a pandas DataFrame.
    Args:
      self (GtdbSpreadsheetParser): The GtdbSpreadsheetParser object.
    Returns:
      pandas.DataFrame: The contents of the TSV file.
    Side Effects:
      Updates the 'accession' column to remove the database prefixes.
      Creates a new column called 'species' which is formed from the 'gtdb_taxonomy' column.
    Examples:
      >>> parser = GtdbSpreadsheetParser()
      >>> spreadsheet = parser.read_in_gtdb_spreadsheet()
    """
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
        """
    Prints summary statistics.
    Args:
      self (GtdbSpreadsheetParser): The GtdbSpreadsheetParser object.
    Side Effects:
      Prints summary statistics to the console.
    Examples:
      >>> parser = GtdbSpreadsheetParser()
      >>> parser.print_summary_statistics()
      Starting Assemblies:    <number>
      After checkm_completeness:    <number>
      After checkm_contamination:    <number>
      After contig_count:    <number>
      After include_novel_species:    <number>
      After include_derived_samples:    <number>
      After removing species with low numbers of genomes:    <number>
      No. Genus:    <number>
      Starting species:    <number>
      Final Species:    <number>
    """
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

    def get_representative_genomes(self, input_spreadsheet_df):
      """
      Filters the input spreadsheet to get only the rows with 'ncbi_refseq_category' equal to 'representative genome',
      and returns a dataframe with the 'assembly_accession' column.
      Args:
        input_spreadsheet_df (DataFrame): The input spreadsheet dataframe.
      Returns:
        DataFrame: A dataframe with the 'assembly_accession' column.
      Examples:
        >>> get_representative_genomes(input_spreadsheet_df)
        Representative genomes dataframe
      """
      filtered_input_spreadsheet_df = input_spreadsheet_df[input_spreadsheet_df['ncbi_refseq_category'] == 'representative genome']
      representative_accessions = pandas.DataFrame(columns=['assembly_accession'])
      representative_accessions['assembly_accession'] = filtered_input_spreadsheet_df['accession']
      return representative_accessions

    def save_representative_genome_accessions_to_file(self, output_filename):
      """
      Saves the representative genome accessions to a file.
      Args:
        output_filename (str): The path to the output file.
      Side Effects:
        Saves the representative genome accessions to a file.
      Examples:
        >>> save_representative_genome_accessions_to_file("representative_genomes.tsv")
        Saves the representative genome accessions to a file
      """
      self.representative_genomes.to_csv(output_filename, index=False)
