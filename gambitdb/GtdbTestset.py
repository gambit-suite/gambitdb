#GtdbTestset(options.species_taxon_file,
            # options.assembly_metadata_file,
            # options.gtdb_metadata_file,
            # options.output_assembly_list_filename,
            # options.output_assembly_to_species_filename,
            # options.max_genomes_per_species,
            # options.debug,
            # options.verbose ).generate_test_set()
# This is a class that will read in the species_taxon.csv, genome_assembly_metadata.csv, and GTDB spreadsheet. 
# It will read in the GTDB spreadsheet using pandas, and a new dataframe (testset) produced with two columns, the accession and the species name which is contained in the gtdb_taxonomy column as ';s__XXXX XXXXX;' where XXXX XXXXX is the species name. Only species names found in species_taxon.csv should be included in this dataframe. Any accession numbers appearing in genome_assembly_metadata.csv will be removed from the testset. 
# The testset dataframe will then be grouped by the species and filtered so that only the first n genomes are kept, where n is the max_genomes_per_species parameter. These accessions will be written to the output_assembly_list_filename file, with one assembly per line. Additionally the whole data frame will be written to the output_assembly_to_species_filename file, with one assembly per line, and the species name in the second column.
# 


# The species_taxon.csv file has the following format:
# species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid,report,diameter,ngenomes
# 1,Fusobacterium nucleatum,species,9341,1,1,1,0.3559,7
# 2,Xanthomonas oryzae,species,9342,2,2,1,0.0,1
# 3,Helicobacter pylori,species,9343,3,3,1,0.6962,135
# 5,Phytoplasma sp000009845,species,9345,5,5,1,0.5146,6
#  
# The genome_assembly_metadata.csv file has the following format:
# assembly_accession,uuid,species_taxid,species
# GCA_000007325.1,GCA_000007325.1,1,Fusobacterium nucleatum
# GCA_000007385.1,GCA_000007385.1,2,Xanthomonas oryzae subspecies 1
# GCA_000008525.1,GCA_000008525.1,3,Helicobacter pylori

# The GTDB spreadsheet has the following format:
# accession       ambiguous_bases checkm_completeness     checkm_contamination    checkm_marker_count     checkm_marker_lineage   checkm_marker_set_count        checkm_strain_heterogeneity     coding_bases    coding_density  contig_count    gc_count        gc_percentage genome_size      gtdb_genome_representative      gtdb_representative     gtdb_taxonomy   gtdb_type_designation_ncbi_taxa gtdb_type_designation_ncbi_taxa_sources        gtdb_type_species_of_genus      l50_contigs     l50_scaffolds   longest_contig  longest_scaffold      lsu_23s_contig_len       lsu_23s_count   lsu_23s_length  lsu_23s_query_id        lsu_5s_contig_len       lsu_5s_count    lsu_5s_length lsu_5s_query_id  lsu_silva_23s_blast_align_len   lsu_silva_23s_blast_bitscore    lsu_silva_23s_blast_evalue      lsu_silva_23s_blast_perc_identity      lsu_silva_23s_blast_subject_id  lsu_silva_23s_taxonomy  mean_contig_length      mean_scaffold_length    mimag_high_quality     mimag_low_quality       mimag_medium_quality    n50_contigs     n50_scaffolds   ncbi_assembly_level     ncbi_assembly_name    ncbi_assembly_type       ncbi_bioproject ncbi_biosample  ncbi_contig_count       ncbi_contig_n50 ncbi_country    ncbi_date       ncbi_genbank_assembly_accession        ncbi_genome_category    ncbi_genome_representation      ncbi_isolate    ncbi_isolation_source   ncbi_lat_lon   ncbi_molecule_count     ncbi_ncrna_count        ncbi_organism_name      ncbi_protein_count      ncbi_refseq_category    ncbi_rrna_count        ncbi_scaffold_count     ncbi_scaffold_l50       ncbi_scaffold_n50       ncbi_scaffold_n75       ncbi_scaffold_n90     ncbi_seq_rel_date        ncbi_spanned_gaps       ncbi_species_taxid      ncbi_ssu_count  ncbi_strain_identifiers ncbi_submitter  ncbi_taxid     ncbi_taxonomy   ncbi_taxonomy_unfiltered        ncbi_total_gap_length   ncbi_total_length       ncbi_translation_table  ncbi_trna_count        ncbi_type_material_designation  ncbi_ungapped_length    ncbi_unspanned_gaps     ncbi_wgs_master protein_count   scaffold_count ssu_contig_len  ssu_count       ssu_gg_blast_align_len  ssu_gg_blast_bitscore   ssu_gg_blast_evalue     ssu_gg_blast_perc_identity     ssu_gg_blast_subject_id ssu_gg_taxonomy ssu_length      ssu_query_id    ssu_silva_blast_align_len       ssu_silva_blast_bitscore       ssu_silva_blast_evalue  ssu_silva_blast_perc_identity   ssu_silva_blast_subject_id      ssu_silva_taxonomy      total_gap_length       trna_aa_count   trna_count      trna_selenocysteine_count
# GB_GCA_000006155.2      1916    93.12   0       1171    g__Bacillus (UID902)    324     0       4305660 80.1789924134926        426   1869446  35.10140364607  5370060 RS_GCF_000742895.1      f       d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A anthracis       not type material       none    f       63      1       181677  5093554 none    0       none  none     none    0       none    none    none    none    none    none    none    none    12506   1790020 f       f       t       25045 5093554  Chromosome      ASM615v2        na      PRJNA299        SAMN02435829    426     25045   none    2002-05-09      GCA_000006155.2none    full    none    none    none    3       none    Bacillus anthracis str. A2012   297     na      none    3       1       50935545093554 5093554 2002/05/09      423     1392    none    A2012   J. Craig Venter Institute       191218  d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus anthracis    d__Bacteria;x__Terrabacteria group;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;x__Bacillus cereus group;s__Bacillus anthracis;x__Bacillus anthracis str. A2012   42300   5370060none    none    none    5327760 0       AAAC00000000.1  5745    3       none    0       none    none    none    none    none    none  none     none    none    none    none    none    none    none    42300   16      31      0
# GB_GCA_000007325.1      1       99.95   0       149     k__Bacteria (UID2329)   89      0       1973459 90.7546102552311        1     590411   27.1515875610888        2174500 GB_GCA_000007325.1      t       d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__Fusobacterium nucleatum      type strain of species  LPSN    t       1       1       21745002174500 2174500 5       2898    AE009951.2      2174500 5       107     AE009951.2      2898    5352    0       100     AE009951.1073886.1076812       Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatum subsp. nucleatum ATCC 25586       2174500 2174500 t       f       f       2174500 2174500 Complete Genome ASM732v1        na      PRJNA295      SAMN02603417     none    none    none    2002-04-09      GCA_000007325.1 none    full    none    none    none    1       none    Fusobacterium nucleatum subsp. nucleatum ATCC 25586    none    na      none    1       1       2174500 2174500 2174500 2002/04/09      0     851      none    ATCC 25586      Integrated Genomics     190304  d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__Fusobacterium nucleatum        d__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__Fusobacterium nucleatum;sb__Fusobacterium nucleatum subsp. nucleatum;x__Fusobacterium nucleatum subsp. nucleatum ATCC 25586   0       2174500 none    none    none    2174500 0       none    2022    1       2174500 5       none    none  none     none    none    none    1517    AE009951.2      1517    2802    0       100     CP028109.679747.681305  Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium nucleatum subsp. nucleatum ATCC 23726       0       19    47       0
import logging
import pandas as pd
import sys

class GtdbTestset:
    def __init__(self, species_taxon_file,
            assembly_metadata_file,
            gtdb_metadata_file,
            output_assembly_list_filename,
            output_assembly_to_species_filename,
            max_genomes_per_species,
            debug,
            verbose):
        self.species_taxon_file = species_taxon_file
        self.assembly_metadata_file = assembly_metadata_file
        self.gtdb_metadata_file = gtdb_metadata_file
        self.output_assembly_list_filename = output_assembly_list_filename
        self.output_assembly_to_species_filename = output_assembly_to_species_filename
        self.max_genomes_per_species = max_genomes_per_species
        self.debug = debug
        self.verbose = verbose

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logging.StreamHandler(sys.stdout))
        if self.debug:
            self.logger.setLevel(logging.DEBUG)
        if self.verbose:
            self.logger.setLevel(logging.INFO)

        self.species_taxon_df = None
        self.assembly_metadata_df = None
        self.gtdb_metadata_df = None
        self.assembly_list_df = None
        self.assembly_to_species_df = None

    def generate_test_set(self):
        self.species_taxon_df = pd.read_csv(self.species_taxon_file, sep=",")
        self.assembly_metadata_df = pd.read_csv(self.assembly_metadata_file, sep=",")
        self.gtdb_metadata_df = pd.read_csv(self.gtdb_metadata_file, sep="\t")

        self.logger.info("Filtering assembly metadata file to only include genomes with a GTDB taxonomy")
        self.gtdb_metadata_df = self.gtdb_metadata_df[self.gtdb_metadata_df["gtdb_taxonomy"].notnull()]
        # update all values in the accession column to remove the database prefixes
        self.gtdb_metadata_df['accession'] = self.gtdb_metadata_df['accession'].str.replace('GB_', '')
        self.gtdb_metadata_df['accession'] = self.gtdb_metadata_df['accession'].str.replace('RS_', '')
        # create a new column called species which will be formed from the gtdb_taxonomy column where the text after s__ is the species
        self.gtdb_metadata_df['species'] = self.gtdb_metadata_df['gtdb_taxonomy'].str.extract(r's__([a-zA-Z0-9_\-\s]+)', expand=False)

        # copy the accession column and the species column to a new data frame stored in self.assembly_to_species_df. rename the accession column to assembly_accession
        self.assembly_to_species_df = self.gtdb_metadata_df[['accession', 'species']].copy()
        self.assembly_to_species_df.rename(columns={'accession': 'assembly_accession'}, inplace=True)

        # get a list of all the species in the species_taxon_df, where the rank is 'species' and the species name does not contain 'subspecies'
        species_list = self.species_taxon_df[(self.species_taxon_df['rank'] == 'species') & (~self.species_taxon_df['name'].str.contains('subspecies'))]['name'].tolist()

        # remove any rows from the assembly_to_species_df where the species is not in the species_list
        self.assembly_to_species_df = self.assembly_to_species_df[self.assembly_to_species_df['species'].isin(species_list)]

        # get a list of all assembly_accession from assembly_metadata_df
        assembly_accession_list = self.assembly_metadata_df['assembly_accession'].tolist()
        # remove any rows from the assembly_to_species_df where the assembly_accession is in the assembly_accession_list
        self.assembly_to_species_df = self.assembly_to_species_df[~self.assembly_to_species_df['assembly_accession'].isin(assembly_accession_list)]

        self.assembly_to_species_df = self.assembly_to_species_df.groupby('species').head(self.max_genomes_per_species)

      
        self.create_output_files()

    def create_output_files(self):
        self.assembly_to_species_df.to_csv(self.output_assembly_to_species_filename, sep=",", index=False)
        self.logger.info("Created output file: %s", self.output_assembly_to_species_filename)

        # output the assembly_accession column from self.assembly_to_species_df to a csv file
        self.assembly_list_df = self.assembly_to_species_df[['assembly_accession']].copy()
        self.assembly_list_df.to_csv(self.output_assembly_list_filename, sep=",", index=False, header=False)
        self.logger.info("Created output file: %s", self.output_assembly_list_filename)
