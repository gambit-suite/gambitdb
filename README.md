# GambitDB
This is a tool to generate a Gambit database. Its primary input is a spreadsheet from GTDB such as this [release](https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz).


# Installation
This software can be installed using pip:
```
pip install .
```
or
```
pip install git+
```

Additional dependancies are:
* GAMBIT
* ncbi-genome-download
* GAMBITtools

## Docker
You can use docker to run the software. To build from scratch run:
```
docker build -t gambitdb .
```
then to run one of the scripts:
```
docker run -v $(pwd):/data gambitdb gambitdb -h
```

# Usage
## One step database
To create a database in one step, execute the `run_gambitdb_bacteria.sh` script:

```bash
./run_gambitdb_bacteria.sh /path/to/working_directory /path/to/gtdb_spreadsheet.tsv num_cores
```
This script will create a GAMBIT database from a GTDB spreadsheet. It parses the spreadsheet, downloads data with ncbi-genome-download and outputs a GAMBIT database.

## gambitdb-gtdb

```
usage: gambitdb-gtdb [options]

Given a GTDB metadata spreadsheet, output a list of accessions to download, a species taxonid file and a genome metadata file

positional arguments:
  gtdb_metadata_spreadsheet
                        GTDB metadata file such as bac120_metadata_r214.tsv

options:
  -h, --help            show this help message and exit
  --checkm_completeness CHECKM_COMPLETENESS, -b CHECKM_COMPLETENESS
                        Minimum checkm completeness of the genome [0-100] (default: 97.0)
  --checkm_contamination CHECKM_CONTAMINATION, -c CHECKM_CONTAMINATION
                        Maximum checkm contamination of the genome [0-100] (default: 2.0)
  --max_contigs MAX_CONTIGS, -d MAX_CONTIGS
                        Maximum number of contigs. Please note some species systematically assemble poorly with short read data. (default: 100)
  --include_derived_samples, -e
                        Include mixed samples from metagenomes, environment, single cell (default: False)
  --include_novel_species, -f
                        Include novel species called sp12345. The genus must be known (default: False)
  --minimum_genomes_per_species MINIMUM_GENOMES_PER_SPECIES, -i MINIMUM_GENOMES_PER_SPECIES
                        Minimum number of genomes in a species, otherwise exclude the species (default: 2)
  --species_filter SPECIES_FILTER, -j SPECIES_FILTER
                        Only include species that match this string (default: )
  --species_taxon_output_filename SPECIES_TAXON_OUTPUT_FILENAME, -s SPECIES_TAXON_OUTPUT_FILENAME
                        Output filename for with the taxonomy (default: species_taxa.csv)
  --genome_assembly_metadata_output_filename GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME, -g GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME
                        Genome metadata (default: assembly_metadata.csv)
  --accessions_output_filename ACCESSIONS_OUTPUT_FILENAME, -a ACCESSIONS_OUTPUT_FILENAME
                        Genome accessions for download (default: accessions_to_download.csv)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb
```
usage: gambitdb [options]

Generate a Gambit database

positional arguments:
  assembly_directory    A directory containing assemblies in FASTA format
  genome_assembly_metadata
                        A CSV containing the assembly file and path, and the species taxon ID
  species_taxon_filename
                        CSV containing species taxonomy, may be generated automatically from assembly metadata if missing

options:
  -h, --help            show this help message and exit
  --species_to_remove SPECIES_TO_REMOVE, -x SPECIES_TO_REMOVE
                        Optional file containing a list of species to remove (1 per line) (default: None)
  --accessions_to_remove ACCESSIONS_TO_REMOVE, -y ACCESSIONS_TO_REMOVE
                        Optional file containing a list of accession numbers to remove (1 per line) (default: None)
  --output_directory OUTPUT_DIRECTORY, -d OUTPUT_DIRECTORY
                        Output directory (default: output_dir)
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: database.gs)
  --database_output_filename DATABASE_OUTPUT_FILENAME, -g DATABASE_OUTPUT_FILENAME
                        Output filename for core database (default: database.gdb)
  --accession_removed_output_filename ACCESSION_REMOVED_OUTPUT_FILENAME, -c ACCESSION_REMOVED_OUTPUT_FILENAME
                        Output filename for a list of accessions removed (default: accessions_removed.csv)
  --species_removed_output_filename SPECIES_REMOVED_OUTPUT_FILENAME, -w SPECIES_REMOVED_OUTPUT_FILENAME
                        Output filename for a list of species removed (default: species_removed.csv)
  --species_taxon_output_filename SPECIES_TAXON_OUTPUT_FILENAME, -t SPECIES_TAXON_OUTPUT_FILENAME
                        Output filename for a list of species taxon IDs (default: species_taxon.csv)
  --genome_assembly_metadata_output_filename GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME, -m GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME
                        Output filename for a list of genome assembly metadata (default: genome_assembly_metadata.csv)
  --kmer KMER, -k KMER  Length of the k-mer to use (default: 11)
  --kmer_prefix KMER_PREFIX, -f KMER_PREFIX
                        Kmer prefix (default: ATGAC)
  --minimum_ngenomes MINIMUM_NGENOMES, -n MINIMUM_NGENOMES
                        Minimum number of genomes for a species to be included (default: 1)
  --cpus CPUS, -p CPUS  Number of cpus to use (default: 1)
  --small_cluster_ngenomes SMALL_CLUSTER_NGENOMES
                        Minimum number of genomes for a species to be included in a small cluster, along with --small_cluster_diameter (default: 4)
  --small_cluster_diameter SMALL_CLUSTER_DIAMETER
                        Maximum diameter of a species to be included in a small cluster along with --small_cluster_ngenomes (default: 0.7)
  --maximum_diameter MAXIMUM_DIAMETER
                        The maximum diameter to allow before attempting to split a species into subspecies (default: 0.7)
  --minimum_cluster_size MINIMUM_CLUSTER_SIZE
                        After splitting a species into subspecies, this is the minimum number of genomes which must be present, otherwise the genome is removed. (default: 2)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-create
```
usage: gambitdb-create [options]

Generate a GAMBIT database. Requires preprocessed input files

positional arguments:
  genome_assembly_metadata
                        A CSV containing the assembly file and path, and the species taxon ID
  species_taxon_filename
                        A CSV containing species taxonomy
  signatures_filename   A signatures .h5 file created by gambit signatures

options:
  -h, --help            show this help message and exit
  --db_key DB_KEY       Unique key for database, no spaces (default: organisation/database)
  --db_version DB_VERSION
                        Unique version, x.y.z (default: 1.0.0)
  --db_author DB_AUTHOR
                        Name of person who created the database (default: Jane Doe)
  --db_date DB_DATE     Date database was created as YYYY-MM-DD (default: 2022-12-31)
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: database.gs)
  --database_output_filename DATABASE_OUTPUT_FILENAME, -g DATABASE_OUTPUT_FILENAME
                        Output filename for core database (default: database.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```



# Scripts for working with existing databases
## gambitdb-apply-patch
```
usage: gambitdb-apply-patch [options]

Given two GAMBIT signatures files, merge them and return a new file.

positional arguments:
  signatures_main_filename
                        A signatures .h5 file created by gambit signatures
  database_main_filename
                        An SQLite database
  signatures_patch_filename
                        A signatures .h5 file created by gambit signatures
  database_patch_filename
                        An SQLite database

options:
  -h, --help            show this help message and exit
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: patched_database.gs)
  --database_output_filename DATABASE_OUTPUT_FILENAME, -d DATABASE_OUTPUT_FILENAME
                        Output filename for database (default: patched_database.gdb)
  --signatures_main_removed_filename SIGNATURES_MAIN_REMOVED_FILENAME
                        Output filename for genome signatures with patched genomes removed (default: main_database_removed.gs)
  --database_main_removed_filename DATABASE_MAIN_REMOVED_FILENAME
                        Output filename for database with patched genomes removed (default: main_database_removed.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```


## gambitdb-compress
```
usage: gambitdb-compress [options]

Compresses a GAMBIT database by filtering out genomes based on specified criteria.

positional arguments:
  signatures_filename   A signatures .h5 file created by gambit signatures
  database_filename     An sqlite database file created by gambit

options:
  -h, --help            show this help message and exit
  --min_species_genomes MIN_SPECIES_GENOMES, -g MIN_SPECIES_GENOMES
                        Minimum number of genomes in a species to consider, ignore the species below this (default: 10)
  --max_species_genomes MAX_SPECIES_GENOMES, -t MAX_SPECIES_GENOMES
                        Max number of genomes in a species to consider, ignore all others above this (default: 100)
  --core_proportion CORE_PROPORTION, -c CORE_PROPORTION
                        Proportion of genomes a kmer must be in for a species to be considered core (default: 1)
  --cpus CPUS, -p CPUS  Number of cpus to use (default: 1)
  --num_genomes_per_species NUM_GENOMES_PER_SPECIES, -r NUM_GENOMES_PER_SPECIES
                        Number of genomes to keep for a species (0 means keep all) (default: 1)
  --keep_species_under_minima
                        Keep species under minima rather than removing them (default: False)
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: filtered_database.gs)
  --database_output_filename DATABASE_OUTPUT_FILENAME, -d DATABASE_OUTPUT_FILENAME
                        Output filename for genome database (default: filtered_database.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-remove-genome-signatures
```
usage: gambitdb-remove-genome-signatures [options]

Given a Gambit signatures file, remove a list of genomes from it and return a new file.

positional arguments:
  signatures_filename   A signatures .h5 file created by gambit signatures
  genomes_to_remove_filename
                        One accession per line in a file

options:
  -h, --help            show this help message and exit
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: filtered_database.gs)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-repair-db
```
usage: gambitdb-repair-db [options]

Given a GAMBIT database, repair it

positional arguments:
  database_main_filename
                        An SQLite database

options:
  -h, --help            show this help message and exit
  --database_output_filename DATABASE_OUTPUT_FILENAME, -d DATABASE_OUTPUT_FILENAME
                        Output filename for database (default: fixed_database.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-rebuild-signatures
```
usage: gambitdb-rebuild-signatures [options]

Given a signatures file, rebuild the database.

positional arguments:
  signatures_filename   A signatures .h5 file created by gambit signatures

options:
  -h, --help            show this help message and exit
  --db_key DB_KEY       Unique key for database, no spaces (default: organisation/database)
  --db_version DB_VERSION
                        Unique version, x.y.z (default: 1.0.0)
  --db_author DB_AUTHOR
                        Name of person who created the database (default: Jane Doe)
  --db_date DB_DATE     Date database was created as YYYY-MM-DD (default: 2022-12-31)
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: database.gs)
  --database_output_filename DATABASE_OUTPUT_FILENAME, -g DATABASE_OUTPUT_FILENAME
                        Output filename for core database (default: database.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-iterative-build
```
usage: gambitdb-iterative-build [options]

Iteratively add species to a database

positional arguments:
  signatures_main_filename
                        A signatures .h5 file created by gambit signatures
  database_main_filename
                        An SQLite database
  gtdb_metadata_spreadsheet
                        GTDB database

options:
  -h, --help            show this help message and exit
  --species_added SPECIES_ADDED
                        file containing a list of species added (default: species_added)
  --min_genomes MIN_GENOMES
                        minimum genomes a species must have (default: 2)
  --species_to_ignore SPECIES_TO_IGNORE
                        file containing a list of species to ignore (default: None)
  --rank RANK, -r RANK  taxonomic rank (genus/species) (default: species)
  --cpus CPUS, -p CPUS  Number of cpus to use (default: 1)
```

# Scripts provided for testing
These are scripts which allow you to access functionality deep within the GAMBITdb software and are mostly for advanced usage, debugging or restarting a failed database run partway through.

## gambitdb-curate
```
usage: gambitdb-curate [options]

Given a species taxon file, and a genome file with metadata, curate the data and produce new files

positional arguments:
  species_taxon_filename
                        CSV containing species taxonomy and diameters, ngenomes - output from gambitdb-diameters
  genome_assembly_metadata
                        A CSV containing the assembly file and path, and the species taxon ID - output from gambitdb-gtdb
  assembly_directory    A directory containing assemblies in FASTA format
  pairwise_distances_filename
                        A pairwise distance file between each assembly

options:
  -h, --help            show this help message and exit
  --species_to_remove SPECIES_TO_REMOVE, -s SPECIES_TO_REMOVE
                        Optional file containing a list of species to remove (1 per line) (default: None)
  --accessions_to_remove ACCESSIONS_TO_REMOVE, -a ACCESSIONS_TO_REMOVE
                        Optional file containing a list of accession numbers to remove (1 per line) (default: None)
  --species_taxon_output_filename SPECIES_TAXON_OUTPUT_FILENAME, -b SPECIES_TAXON_OUTPUT_FILENAME
                        Output filename for genome signatures (default: species_taxon_curated.csv)
  --genome_assembly_metadata_output_filename GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME, -g GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME
                        Output filename for core database (default: genome_assembly_metadata_curated.csv)
  --accession_removed_output_filename ACCESSION_REMOVED_OUTPUT_FILENAME, -c ACCESSION_REMOVED_OUTPUT_FILENAME
                        Output filename for a list of accessions removed (default: accessions_removed.csv)
  --species_removed_output_filename SPECIES_REMOVED_OUTPUT_FILENAME, -d SPECIES_REMOVED_OUTPUT_FILENAME
                        Output filename for a list of species removed (default: species_removed.csv)
  --minimum_ngenomes MINIMUM_NGENOMES, -n MINIMUM_NGENOMES
                        Minimum number of genomes for a species to be included (default: 2)
  --small_cluster_ngenomes SMALL_CLUSTER_NGENOMES
                        Minimum number of genomes for a species to be included in a small cluster, along with --small_cluster_diameter (default: 4)
  --small_cluster_diameter SMALL_CLUSTER_DIAMETER
                        Maximum diameter of a species to be included in a small cluster along with --small_cluster_ngenomes (default: 0.7)
  --maximum_diameter MAXIMUM_DIAMETER
                        The maximum diameter to allow before attempting to split a species into subspecies (default: 0.7)
  --minimum_cluster_size MINIMUM_CLUSTER_SIZE
                        After splitting a species into subspecies, this is the minimum number of genomes which must be present, otherwise the genome is removed. (default: 2)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-diameters
```
usage: gambitdb-diameters [options]

Given files containing assembly metadata, pairwise distances and species taxon information output a new species file with diameters, and a min-inter file

positional arguments:
  genome_assembly_metadata
                        A CSV containing the assembly file and path, and the species taxon ID
  pairwise_distances_filename
                        A pairwise distance file between each assembly
  species_taxon_filename
                        A CSV containing species taxon IDs

options:
  -h, --help            show this help message and exit
  --species_taxon_output_filename SPECIES_TAXON_OUTPUT_FILENAME, -s SPECIES_TAXON_OUTPUT_FILENAME
                        Output filename for the modified species taxon IDs plus diameters (default: species_data_diameters.csv)
  --min_inter_output_filename MIN_INTER_OUTPUT_FILENAME, -i MIN_INTER_OUTPUT_FILENAME
                        Output filename for min inter values (default: min_inter.csv)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-gtdb-testset
```
usage: gambitdb-gtdb-testset [options]

Check a gambit database built using GTDB against genomes which were not used in the training set. Produces a list of genomes to use and their predicted species (GTDB)

positional arguments:
  species_taxon_file    Species taxon file produced by gambitdb, e.g. species_taxon.csv
  assembly_metadata_file
                        Assembly metadata file produced by gambitdb, e.g. genome_assembly_metadata.csv
  gtdb_metadata_file    GTDB spreadsheet

options:
  -h, --help            show this help message and exit
  --output_assembly_list_filename OUTPUT_ASSEMBLY_LIST_FILENAME, -a OUTPUT_ASSEMBLY_LIST_FILENAME
                        Output assemblies for download filename. These can be used with ncbi-genome-downloader (default: assemblies_for_download.txt)
  --output_assembly_to_species_filename OUTPUT_ASSEMBLY_TO_SPECIES_FILENAME, -b OUTPUT_ASSEMBLY_TO_SPECIES_FILENAME
                        Output assembly to species filename (default: species_to_assembly.csv)
  --max_genomes_per_species MAX_GENOMES_PER_SPECIES, -c MAX_GENOMES_PER_SPECIES
                        Max genomes per species (default: 5)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```
## gambitdb-pairwise-table
```
usage: gambitdb-pairwise-table [options]

Given a directory of assemblies in FASTA format, generate a table of pairwise distances

positional arguments:
  assembly_directory    A directory containing assemblies in FASTA format

options:
  -h, --help            show this help message and exit
  --accessions_to_remove ACCESSIONS_TO_REMOVE, -a ACCESSIONS_TO_REMOVE
                        Optional file containing a list of accession numbers to remove/ignore if found in assembly directory (1 per line) (default: None)
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: signatures.h5)
  --distance_table_output_filename DISTANCE_TABLE_OUTPUT_FILENAME, -p DISTANCE_TABLE_OUTPUT_FILENAME
                        Output filename for pairwise distance table (default: pw-dists.csv)
  --kmer KMER, -k KMER  Length of the k-mer to use (default: 11)
  --kmer_prefix KMER_PREFIX, -f KMER_PREFIX
                        Kmer prefix (default: ATGAC)
  --cpus CPUS, -c CPUS  Number of cpus to use (default: 1)
  --debug               Turn on debugging (default: False)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-merge-signatures
```
usage: gambitdb-merge-signatures [options]

Given two Gambit signatures files, merge them and return a new file.

positional arguments:
  signatures_main_filename
                        A signatures .h5 file created by gambit signatures
  signatures_patch_filename
                        A patch signatures .h5 file created by gambit signatures

options:
  -h, --help            show this help message and exit
  --signatures_output_filename SIGNATURES_OUTPUT_FILENAME, -s SIGNATURES_OUTPUT_FILENAME
                        Output filename for genome signatures (default: merged_database.gs)
  --verbose, -v         Turn on verbose output (default: False)
```


## Results

The results of the analysis will be output to '/path/to/working_directory/final'

## Contributing

Contributions to this project are welcome. To contribute, please fork the repository and submit a pull request.

## License

This project is licensed under the GNU GPL 3 License - see the [LICENSE](LICENSE) file for details.