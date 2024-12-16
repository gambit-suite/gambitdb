# GambitDB
This is a tool to generate a Gambit database. Its primary input is a spreadsheet from GTDB such as this [release](https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_metadata_r214.tar.gz).


# Installation
This software can be installed using pip:
```
pip install .
```
or
```
pip install git+https://github.com/gambit-suite/gambitdb.git
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
### Description
There is a single shell script for creating a database from a GTDB spreadsheet (one for bacteria, one for archaea). This script will download the data, create a GAMBIT database and signatures files, then check the recall against the downloaded files for QC purposes. It is the easiest way to create a database.

### Usage
To create a database in one step, execute the `run_gambitdb_bacteria.sh` script:

```bash
./run_gambitdb_bacteria.sh /path/to/working_directory /path/to/gtdb_spreadsheet.tsv num_cores
```
This script will create a GAMBIT database from a GTDB spreadsheet. It parses the spreadsheet, downloads data with ncbi-genome-download and outputs a GAMBIT database.

## gambitdb-gtdb

### Description
This script will parse a GTDB spreadsheet (see https://data.gtdb.ecogenomic.org/releases) and output a list of accessions to download, a species taxon file and a genome metadata file. It is the first step in creating a database.
The script provides several options for customization, including the ability to set a maximum number of contigs, include derived samples from metagenomes, environment, single cell, include novel species, set a minimum number of genomes in a species, and specify output filenames for the taxonomy, genome metadata, and genome accessions for download.

### Usage
To run this script, use the following command:
```
gambitdb-gtdb /path/to/gtdb_spreadsheet.tsv
```

The parameters for the script are:
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
### Description
The `gambitdb` script is used to generate a Gambit database. It requires a directory containing assemblies in FASTA format and a CSV file containing the assembly file and path, and the species taxon ID. Optionally, it can also take a CSV containing species taxonomy. The script provides several options for customization, including the ability to specify species and accessions to remove, output directory, output filenames, k-mer length, k-mer prefix, minimum number of genomes for a species to be included, number of CPUs to use, and parameters for including a species in a small cluster.

### Usage
To run this script, use the following command:
```bash
gambitdb [options] <assembly_directory> <genome_assembly_metadata> <species_taxon_filename> 
```

The parameters for the script are:
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
### Description
The `gambitdb-create` script is used to generate a GAMBIT database. It requires preprocessed input files including a CSV containing the assembly file and path, and the species taxon ID, a CSV containing species taxonomy, and a signatures .h5 file created by gambit signatures.

### Usage
To run this script, use the following command:

```
gambitdb-create [options] genome_assembly_metadata species_taxon_filename signatures_filename 
```

The parameters for the script are:
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

### Description
The `gambitdb-apply-patch` script is used to apply a patch to an existing GAMBIT database. It requires a signatures .h5 file created by gambit signatures, an SQLite database, a signatures .h5 file created by gambit signatures, and an SQLite database. The script provides several options for customization, including the ability to specify output filenames and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-apply-patch [options] <signatures_main_filename> <database_main_filename> <signatures_patch_filename> <database_patch_filename>
```

The parameters for the script are:
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
### Description
The `gambitdb-compress` script is used to compress a GAMBIT database by filtering out genomes based on specified criteria. It requires a signatures .h5 file created by gambit signatures and an SQLite database. The script provides several options for customization, including the ability to specify output filenames, the minimum number of genomes in a species to consider, the maximum number of genomes in a species to consider, the proportion of genomes a kmer must be in for a species to be considered core, the number of cpus to use, the number of genomes to keep for a species, and whether to keep species under minima rather than removing them.

### Usage
To run this script, use the following command:
```
gambitdb-compress [options] <signatures_filename> <database_filename>
```

The parameters for the script are:
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
### Description
The `gambitdb-remove-genome-signatures` script is used to remove a list of genomes from a GAMBIT database. It requires a signatures .h5 file created by gambit signatures and a list of genomes to remove. The script provides several options for customization, including the ability to specify output filenames and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-remove-genome-signatures [options] <signatures_filename> <genomes_to_remove_filename>
```

The parameters for the script are:
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
### Description
The `gambitdb-repair-db` script is used to repair a GAMBIT database. It requires an SQLite database. The script provides several options for customization, including the ability to specify output filenames and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-repair-db [options] <database_main_filename>
```

The parameters for the script are:
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
### Description
The `gambitdb-rebuild-signatures` script is used to rebuild a GAMBIT database. It requires a signatures .h5 file created by gambit signatures. The script provides several options for customization, including the ability to specify output filenames, database key, database version, database author, database date, and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-rebuild-signatures [options] <signatures_filename>
```

The parameters for the script are:
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
### Description
The `gambitdb-iterative-build` script is used to iteratively build a GAMBIT database. It requires a signatures .h5 file created by gambit signatures, an SQLite database, and a GTDB spreadsheet. The script provides several options for customization, including the ability to specify output filenames, the minimum number of genomes a species must have, a list of species to ignore, the taxonomic rank to use, and the number of cpus to use.

### Usage
To run this script, use the following command:
```
gambitdb-iterative-build [options] <signatures_main_filename> <database_main_filename> <gtdb_metadata_spreadsheet>
```

The parameters for the script are:
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
### Description
The `gambitdb-curate` script is used to curate a GAMBIT database. It requires a species taxon file, a genome metadata file, a directory containing assemblies in FASTA format, and a pairwise distance file between each assembly. The script provides several options for customization, including the ability to specify species and accessions to remove, output filenames, the minimum number of genomes for a species to be included, the number of cpus to use, the number of genomes for a species to be included in a small cluster, the maximum diameter of a species to be included in a small cluster, the maximum diameter to allow before attempting to split a species into subspecies, and the minimum number of genomes which must be present after splitting a species into subspecies.

### Usage
To run this script, use the following command:
```
gambitdb-curate [options] <species_taxon_filename> <genome_assembly_metadata> <assembly_directory> <pairwise_distances_filename>
```

The parameters for the script are:
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

### Description
The `gambitdb-diameters` script is used to calculate the diameters of species in a GAMBIT database. It requires a species taxon file, a pairwise distance file between each assembly, and a CSV containing species taxon IDs. The script provides several options for customization, including the ability to specify output filenames and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-diameters [options] <genome_assembly_metadata> <pairwise_distances_filename> <species_taxon_filename>
```

The parameters for the script are:
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

### Description
The `gambitdb-gtdb-testset` script is used to check a GAMBIT database built using GTDB against genomes which were not used in the training set. It requires a species taxon file, a genome metadata file, and a GTDB spreadsheet. The script provides several options for customization, including the ability to specify output filenames, the maximum number of genomes per species, and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-gtdb-testset [options] <species_taxon_file> <assembly_metadata_file> <gtdb_metadata_file>
```

The parameters for the script are:
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
### Description
The `gambitdb-pairwise-table` script is used to generate a table of pairwise distances between assemblies. It requires a directory containing assemblies in FASTA format. The script provides several options for customization, including the ability to specify output filenames, the k-mer length, the k-mer prefix, the number of cpus to use, and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-pairwise-table [options] <assembly_directory>
```

The parameters for the script are:
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
### Description
The `gambitdb-merge-signatures` script is used to merge two GAMBIT signatures files. It requires two signatures .h5 files created by gambit signatures. The script provides several options for customization, including the ability to specify output filenames and turn on verbose output.

### Usage
To run this script, use the following command:
```
gambitdb-merge-signatures [options] <signatures_main_filename> <signatures_patch_filename>
```

The parameters for the script are:
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

## gambitdb-update-taxa-report
### Description
The `gambitdb-update-taxa-report` script updates report flags and taxonomic rankgins in a GAMBIT database. The script specifically handles subspeces designations and report flag settings for different taxonomic rankings. For our use cases, we set subspecies report to 0, and species & genus to 1. 

### Usage
To run this script, use the following command:
```
gambitdb-update-flags [options] <database_filename>
```

The parameters for the script are:
```
usage: update-report-flags [options]

Update report flags in GAMBIT database

positional arguments:
  database_filename     An SQLite database

options:
  -h, --help            show this help message and exit
  --database_output_filename DATABASE_OUTPUT_FILENAME, -d DATABASE_OUTPUT_FILENAME
                        Output filename for database (default: updated_database.gdb)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-merge-duplicates
### Description
The `gambitdb-merge-duplicates` script identifies and merges duplicate taxa entries in a gambit db. The script will preserve entries with assigned genomes and performs a safe merge operation by creating a copy of the database and will generate a summary of all merges. 

To run this script:
```
gambitdb-merge-duplicates [options] <database_filename>
```

The parameters for the script are:
```
usage: merge-duplicate-taxa [options]

Merge duplicate taxa entries in GAMBIT database, keeping entries with assigned genomes

positional arguments:
  database_filename     A GAMBIT SQLite database

options:
  -h, --help            show this help message and exit
  --database_output_filename DATABASE_OUTPUT_FILENAME, -d DATABASE_OUTPUT_FILENAME
                        Output filename for database (default: merged_duplicates_database.gdb)
  --summary_filename SUMMARY_FILENAME, -s SUMMARY_FILENAME
                        Output filename for merge summary (default: merged_duplicates_summary.csv)
  --verbose, -v         Turn on verbose output (default: False)
```

## gambitdb-fungi
### Description
The `gambitdb-fungi` script processes NCBI RefSeq fungal genome data to create a GAMBIT database. It takes a RefSeq assembly summary file as input, filters genomes based on quality criteria, downloads assemblies, and generates the necessary metadata files. The script provides extensive options for customization, including the ability to set contig limits, minimum genomes per species, taxonomic grouping levels, and handling of atypical or metagenome-derived genomes. Default max contigs set to 100,000 to not filter out refseq matches.

### Usage
To run this script, use the following command:
```
gambitdb-fungi [options] <fungi_metadata_spreadsheet>
```

```
usage: gambitdb-fungi [options]

Given a Fungi RefSeq metadata spreadsheet, output a list of accessions to download, a species taxonid file and a genome metadata file

positional arguments:
  fungi_metadata_spreadsheet    Fungi metadata file

options:
  -h, --help            show this help message and exit
  --max_contigs MAX_CONTIGS, -d MAX_CONTIGS
                        Maximum number of contigs (default: 100000)
  --minimum_genomes_per_species MINIMUM_GENOMES_PER_SPECIES, -i MINIMUM_GENOMES_PER_SPECIES
                        Minimum number of genomes in a species (default: 2)
  --exclude_atypical    Exclude atypical genomes when searching NCBI (default: True)
  --is_metagenome_derived IS_METAGENOME_DERIVED
                        Check if the genome metagenome derived (MAGs) (default: metagenome_derived_exclude)
  --parent_taxonomy PARENT_TAXONOMY
                        Parent taxonomy to search for related to input taxa (choices: genus, family, order, class, phylum, kingdom) (default: genus)
  --genome_assembly_metadata_output_filename GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME, -g GENOME_ASSEMBLY_METADATA_OUTPUT_FILENAME
                        Genome metadata (default: assembly_metadata.csv)
  --species_taxon_output_filename SPECIES_TAXON_OUTPUT_FILENAME, -a SPECIES_TAXON_OUTPUT_FILENAME
                        Species data (default: species_taxon.csv)
  --filtered_out_genomes_filename FILTERED_OUT_GENOMES_FILENAME
                        Filtered out genomes (default: filtered_out_genomes.csv)
  --output_fasta_directory OUTPUT_FASTA_DIRECTORY, -o OUTPUT_FASTA_DIRECTORY
                        Output directory for fastas (default: .)
  --debug              Turn on debugging (default: False)
  --verbose, -v        Turn on verbose output (default: False)
```

## gambitdb-fungi-analyze
### Description
This script analyzes species distances and identifies overlaps in genome data AFTER the gambitdb run. This tool is useful for understanding where overlaps exists even after the initial gambitdb end-end script to see where there are still fundamental red-flags with the data. 

### Usage 
To run this script use the following command:
```
gambitdb-fungi-analyze [options] <genomes_path> <species_path> <distances_path>
```

The parameters for the script are:

```
usage: gambitdb-fungi-analyze [options]

Analyze species distances and identify overlaps

positional arguments:
  genomes_path          Path to curated genomes CSV
  species_path          Path to curated taxa CSV  
  distances_path        Path to pairwise distances CSV

options:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory (default: species_analysis_output)
  --min-distance-threshold MIN_DISTANCE_THRESHOLD, -t MIN_DISTANCE_THRESHOLD
                        Minimum distance threshold for clustering (default: 0.7)
  --debug               Enable debug logging (default: False)
  --verbose, -v         Enable verbose output (default: False)
```

## gambitdb-fungi-fix-genera
### Description
The `gambitdb-fungi-fix-genera` script fixes species diameters in a gambit db by examning and updating genera with zero-vaue distance thresholds that can be caused by subspeciation and genera introduction. This script will create a working copy with updated genus thresholds.

### Usage
To run this script, use the following command:

```
gambitdb-fungi-fix-genera [options] <source_database_filename>
```

The parameters for the script are:
```
usage: gambitdb-fungi-fix-genera [options]

Apply gambit genus diameters from one database to another

positional arguments:
  source_database_filename     Source gambit database file containing genus diameters

options:
  -h, --help            show this help message and exit
  --output_filename OUTPUT_FILENAME, -d OUTPUT_FILENAME
                        Output filename for database with genus diameters applied (default: fixed_database.gdb)
```
## Contributing

Contributions to this project are welcome. To contribute, please fork the repository and submit a pull request.

## License

This project is licensed under the GNU GPL 3 License - see the [LICENSE](LICENSE) file for details.
