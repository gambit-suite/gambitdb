#!/bin/bash
# run_gambitdb_fungi.sh ~/path/to/workin/dir 1
# The file should be in the directory. Nothing else should be there as it creates/deletes folders.
BASEDIR=$1
INPUTFILE=$2
CORES=$2

# Set variables
DATE=$(date +%Y-%m-%d)

# Step 1 - Download accessions from ncbi with metadata using datasets
# datasets summary genome taxon 4751 --as-json-lines --assembly-source 'GenBank' --exclude-atypical | dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,assmstats-number-of-contigs,assmstats-gc-percent,assmstats-contig-n50,organism-name,organism-tax-id | tee ${BASEDIR}/${DATE}-ncbi-fungi.tsv


./scripts/gambitdb-fungi $INPUTFILE -o fungi_data

./scripts/gambitdb -v --cpus 4 -d intermediate_files/ \
 ~/gambitdbdata fungi_outputs/assembly_metadata_with_species.csv \
  fungi_outputs/species_merged_taxon.csv
