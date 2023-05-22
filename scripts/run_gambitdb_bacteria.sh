#!/bin/bash
# run_gambitdb_bacteria.sh ~/path/to/workin/dir /path/to/gtdb_file.tsv 1
# The file should be in the directory. Nothing else should be there as it creates/deletes folders.
BASEDIR=$1
INPUTFILE=$2
CORES=$3

# Take the GTDB spreadsheet and filter the genomes to remove the ones that are lower quality
gambitdb-gtdb --checkm_completeness 95 --checkm_contamination 5 --include_novel_species --max_contigs 300 -s $BASEDIR/species_taxa.csv -g $BASEDIR/assembly_metadata.csv -a $BASEDIR/accessions_to_download.csv $INPUTFILE 

# Overwrite the existing base directory files and setup the folders
cd $BASEDIR/
mkdir $BASEDIR/fasta
rm -rf $BASEDIR/intermediate_files
mkdir $BASEDIR/intermediate_files
rm -rf $BASEDIR/final
mkdir $BASEDIR/final

# Download genomes from NCBI
# This software looks up different databases so we have to try for refseq and genbank
ncbi-genome-download -r 3 -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s genbank -p $CORES bacteria 
ncbi-genome-download -r 3 -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s refseq -p $CORES bacteria 

# These files have a strain ID embedded so we need to remove them
cd fasta
for file in *_genomic.fna.gz; do
    mv "$file" "$(echo "$file" | awk -F '_' '{print $1"_"$2".fna.gz"}')"
done

# Build the signatures and database files. This is the key part of the process
gambitdb -v --cpus $CORES -d $BASEDIR/intermediate_files $BASEDIR/fasta $BASEDIR/assembly_metadata.csv $BASEDIR/species_taxa.csv
gambitdb-create --database_output_filename $BASEDIR/final/database.gdb --signatures_output_filename $BASEDIR/final/database.gs $BASEDIR/intermediate_files/genome_assembly_metadata.csv $BASEDIR/intermediate_files/species_taxon.csv $BASEDIR/intermediate_files/database.gs

# Check the database works by running a query on the input genomes and looking at how they are 
# classified.
ls $BASEDIR/fasta/ | grep fna.gz > $BASEDIR/intermediate_files/assembly_filenames.txt
gambit -d $BASEDIR/final query -c $CORES  -l $BASEDIR/intermediate_files/assembly_filenames.txt --ldir $BASEDIR/fasta  -o $BASEDIR/results.csv
gambitdb-database-recall -o $BASEDIR/recall_results.txt $BASEDIR/intermediate_files/genome_assembly_metadata.csv $BASEDIR/results.csv
