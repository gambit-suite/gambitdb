#!/bin/bash
# run_gambitdb.sh ~/path/to/workin/dir gtdb_file.tsv 1
# The file should be in the directory. Nothing else should be there as it creates/deletes folders.
BASEDIR=$1
INPUTFILE=$2
CORES=$3

cd ~/code/gambitdb/
rm output_logfile.txt
 ./scripts/gambitdb-gtdb  -s $BASEDIR/species_taxa.csv -g $BASEDIR/assembly_metadata.csv -a $BASEDIR/accessions_to_download.csv $BASEDIR/$INPUTFILE 
cd $BASEDIR/
mkdir $BASEDIR/fasta
rm -rf $BASEDIR/intermediate_files
mkdir $BASEDIR/intermediate_files
rm -rf $BASEDIR/final
mkdir $BASEDIR/final

# This software looks up different databases so we have to try for refseq and genbank
ncbi-genome-download -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s genbank -p $CORES bacteria 
ncbi-genome-download -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s refseq -p $CORES bacteria 

# These files have a strain ID embedded so we need to remove them
cd fasta
for file in *_genomic.fna.gz; do
    mv "$file" "$(echo "$file" | awk -F '_' '{print $1"_"$2".fna.gz"}')"
done

cd ~/code/gambitdb/
./scripts/gambitdb -v --cpus $CORES -d $BASEDIR/intermediate_files $BASEDIR/fasta $BASEDIR/assembly_metadata.csv $BASEDIR/species_taxa.csv
./scripts/gambit-create --database_output_filename $BASEDIR/final/database.gdb --signatures_output_filename $BASEDIR/final/database.gs $BASEDIR/intermediate_files/genome_assembly_metadata.csv $BASEDIR/intermediate_files/species_taxon.csv $BASEDIR/intermediate_files/database.gs
rm results.csv
gambit -d $BASEDIR/final query -o $BASEDIR/results.csv $BASEDIR/fasta/*.gz

./scripts/gambitdb-database-recall -o $BASEDIR/recall_results.txt $BASEDIR/assembly_metadata.csv $BASEDIR/results.csv

