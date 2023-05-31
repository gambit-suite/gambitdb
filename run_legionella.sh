cd ~/code/gambitdb/
 BASEDIR=~/data/gambitdb/Legionella_20230531_small
 ./scripts/gambitdb-gtdb  -s $BASEDIR/species_taxa.csv -g $BASEDIR/assembly_metadata.csv -a $BASEDIR/accessions_to_download.csv $BASEDIR/Legionella_r214.tsv  
cd $BASEDIR/
mkdir fasta
mkdir intermediate_files

# This software looks up different databases so we have to try for refseq and genbank
ncbi-genome-download -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s genbank -p 4 bacteria 
ncbi-genome-download -A accessions_to_download.csv -F fasta  -o fasta --flat-output -s refseq -p 4 bacteria 

# These files have a strain ID embedded so we need to remove them
cd fasta
for file in *_genomic.fna.gz; do
    mv "$file" "$(echo "$file" | awk -F '_' '{print $1"_"$2".fna.gz"}')"
done

cd ~/code/gambitdb/
mkdir $BASEDIR/intermediate_files
rm $BASEDIR/intermediate_files/database.gdb
rm $BASEDIR/intermediate_files/database.gs

rm -rf $BASEDIR/final
mkdir $BASEDIR/final
./scripts/gambitdb -d $BASEDIR/intermediate_files $BASEDIR/fasta $BASEDIR/assembly_metadata.csv $BASEDIR/species_taxa.csv
./scripts/gambit-create --database_output_filename $BASEDIR/final/database.gdb --signatures_output_filename $BASEDIR/final/database.gs $BASEDIR/intermediate_files/genome_assembly_metadata.csv $BASEDIR/intermediate_files/species_taxon.csv $BASEDIR/intermediate_files/database.gs
rm results.csv
gambit -d $BASEDIR/final query -o results.csv $BASEDIR/fasta/*.gz

./scripts/gambitdb-database-recall  $BASEDIR/assembly_metadata.csv results.csv

#### report in species taxon csv, default to 1 (0 means uncallable)
# species exists rather than being removed (if in species file) so has report

# assemblies removed from metadata pop up in the pairwise comparison file and cause errors

