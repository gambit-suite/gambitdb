#!/bin/bash

# This script is used to run the gambitdb-fungi script on a large input file
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1
lines_per_file=100
prefix="fungi_part_"

echo "Splitting input file..."

# Get the headers (first two lines starting with # for refseq)
header1=$(grep "^##" "$input_file")
header2=$(grep "^#[^#]" "$input_file")

total_lines=$(grep -v "^#" "$input_file" | wc -l)

grep -v "^#" "$input_file" > temp_content.txt

split_counter=1
for i in $(seq 0 $lines_per_file $total_lines); do
    output_file="${prefix}${split_counter}.txt"
    
    # Add headers to each file
    echo "$header1" > "$output_file"
    echo "$header2" >> "$output_file"
    
    # Add the next batch of lines
    head -n $lines_per_file <(tail -n +$((i + 1)) temp_content.txt) >> "$output_file"
    
    split_counter=$((split_counter + 1))
done

# Clean up temporary file
rm temp_content.txt

total_parts=$((split_counter - 1))
echo "Split complete. Created $total_parts files."

# Process each split file
echo "Processing split files..."
for i in $(seq 1 $total_parts); do
    echo "Processing part $i of $total_parts..."
    python3 ./scripts/gambitdb-fungi "${prefix}${i}.txt" \
        -g "assembly_metadata_${i}.csv" \
        -a "species_taxon_${i}.csv" \
        --filtered_out_genomes_filename "filtered_out_genomes_${i}.csv" \
        -o fungi_data
done

echo "All parts processed. Concatenating results..."

# Concatenate assembly metadata files
echo "Concatenating assembly metadata files..."
head -n 1 "assembly_metadata_1.csv" > "final_assembly_metadata.csv"
for i in $(seq 1 $total_parts); do
    if [ -f "assembly_metadata_${i}.csv" ]; then
        tail -n +2 "assembly_metadata_${i}.csv" >> "final_assembly_metadata.csv"
    fi
done

# Concatenate species taxon files
echo "Concatenating species taxon files..."
head -n 1 "species_taxon_1.csv" > "final_species_taxon.csv"
for i in $(seq 1 $total_parts); do
    if [ -f "species_taxon_${i}.csv" ]; then
        tail -n +2 "species_taxon_${i}.csv" >> "final_species_taxon.csv"
    fi
done

# Concatenate filtered out genomes files
echo "Concatenating filtered out genomes files..."
head -n 1 "filtered_out_genomes_1.csv" > "final_filtered_out_genomes.csv"
for i in $(seq 1 $total_parts); do
    if [ -f "filtered_out_genomes_${i}.csv" ]; then
        tail -n +2 "filtered_out_genomes_${i}.csv" >> "final_filtered_out_genomes.csv"
    fi
done

# Clean up intermediate files
echo "Cleaning up intermediate files..."
for i in $(seq 1 $total_parts); do
    rm -f "${prefix}${i}.txt"
    rm -f "assembly_metadata_${i}.csv"
    rm -f "species_taxon_${i}.csv"
    rm -f "filtered_out_genomes_${i}.csv"
done


echo "Running signature gambitdb processing..."
./scripts/gambitdb -v --cpus 4 \
    -d fungi_db_files/ \
    ./fungi_data \
    final_assembly_metadata.csv \
    final_species_taxon.csv

echo "Running signature gambitdb-create..."
./scripts/gambitdb-create --database_output_filename ./final/database.gdb \ 
    --signatures_output_filename ./final/database.gs \ 
    fungi_db_files/genome_assembly_metadata.csv \ 
    fungi_db_files/species_taxon.csv \
    ./fungi_db_files/database.gs


echo "Process complete. Final files created:"
echo "- final_assembly_metadata.csv"
echo "- final_species_taxon.csv"
echo "- final_filtered_out_genomes.csv"
echo "- fungi_db_files/"
echo "- final/database.gdb"