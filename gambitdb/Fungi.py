import pandas as pd
import requests
import time
from typing import Dict, List, Set
from dataclasses import dataclass
import json
import logging

@dataclass
class GenomeMetadata:
    accession: str
    contig_count: int
    species_taxid: str
    organism_name: str
    assembly_level: str
    bioproject: str
    
class FungiParser:
    """
    Reads in a Fungi spreadsheet, parses each row and outputs a modified spreadsheet.
    """
    def __init__(self, 
                 fungi_metadata_spreadsheet: str, 
                 max_contigs: int, 
                 minimum_genomes_per_species: int,
                 genome_assembly_metadata_output_filename: str, 
                 accessions_output_filename: str,
                 verbose: bool = False,
                 debug: bool = False):
        """
        Initializes the FungiParser class.
        Args:
          fungi_metadata_spreadsheet (str): The path to the Fungi metadata spreadsheet.
          max_contigs (int): The maximum number of contigs to include a genome.
          minimum_genomes_per_species (int): The minimum number of genomes per species to include a species.
          species_taxon_output_filename (str): The path to the output file for species taxon information.
          genome_assembly_metadata_output_filename (str): The path to the output file for genome assembly metadata.
          accessions_output_filename (str): The path to the output file for accessions.
        """
        
        self.logger = logging.getLogger(__name__)
        self.debug = debug
        self.verbose = verbose
        
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)
        
        # Add a console handler to the logger
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG if self.verbose else logging.ERROR)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        self.max_contigs = max_contigs
        self.minimum_genomes_per_species = minimum_genomes_per_species
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.accessions_output_filename = accessions_output_filename
        
        #Refseq file is weird so getting the column names from the header, can change this later
        with open(fungi_metadata_spreadsheet) as f:
            
            f.readline()
            
            header = f.readline().strip()
            
            columns = header[1:].split('\t')
        
        self.df = pd.read_csv(
            fungi_metadata_spreadsheet, 
            sep='\t',
            names=columns,
            comment='#',
            skiprows=2,
            dtype={'species_taxid': str}
        )
        
        self.valid_species: Set[str] = set()
        self.valid_genomes: List[GenomeMetadata] = []
        self.species_genome_counts: Dict[str, int] = {}
        
    def process_genome_report(self, report: dict) -> GenomeMetadata:
        """
        Process a single genome report into GenomeMetadata from NCBI API
        """
        assembly_info = report.get('assembly_info', {})
        organism_info = report.get('organism', {})
        
        return GenomeMetadata(
            accession=report['accession'],
            contig_count=report['assembly_stats']['number_of_contigs'],
            species_taxid=str(organism_info.get('tax_id', '')),
            organism_name=organism_info.get('organism_name', ''),
            assembly_level=assembly_info.get('assembly_level', ''),
            bioproject=assembly_info.get('bioproject_accession', ''),
        )

    def get_genome_data_for_taxon(self, taxon_id: str) -> List[GenomeMetadata]:
        """
        Query NCBI Datasets API for genome data for a given taxon ID
        """
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon_id}/dataset_report"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            
            if 'reports' not in data:
                self.logger.warning(f"No genomes found for taxon {taxon_id}")
                return []
            
            reports = data['reports']
            total_genomes = len(reports)
            
            if total_genomes < self.minimum_genomes_per_species:
                self.logger.debug(f"Taxon {taxon_id} only has {total_genomes} genomes (minimum required: {self.minimum_genomes_per_species})")
                return []
            
            self.logger.debug(f"Found {total_genomes} genomes for taxon {taxon_id}, processing...")
            
            genome_list = []
            for report in reports:
                if (report.get('accession') and 
                    report.get('assembly_stats', {}).get('number_of_contigs') is not None):
                    
                    genome = self.process_genome_report(report)
                    genome_list.append(genome)
            
            return genome_list
            
        #Errors I would expect to see, not all fields always present in responst from NCBI API
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Error querying NCBI API for taxon {taxon_id}: {e}")
            return []
        except json.JSONDecodeError as e:
            self.logger.error(f"Error parsing JSON response for taxon {taxon_id}: {e}")
            return []
        except Exception as e:
            self.logger.error(f"Unexpected error processing taxon {taxon_id}: {e}")
            return []

    def process_fungi_data(self):
        """
        Process each unique species taxid through the NCBI API
        """
        
        #Grab unique species here
        unique_taxids = self.df['species_taxid'].unique()
        self.logger.debug(f"Processing {len(unique_taxids)} unique species taxids...")

        for taxon_id in unique_taxids:
            self.logger.debug(f"Processing taxon ID: {taxon_id}")
            
            # Get genome data from NCBI API
            genomes = self.get_genome_data_for_taxon(taxon_id)
            
            # Can probably handle this logic earlier
            if genomes:
                valid_genomes = [
                    genome for genome in genomes 
                    if genome.contig_count <= self.max_contigs
                ]
                
                if len(valid_genomes) >= self.minimum_genomes_per_species:
                    self.logger.debug(f"Taxon {taxon_id} has {len(valid_genomes)} valid genomes after contig filtering")
                    self.valid_species.add(taxon_id)
                    self.valid_genomes.extend(valid_genomes)
                    self.species_genome_counts[taxon_id] = len(valid_genomes)
                else:
                    self.logger.debug(f"Taxon {taxon_id} only has {len(valid_genomes)} valid genomes after contig filtering")
            
            # I added this sleep to avoid hitting the API too often
            time.sleep(0.1)
    
    def write_outputs(self):
        """
        Write the filtered data to output files
        """
        self.logger.debug("Writing output files for FungiParser")
        
        # Write genome assembly metadata
        genome_data = pd.DataFrame([
            {
                'accession': genome.accession,
                'contig_count': genome.contig_count,
                'species_taxid': genome.species_taxid,
                'organism_name': genome.organism_name,
                'assembly_level': genome.assembly_level,
                'bioproject': genome.bioproject
            }
            for genome in self.valid_genomes
        ])
        genome_data.to_csv(self.genome_assembly_metadata_output_filename, sep='\t', index=False)
        self.logger.debug(f"Wrote {len(self.valid_genomes)} genomes to {self.genome_assembly_metadata_output_filename}")
        
        # Write accessions output
        accessions_data = pd.DataFrame({
            'accession': [genome.accession for genome in self.valid_genomes]
        })
        accessions_data.to_csv(self.accessions_output_filename, sep='\t', index=False)
        self.logger.debug(f"Wrote {len(self.valid_genomes)} accessions to {self.accessions_output_filename}")
        
        
    def generate_spreadsheets(self):
      """Main method to process data and generate output files"""
      self.process_fungi_data()
      self.write_outputs()