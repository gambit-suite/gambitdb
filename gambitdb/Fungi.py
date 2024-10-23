import pandas as pd
import requests
import time
import sys
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
    parent_taxid: int
    rank: str
    
class FungiParser:
    """
    Reads in a Fungi spreadsheet, parses each row and outputs a modified spreadsheet.
    """
    def __init__(self, 
                 fungi_metadata_spreadsheet: str, 
                 max_contigs: int, 
                 minimum_genomes_per_species: int,
                 genome_assembly_metadata_output_filename: str, 
                 taxon_output_filename: str,
                 exclude_atypical: bool,
                 is_metagenome_derived: str,
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
          taxon_output_filename (str): The path to the output file for accessions.
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
        self.taxon_output_filename = taxon_output_filename
        
        #Datasets API filters
        self.exclude_atypical = str(exclude_atypical).lower() #necessary formattting for datasets API
        self.is_metagenome_derived = is_metagenome_derived
        
        
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
        
    def get_parent_taxon(self, taxon_id: str) -> str:
        """
        Get the parent taxon ID for a given taxon ID.
        """
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxon_id}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            parent_taxon = data['taxonomy_nodes'][0]['taxonomy']['lineage'][-1]
            rank = data['taxonomy_nodes'][0]['taxonomy']['rank']
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Error querying NCBI API for taxon {taxon_id}: {e}")
            sys.exit(1)
        except json.JSONDecodeError as e:
            self.logger.error(f"Error parsing JSON response for taxon {taxon_id}: {e}")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"Unexpected error processing taxon {taxon_id}: {e}")
            sys.exit(1)
        
        return parent_taxon, rank.lower()
        
    def process_genome_report(self, report: dict, parent_taxid: int, rank: str) -> GenomeMetadata:
        """
        Process a single genome report into GenomeMetadata from NCBI API
        """

        organism_info = report.get('organism', {})
        
        return GenomeMetadata(
            accession=report['accession'],
            contig_count=report['assembly_stats']['number_of_contigs'],
            species_taxid=str(organism_info.get('tax_id', '')),
            organism_name=organism_info.get('organism_name', ''),
            parent_taxid=parent_taxid,
            rank=rank
        )

    def get_genome_data_for_taxon(self, taxon_id: str, parent_taxid: int, rank: str) -> List[GenomeMetadata]:
        """
        Query NCBI Datasets API for genome data for a given taxon ID.
        Handles pagination for large result sets and performs validation of total genomes found.
        """
        
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon_id}/dataset_report"
        genome_list = []
        page_token = None
        
        try:
            
            params = {
            'filters.exclude_atypical': self.exclude_atypical,
            'filters.is_metagenome_derived': self.is_metagenome_derived
            }
            
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            
            if 'reports' not in data:
                self.logger.warning(f"No genomes found for taxon {taxon_id}")
                return []
            
            total_count = data.get('total_count', len(data['reports']))
            if total_count < self.minimum_genomes_per_species:
                self.logger.debug(f"Taxon {taxon_id} has {total_count} total genomes (minimum required: {self.minimum_genomes_per_species})")
                return []
                
            self.logger.debug(f"Found {total_count} total genomes for taxon {taxon_id}, processing...")
            
            reports = data['reports']
            for report in reports:
                if (report.get('accession') and 
                    report.get('assembly_stats', {}).get('number_of_contigs') is not None):
                    genome = self.process_genome_report(report, parent_taxid, rank)
                    genome_list.append(genome)
                    
            while True:
                page_token = data.get('next_page_token')
                if not page_token:
                    break
                    
                self.logger.debug(f"Fetching next page for taxon {taxon_id}")
                
                paginated_url = f"{url}?page_token={page_token}"
                response = requests.get(paginated_url)
                response.raise_for_status()
                data = response.json()
                
                reports = data['reports']
                for report in reports:
                    if (report.get('accession') and 
                        report.get('assembly_stats', {}).get('number_of_contigs') is not None):
                        genome = self.process_genome_report(report, parent_taxid, rank)
                        genome_list.append(genome)
                        
        except requests.exceptions.RequestException as e:
            self.logger.error(f"Error querying NCBI API for taxon {taxon_id}: {e}")
            return genome_list if len(genome_list) >= self.minimum_genomes_per_species else []
        except json.JSONDecodeError as e:
            self.logger.error(f"Error parsing JSON response for taxon {taxon_id}: {e}")
            return genome_list if len(genome_list) >= self.minimum_genomes_per_species else []
        except Exception as e:
            self.logger.error(f"Unexpected error processing taxon {taxon_id}: {e}")
            return genome_list if len(genome_list) >= self.minimum_genomes_per_species else []
            
        return genome_list

    def process_fungi_data(self):
        """
        Process each unique species taxid through the NCBI API
        """
        
        unique_taxids = self.df['species_taxid'].unique()
        self.logger.debug(f"Processing {len(unique_taxids)} unique species taxids...")

        for taxon_id in unique_taxids:
            self.logger.debug(f"Processing taxon ID: {taxon_id}")
            
            # Get parent taxon ID & rank
            parent_taxon,rank = self.get_parent_taxon(taxon_id)
            
            # Get genome data from NCBI API
            genomes = self.get_genome_data_for_taxon(taxon_id, parent_taxon, rank)
            
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
        
        #uuid,assembly_filename,species_taxid,assembly_accession
        genome_data = pd.DataFrame([
            {
                'uuid': genome.accession,
                'assembly_filename': f'{genome.accession}_.fasta',
                'species_taxid': genome.species_taxid,
                'accession': genome.accession
            }
            for genome in self.valid_genomes
        ])
        genome_data.to_csv(self.genome_assembly_metadata_output_filename, index=False)
        self.logger.debug(f"Wrote {len(self.valid_genomes)} genomes to {self.genome_assembly_metadata_output_filename}")
        
        #species_taxid,name,rank,parent_taxid,ncbi_taxid,gambit_taxid
        taxon_data = pd.DataFrame([
            {
                'species_taxid': genome.species_taxid,
                'name': genome.organism_name,
                'rank': genome.rank,
                'parent_taxid': genome.parent_taxid,
                'ncbi_taxid': genome.species_taxid,
                'gambit_taxid': genome.species_taxid
            }
            for genome in self.valid_genomes
            ])
        taxon_data.to_csv(self.taxon_output_filename, index=False)
        self.logger.debug(f"Wrote {len(self.valid_genomes)} accessions to {self.taxon_output_filename}")
        
        
    def generate_spreadsheets(self):
      """Main method to process data and generate output files"""
      self.process_fungi_data()
      self.write_outputs()