import pandas as pd
from typing import Dict, List, Set, Optional
from dataclasses import dataclass
import json
import logging
from pathlib import Path
import zipfile
import shutil
import os
from gambitdb.NCBIDatasets import NCBIDatasetClient

@dataclass
class GenomeMetadata:
    accession: str
    assembly_name: str
    assembly_source: str
    contig_count: int
    species_taxid: str
    organism_name: str
    parent_taxid: int
    rank: str
    assembly_filename: str = None
    
@dataclass
class FilteredOutGenome:
    accession: str
    species_taxid: str
    organism_name: str
    reason: str
    
    
class FungiParser:
    """
    Reads in a Fungi spreadsheet, parses each row and outputs a modified spreadsheet.
    """
    def __init__(self, 
                 fungi_metadata_spreadsheet: str, 
                 max_contigs: int, 
                 minimum_genomes_per_species: int,
                 genome_assembly_metadata_output_filename: str,
                 output_fasta_directory: str, 
                 taxon_output_filename: str,
                 filtered_out_genomes_filename: str,
                 exclude_atypical: bool,
                 is_metagenome_derived: str,
                 parent_taxonomic_level: str,
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
        # Keep track of filtered out genomes
        self.filtered_genomes: List[FilteredOutGenome] = []
        
        self.max_contigs = max_contigs
        self.minimum_genomes_per_species = minimum_genomes_per_species
        self.genome_assembly_metadata_output_filename = genome_assembly_metadata_output_filename
        self.taxon_output_filename = taxon_output_filename
        self.filtered_genomes_output = filtered_out_genomes_filename
        self.output_fasta_directory = output_fasta_directory
        
        #Immediately set path
        self.genome_assembly_outdir = Path(os.getcwd()) / self.output_fasta_directory.lstrip('/')
        
        #Set up the NCBI Datasets API client
        self.api_client = NCBIDatasetClient(rate_limit=0.1)
        
        #Datasets API filters
        self.exclude_atypical = str(exclude_atypical).lower() #necessary formattting for datasets API
        self.is_metagenome_derived = is_metagenome_derived
        self.genome_report_params = {
            'filters.exclude_atypical': self.exclude_atypical,
            'filters.is_metagenome_derived': self.is_metagenome_derived
            }
        self.parent_taxonomic_level = str(parent_taxonomic_level).lower()
        self.prefered_assembly_database = "SOURCE_DATABASE_REFSEQ"
        
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
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxon_id}/dataset_report"
        
        try:
            taxon_data = self.api_client.get_json(url)
            parent_taxon = taxon_data['reports'][0]['taxonomy']['classification'][self.parent_taxonomic_level]['id']
            rank = taxon_data['reports'][0]['taxonomy']['rank'].lower()
        except Exception as e:
            self.logger.error(f"Unexpected error processing taxon {taxon_id}: {e}")
            return None, None
        
        return parent_taxon, rank.lower()
        
    def process_genome_report(self, report: dict, taxon_id: int, parent_taxid: int, rank: str) -> GenomeMetadata:
        """
        Process a single genome report into GenomeMetadata from NCBI API
        """

        organism_info = report.get('organism', {})
        
        return GenomeMetadata(
            accession=report['accession'],
            assembly_name=report['assembly_info']['assembly_name'],
            assembly_source=report['source_database'],
            contig_count=report['assembly_stats']['number_of_contigs'],
            species_taxid=str(taxon_id),
            organism_name=' '.join(organism_info.get('organism_name', '').split()[:2]),
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
            
            genome_data = self.api_client.get_json(url, params=self.genome_report_params)
            
            if 'reports' not in genome_data:
                self.logger.warning(f"No genomes found for taxon {taxon_id}")
                self.logger.error(f"No genomes found for taxon {taxon_id}")
                return []
            
            total_count = genome_data.get('total_count', len(genome_data['reports']))
            if total_count < self.minimum_genomes_per_species:
                self.logger.debug(f"Taxon {taxon_id} has {total_count} total genomes (minimum required: {self.minimum_genomes_per_species})")
                return []
                
            self.logger.debug(f"Found {total_count} total genomes for taxon {taxon_id}, processing...")
            
            reports = genome_data['reports']
            for report in reports:
                if (report.get('accession') and 
                    report.get('assembly_stats', {}).get('number_of_contigs') is not None):
                    genome = self.process_genome_report(report, taxon_id, parent_taxid, rank)
                    genome_list.append(genome)
                    
            while True:
                for report in genome_data.get('reports', []):
                    if (report.get('accession') and 
                        report.get('assembly_stats', {}).get('number_of_contigs') is not None):
                        if genome := self.process_genome_report(report, taxon_id, parent_taxid, rank):
                            genome_list.append(genome)
                
                # Handle pagination
                page_token = genome_data.get('next_page_token')
                if not page_token:
                    break
                    
                self.logger.debug(f"Fetching next page for taxon {taxon_id}")
                genome_data = self.api_client.get_json(url, params={**self.genome_report_params, 'page_token': page_token})
                        
        except Exception as e:
            self.logger.error(f"Unexpected error processing genome data for taxon {taxon_id}: {e}")
            return genome_list if len(genome_list) >= self.minimum_genomes_per_species else []
            
        return genome_list
    

    def process_fungi_data(self):
        """Process each unique species taxid through the NCBI API"""
        unique_taxids = self.df['species_taxid'].unique()
        self.logger.debug(f"Processing {len(unique_taxids)} unique species taxids...")

        for taxon_id in unique_taxids:
            self.logger.debug(f"Processing taxon ID: {taxon_id}")
            
            # Get parent taxon ID & rank
            parent_info = self.get_parent_taxon(taxon_id)
            if not parent_info:
                continue
            parent_taxon, rank = parent_info
            
            # Get genome data from NCBI API
            genomes = self.get_genome_data_for_taxon(taxon_id, parent_taxon, rank)
            if not genomes:
                continue
            
            # Filter genomes
            valid_genomes = []
            for genome in genomes:
                filter_reason = self._filter_genome(genome)
                if filter_reason:
                    self.filtered_genomes.append(
                        FilteredOutGenome(
                            genome.accession,
                            genome.species_taxid,
                            genome.organism_name,
                            filter_reason
                        )
                    )
                else:
                    valid_genomes.append(genome)
            
            # Apply assembly preferences
            original_count = len(valid_genomes)
            valid_genomes = self._handle_assembly_preferences(valid_genomes)
            
            # Track filtered assemblies
            if len(valid_genomes) < original_count:
                kept_accessions = {genome.accession for genome in valid_genomes}
                filtered_assemblies = [genome for genome in genomes if genome.accession not in kept_accessions]
                
                for genome in filtered_assemblies:
                    self.filtered_genomes.append(
                        FilteredOutGenome(
                            genome.accession,
                            genome.species_taxid,
                            genome.organism_name,
                            f"Non-preferred assembly source: {genome.assembly_source}"
                        )
                    )
            
            # Process valid genomes for species
            if len(valid_genomes) >= self.minimum_genomes_per_species:
                self.logger.debug(
                    f"Taxon {taxon_id} has {len(valid_genomes)} valid genomes after filtering"
                )
                self.valid_species.add(taxon_id)
                self.valid_genomes.extend(valid_genomes)
                self.species_genome_counts[taxon_id] = len(valid_genomes)
            else:
                self.logger.debug(
                    f"Taxon {taxon_id} only has {len(valid_genomes)} valid genomes after filtering"
                )
                for genome in valid_genomes:
                    self.filtered_genomes.append(
                        FilteredOutGenome(
                            genome.accession,
                            genome.species_taxid,
                            genome.organism_name,
                            f"Species has fewer than {self.minimum_genomes_per_species} valid genomes"
                        )
                    )
                    
    def _filter_genome(self, genome: GenomeMetadata) -> Optional[str]:
        """Return filter reason if genome should be filtered, otherwise return None"""
        
        if genome.parent_taxid is None:
            return "Missing parent taxon information"
        # Fix this to calculate fractional genome (genome size of representative/contig count)
        # 
        if genome.contig_count > self.max_contigs:
            return f"Contig count {genome.contig_count} exceeds maximum {self.max_contigs}"
        if 'sp.' in genome.organism_name:
            return "Contains 'sp.' in organism name"
        return None
                    
    def _handle_assembly_preferences(self, genomes):
        """
        Handle assembly source preferences by selecting preferred source when duplicates exist
        """
        
        assembly_groups = {}
        for genome in genomes:
            if genome.assembly_name not in assembly_groups:
                assembly_groups[genome.assembly_name] = []
            assembly_groups[genome.assembly_name].append(genome)
        
        preferred_genomes = []
        for assembly_name, group in assembly_groups.items():
            if len(group) == 1:
                preferred_genomes.append(group[0])
            else:
                refseq_assemblies = [genome for genome in group if genome.assembly_source == self.prefered_assembly_database]
                if refseq_assemblies:
                    preferred_genomes.append(refseq_assemblies[0])
                    self.logger.debug(f"Selected RefSeq assembly {refseq_assemblies[0].accession} "
                                    f"over GenBank alternatives for {assembly_name}")
                else:
                    preferred_genomes.append(group[0])
                    self.logger.debug(f"No RefSeq assembly available for {assembly_name}, "
                                    f"using {group[0].accession}")
        
        return preferred_genomes
                      
    def download_genomes_bulk(self):
        """Download genome assemblies in bulk groups by species_taxid"""
        
        self.logger.debug(f"Starting bulk genome downloads by species")
        self.genome_assembly_outdir.mkdir(exist_ok=True)
        
        # Group genomes by species, excluding already downloaded ones
        species_groups = self._group_genomes_for_download()
        if not species_groups:
            self.logger.debug("No new genomes to download")
            return
        
        # Process each species group
        for species_taxid, genomes in species_groups.items():
            self._download_species_group(species_taxid, genomes)
            
        successful_downloads = len([genome for genome in self.valid_genomes if hasattr(genome, 'assembly_filename')])
        self.logger.debug(f"Completed all bulk downloads. {successful_downloads} genomes successful")
    
    def _group_genomes_for_download(self) -> Dict[str, List[GenomeMetadata]]:
        """Group genomes by species_taxid, checking for existing files"""
        species_groups = {}
        for genome in self.valid_genomes:
            final_fasta = self.genome_assembly_outdir / f"{genome.accession}.fna"
            if final_fasta.exists():
                genome.assembly_filename = str(final_fasta.absolute())
            else:
                species_groups.setdefault(genome.species_taxid, []).append(genome)
        return species_groups
    
    def _download_species_group(self, species_taxid: str, genomes: List[GenomeMetadata]):
        """Download and process a group of genomes for a single species"""
        self.logger.debug(f"Downloading bulk group for species_taxid {species_taxid} ({len(genomes)} genomes)")
        
        request_body = {
            "accessions": [genome.accession for genome in genomes],
            "include_annotation_type": ["GENOME_FASTA"],
            "hydrated": "FULLY_HYDRATED",
            "include_tsv": False
        }
        
        temp_zip = self.genome_assembly_outdir / f"bulk_download_{species_taxid}_temp.zip"
        try:
            url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"
            self.api_client.download_file(
                url=url,
                output_path=temp_zip,
                method='POST',
                json=request_body
            )
            
            self._process_bulk_download(temp_zip, genomes)
            
        except Exception as e:
            self.logger.error(
                f"Error in bulk download for species {species_taxid}: {e}, "
                "removing genomes from valid list"
            )
            self._handle_failed_download(genomes)
        
        finally:
            self._cleanup_download_files(temp_zip)
    
    def _process_bulk_download(self, zip_path: Path, genomes: List[GenomeMetadata]):
        """Extract and process files from bulk download"""
        with zipfile.ZipFile(zip_path) as zip_ref:
            fasta_files = [f for f in zip_ref.namelist() if f.endswith('_genomic.fna')]
            
            for fasta_file in fasta_files:
                accession = fasta_file.split('/')[2]
                final_fasta = self.genome_assembly_outdir / f"{accession}.fna"
                
                try:
                    # Extract and move file
                    zip_ref.extract(fasta_file)
                    os.rename(
                        os.path.join('ncbi_dataset', 'data', accession, os.path.basename(fasta_file)),
                        final_fasta
                    )
                    
                    for genome in genomes:
                        if genome.accession == accession:
                            genome.assembly_filename = str(final_fasta.absolute())
                            break
                
                except (FileNotFoundError, OSError) as e:
                    self.logger.warning(f"Could not process file for {accession}: {e}")
                    self._remove_failed_genome(accession, genomes)
    
    def _handle_failed_download(self, genomes: List[GenomeMetadata]):
        """Handle failed downloads by removing genomes from valid list"""
        for genome in genomes:
            if genome in self.valid_genomes:
                self.valid_genomes.remove(genome)
                self.filtered_genomes.append(
                    FilteredOutGenome(
                        genome.accession,
                        genome.species_taxid,
                        genome.organism_name,
                        "Failed to download genome assembly"
                    )
                )
    
    def _remove_failed_genome(self, accession: str, genomes: List[GenomeMetadata]):
        """Remove a single failed genome from valid genomes list"""
        for genome in genomes:
            if genome.accession == accession and genome in self.valid_genomes:
                self.valid_genomes.remove(genome)
                self.filtered_genomes.append(
                    FilteredOutGenome(
                        genome.accession,
                        genome.species_taxid,
                        genome.organism_name,
                        "Failed to process downloaded genome assembly"
                    )
                )
                break
            
    def _cleanup_download_files(self, temp_zip: Path):
        """Clean up temporary download files"""
        try:
            if temp_zip.exists():
                temp_zip.unlink()
            if os.path.exists('ncbi_dataset'):
                shutil.rmtree('ncbi_dataset', ignore_errors=True)
        except Exception as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")
                
    def write_filtered_genomes(self):
        """
        Write information about filtered genomes to CSV files:
        1. Detailed file with all filtered genomes
        2. Summary file with unique species that were dropped
        """
        
        self.logger.debug(f"Writing filtered genomes information to {self.filtered_genomes_output}")
        
        unique_filtered = {}
        for genome in self.filtered_genomes:
            if genome.accession in unique_filtered:
                existing = unique_filtered[genome.accession]
                if genome.reason != existing['filter_reason']:
                    existing['filter_reason'] = f"{existing['filter_reason']}; {genome.reason}"
            else:
                unique_filtered[genome.accession] = {
                    'assembly_accession': genome.accession,
                    'species_taxid': genome.species_taxid,
                    'organism_name': genome.organism_name,
                    'filter_reason': genome.reason
                }
        
        filtered_data = pd.DataFrame(list(unique_filtered.values()))
        
        filtered_data = filtered_data.sort_values(['species_taxid', 'filter_reason'])
        filtered_data.to_csv(self.filtered_genomes_output, index=False)
        self.logger.debug(f"Wrote information about {len(filtered_data)} filtered genomes (after removing duplicates)")
        
        summary_file = self.filtered_genomes_output.replace('.csv', '_species_summary.csv')
        
        species_summary = (
            filtered_data
            .groupby(['species_taxid', 'organism_name'])
            .agg({
                'assembly_accession': 'count',
                'filter_reason': lambda x: '; '.join(sorted(set(';'.join(x).split('; '))))
            })
            .reset_index()
            .rename(columns={
                'assembly_accession': 'number_of_genomes_filtered',
                'filter_reason': 'filter_reasons'
            })
            .sort_values('number_of_genomes_filtered', ascending=False)
        )
        
        species_summary.to_csv(summary_file, index=False)
        self.logger.debug(f"Wrote summary of {len(species_summary)} unique species that were filtered out")
            
    def write_outputs(self):
        """
        Write the filtered data to output files
        """
        self.logger.debug("Writing output files for FungiParser")

        #uuid,assembly_filename,species,species_taxid,assembly_accession
        genome_data = pd.DataFrame([
            {
                'uuid': genome.accession,
                'assembly_filename': genome.assembly_filename,
                'species': genome.organism_name,
                'species_taxid': genome.species_taxid,
                'assembly_accession': genome.accession
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
        taxon_data = taxon_data.drop_duplicates() # Drop duplicate rows
        taxon_data.to_csv(self.taxon_output_filename, index=False)
        self.logger.debug(f"Wrote {len(self.valid_genomes)} accessions to {self.taxon_output_filename}")
        
        
        # Write filtered out genomes
        self.write_filtered_genomes()
        
    def generate_spreadsheets(self):
        """Main method to process data and generate output files"""
        self.process_fungi_data()
        self.download_genomes_bulk()
        self.write_outputs()