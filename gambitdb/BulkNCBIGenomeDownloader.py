import pandas as pd
import os
import zipfile
import logging
import json
import hashlib
import time
from datetime import datetime
import signal
from pathlib import Path
from typing import List, Optional
from dataclasses import dataclass, asdict
import psutil
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from gambitdb.NCBIDatasets import NCBIDatasetClient

# Added for tracking downloads
@dataclass
class GenomeMetadata:
    uuid: str
    species_taxid: str
    assembly_accession: str
    species: str
    assembly_filename: str = None
    download_status: str = "pending"  # pending, downloading, completed, failed
    download_attempts: int = 0
    download_time: Optional[float] = None
    error_message: Optional[str] = None
    file_size: Optional[int] = None

class BulkNCBIGenomeDownloader:
    """
    Downloader optimized for large-scale genome downloads (100,000+ genomes).
    Includes checkpointing, parallel downloads, and robust error handling.
    """
    def __init__(self, 
                 assembly_metadata_csv: str,
                 output_assembly_csv: str,
                 output_fasta_directory: str,
                 api_key: Optional[str] = None,
                 checkpoint_file: str = "download_progress.json",
                 max_workers: int = 4,
                 max_retries: int = 5,
                 rate_limit: float = 0.5,
                 min_disk_space_gb: float = 100.0,
                 batch_size: int = 50,
                 log_file: Optional[str] = None,
                 verbose: bool = False):
        """
        Args:
            assembly_metadata_csv (str): Path to input CSV with genome metadata
            output_assembly_csv (str): Path to output CSV for updated genome metadata
            output_fasta_directory (str): Directory to save FASTA files
            api_key (str): NCBI API key for authentication
            checkpoint_file (str): Path to save download progress
            max_workers (int): Maximum number of parallel download threads
            max_retries (int): Maximum number of retry attempts per batch
            rate_limit (float): Rate limit for API requests
            min_disk_space_gb (float): Minimum required disk space in GB
            batch_size (int): Number of genomes to download in each batch
            log_file (str): Optional log file path
            verbose (bool): Whether to output debug info
        """
        
        
        self.logger = logging.getLogger(__name__)
        
        if verbose:
            log_level = logging.DEBUG
        else:
            log_level = logging.INFO
            
        self.logger.setLevel(log_level)
        
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        # Add file handler if log_file is provided
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(log_level)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
        
        # Configuration parameters
        self.input_csv = assembly_metadata_csv
        self.output_csv = output_assembly_csv
        self.api_key = api_key
        self.checkpoint_file = checkpoint_file
        self.max_workers = max_workers
        self.max_retries = max_retries
        self.rate_limit = rate_limit
        self.min_disk_space_gb = min_disk_space_gb
        self.batch_size = batch_size
        
        # Setup output directory
        self.output_fasta_directory = Path(os.getcwd()) / output_fasta_directory.lstrip('/')
        self.output_fasta_directory.mkdir(exist_ok=True, parents=True)
        
        # Setup NCBI Datasets client
        self.api_client = NCBIDatasetClient(rate_limit=rate_limit)
        
        # Things we want to track over lifetime of the downloader
        self.start_time = time.time()
        self.genomes = []
        self.genomes_by_accession = {}
        self.processed_count = 0
        self.successful_count = 0
        self.failed_count = 0
        self.exit_requested = False
        
        # Setup signal handling for graceful shutdown
        signal.signal(signal.SIGINT, self._handle_exit_signal)
        signal.signal(signal.SIGTERM, self._handle_exit_signal)
        
        # Lock for thread-safe operations
        self.lock = threading.RLock()
        
        # Load data and restore state from checkpoint
        self._load_data()
        self._restore_checkpoint()
        
    def _handle_exit_signal(self, signum, frame):
        """Handle termination signals gracefully"""
        
        self.logger.warning(f"Received signal {signum}. Initiating graceful shutdown...")
        self.exit_requested = True
        
    def _get_available_disk_space_gb(self) -> float:
        """Return available disk space in GB"""
        
        disk_usage = psutil.disk_usage(self.output_fasta_directory)
        return disk_usage.free / (1024 * 1024 * 1024) 
        
    def _check_disk_space(self) -> bool:
        """Check if enough disk space is available"""
        
        available_space_gb = self._get_available_disk_space_gb()
        if available_space_gb < self.min_disk_space_gb:
            self.logger.error(
                f"Insufficient disk space: {available_space_gb:.2f}GB available, "
                f"{self.min_disk_space_gb}GB required"
            )
            return False
        return True
        
    def _load_data(self):
        """Read the input CSV and prepare genome metadata"""
        
        self.logger.info(f"Loading genome data from {self.input_csv}")
        
        try:
            self.df = pd.read_csv(self.input_csv)
            required_columns = ['uuid', 'species_taxid', 'assembly_accession', 'species']
            missing_columns = [col for col in required_columns if col not in self.df.columns]
            
            if missing_columns:
                raise ValueError(f"Missing required columns in input CSV: {missing_columns}")
                
            # Load genome metadata into data class object
            for _, row in self.df.iterrows():
                genome = GenomeMetadata(
                    uuid=row['uuid'],
                    species_taxid=str(row['species_taxid']),
                    assembly_accession=row['assembly_accession'],
                    species=row['species']
                )
                self.genomes.append(genome)
                self.genomes_by_accession[genome.assembly_accession] = genome
                
            self.logger.info(f"Loaded {len(self.genomes)} genome records")
            
        except Exception as exc:
            self.logger.error(f"Error loading input data: {exc}")
            raise
            
    def _save_checkpoint(self):
        """Save current download progress to checkpoint file
            Sometimes NCBI Datasets API can be slow, so we save progress periodically
        """
        with self.lock:
            try:
                # Convert genome objects to dictionaries
                genome_data = [asdict(genome) for genome in self.genomes]
                
                checkpoint = {
                    'timestamp': datetime.now().isoformat(),
                    'processed_count': self.processed_count,
                    'successful_count': self.successful_count, 
                    'failed_count': self.failed_count,
                    'genomes': genome_data
                }
                
                # Write to a temporary file first, then rename to avoid corruption
                temp_checkpoint = f"{self.checkpoint_file}.tmp"
                with open(temp_checkpoint, 'w') as f:
                    json.dump(checkpoint, f, indent=2)
                    
                os.replace(temp_checkpoint, self.checkpoint_file)
                
            except Exception as exc:
                self.logger.error(f"Error saving checkpoint: {exc}")
                
    def _restore_checkpoint(self):
        """Restore download progress from checkpoint file if it exists"""
        if not os.path.exists(self.checkpoint_file):
            self.logger.info("No checkpoint file found, starting fresh download")
            return
            
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
            
            # Restore state from checkpoint
            self.processed_count = checkpoint.get('processed_count', 0)
            self.successful_count = checkpoint.get('successful_count', 0)
            self.failed_count = checkpoint.get('failed_count', 0)
            
            # Update genome metadata status from checkpoint
            for genome_data in checkpoint.get('genomes', []):
                accession = genome_data.get('assembly_accession')
                if accession in self.genomes_by_accession:
                    genome = self.genomes_by_accession[accession]
                    
                    # Restore genome metadata
                    genome.download_status = genome_data.get('download_status', 'pending')
                    genome.download_attempts = genome_data.get('download_attempts', 0)
                    genome.download_time = genome_data.get('download_time')
                    genome.error_message = genome_data.get('error_message')
                    genome.file_size = genome_data.get('file_size')
                    
                    # Check if assembly file exists and update filename
                    if genome.download_status == 'completed':
                        final_fasta = self.output_fasta_directory / f"{accession}.fna"
                        if final_fasta.exists():
                            genome.assembly_filename = str(final_fasta.absolute())
                        else:
                            # File missing despite completed status
                            genome.download_status = 'pending'
                            genome.assembly_filename = None
                            
            # Log recovery statistics
            completed = sum(1 for genome in self.genomes if genome.download_status == 'completed')
            pending = sum(1 for genome in self.genomes if genome.download_status == 'pending')
            failed = sum(1 for genome in self.genomes if genome.download_status == 'failed')
            
            self.logger.info(f"Restored checkpoint: {completed} completed, {pending} pending, {failed} failed")
            
        except Exception as exc:
            self.logger.error(f"Error restoring from checkpoint: {exc}")
            self.logger.info("Starting fresh download")
    
    def process_genomes(self):
        """Process the list of genomes, downloading any that aren't already present
           If download fails, retry up to max_retries
        """
        self.logger.info(f"Processing {len(self.genomes)} genomes")
        
        # Check which genomes still need to be downloaded (if loading back from checkpoint)
        pending_genomes = self._identify_pending_genomes()
        
        if not pending_genomes:
            self.logger.info("No new genomes to download")
            return
        
        # Let's check the available disk space before starting
        available_space = self._get_available_disk_space_gb()
        self.logger.info(f"Available disk space: {available_space:.2f}GB (minimum required: {self.min_disk_space_gb}GB)")
        
        # Best to create batches of genomes to download
        # This will help with retries and also manage the number of concurrent downloads
        batches = self._create_download_batches(pending_genomes)
        self.logger.info(f"Created {len(batches)} batches of ~{self.batch_size} genomes each")
        
        completed_batches = 0
        total_batches = len(batches)
        
        for batch_epoch, batch in enumerate(batches, 1):
            # Skip already processed batches
            if self.exit_requested:
                self.logger.warning("Exit requested. Saving progress and stopping.")
                break
                
            if not self._check_disk_space():
                self.logger.error("Insufficient disk space. Pausing downloads.")
                break
                
            _ = self._hash_batch(batch)
            self.logger.info(
                f"Processing batch {batch_epoch}/{total_batches} ({len(batch)} genomes) "
                f"[{self.successful_count} success, {self.failed_count} failed]"
            )
            
            success = self._process_genome_batch(batch, batch_epoch)
            
            if success:
                completed_batches += 1
                # Let's save progress every ~5 batches
                if batch_epoch % 5 == 0:
                    self._save_checkpoint()
            
            # And let's log the progress every 10 batches
            if batch_epoch % 10 == 0:
                elapsed = time.time() - self.start_time
                genomes_per_sec = self.processed_count / elapsed if elapsed > 0 else 0
                est_remaining = (len(pending_genomes) - self.processed_count) / genomes_per_sec if genomes_per_sec > 0 else 0
                
                self.logger.info(
                    f"Progress: {self.processed_count}/{len(pending_genomes)} genomes processed "
                    f"({self.successful_count} success, {self.failed_count} failed)"
                )
                self.logger.info(
                    f"Rate: {genomes_per_sec:.2f} genomes/sec, "
                    f"Est. remaining: {est_remaining/3600:.1f} hours"
                )
                self.logger.info(
                    f"Disk space: {self._get_available_disk_space_gb():.2f}GB available"
                )
                
        # Finally we can save the final checkpoint
        self._save_checkpoint()
        self.logger.info(
            f"Completed {completed_batches}/{total_batches} batches "
            f"({self.successful_count}/{len(pending_genomes)} genomes downloaded successfully)"
        )
        
    def _identify_pending_genomes(self) -> List[GenomeMetadata]:
        """Identify genomes that still need to be downloaded"""
        
        pending_genomes = []
        
        for genome in self.genomes:
            final_fasta = self.output_fasta_directory / f"{genome.assembly_accession}.fna"
            
            if final_fasta.exists():
                # File exists, mark as completed
                file_size = os.path.getsize(final_fasta)
                genome.assembly_filename = str(final_fasta.absolute())
                if genome.download_status != 'completed':
                    genome.download_status = 'completed'
                    genome.file_size = file_size
                    self.successful_count += 1
                    self.processed_count += 1
                    self.logger.debug(f"Genome {genome.assembly_accession} already exists ({file_size/1024/1024:.2f}MB)")
                else:
                    genome.file_size = file_size
            elif genome.download_status == 'pending' or (
                  genome.download_status == 'failed' and genome.download_attempts < self.max_retries):
                pending_genomes.append(genome)
                
        return pending_genomes
        
    def _hash_batch(self, batch: List[GenomeMetadata]) -> str:
        """Create a hash of the batch for identification"""
        batch_str = "|".join(sorted([genome.assembly_accession for genome in batch]))
        return hashlib.md5(batch_str.encode()).hexdigest()
        
    def _create_download_batches(self, pending_genomes: List[GenomeMetadata]) -> List[List[GenomeMetadata]]:
        """Create batches of genomes for downloading, prioritizing retries"""
        # Sort by priority: new downloads first, then retries with fewer attempts
        sorted_genomes = sorted(
            pending_genomes, 
            key=lambda genome: (genome.download_status != 'pending', genome.download_attempts)
        )
        
        # Let there be batches
        batches = []
        for i in range(0, len(sorted_genomes), self.batch_size):
            batch = sorted_genomes[i:i + self.batch_size]
            batches.append(batch)
            
        return batches
    
    def _process_genome_batch(self, batch: List[GenomeMetadata], batch_num: int) -> bool:
        """
        Process a single batch of genomes
        Returns True if processed successfully, False otherwise.
        """
        
        batch_id = f"batch_{batch_num}_{self._hash_batch(batch)[:8]}"
        self.logger.info(f"Processing batch {batch_id} ({len(batch)} genomes)")
        
        # Update download state
        for genome in batch:
            genome.download_status = 'downloading'
            genome.download_attempts += 1
        
        temp_dir = self.output_fasta_directory / f"temp_{batch_id}"
        temp_dir.mkdir(exist_ok=True)
        temp_zip = temp_dir / f"download.zip"
        
        start_time = time.time()
        success = False
        
        try:
            # Prepare request for NCBI Datasets API
            request_body = {
                "accessions": [genome.assembly_accession for genome in batch],
                "include_annotation_type": ["GENOME_FASTA"],
                "hydrated": "FULLY_HYDRATED",
                "include_tsv": False
            }
            
            headers = {}
            if self.api_key:
                # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/
                headers['api-key'] = self.api_key
            
            # Download the batch
            url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"
            self.logger.info(f"Downloading batch {batch_id}")
            self.api_client.download_file(
                url=url,
                output_path=str(temp_zip),
                method='POST',
                json=request_body
            )
            
            # Process the downloaded zip file
            self.logger.info(f"Extracting and processing batch {batch_id}")
            success = self._process_bulk_download(temp_zip, batch)
            
            # Update metrics
            download_time = time.time() - start_time
            self.logger.info(
                f"Batch {batch_id} completed in {download_time:.1f}s - "
                f"{len(batch)/download_time:.2f} genomes/sec"
            )
            
        except Exception as exc:
            error_msg = f"Download failed for batch {batch_id}: {str(exc)}"
            self.logger.error(error_msg)
            
            # Mark all genomes in batch as failed
            with self.lock:
                for genome in batch:
                    genome.download_status = 'failed'
                    genome.error_message = error_msg
                    genome.download_time = time.time() - start_time
                    self.failed_count += 1
                    self.processed_count += 1
            
            success = False
            
        finally:
            # Clean up temporary files
            self._cleanup_download_files(temp_zip, temp_dir)
            
        return success
    
    def _process_single_genome(self, accession: str, source_path: Path, dest_path: Path) -> tuple:
        """
        Process a single genome file.
        Returns (success, file_size, error_message)
        """
        
        try:
            # Make sure source file exists
            if not os.path.exists(source_path):
                return False, None, f"Source file not found: {source_path}"
                
            # Copy the file to destination
            shutil.copy2(source_path, dest_path)
            
            # Verify file size
            if not dest_path.exists():
                return False, None, "Failed to create destination file"
                
            file_size = dest_path.stat().st_size
            if file_size < 1000:  # Sanity check - files should be much larger
                return False, None, f"Suspiciously small file: {file_size} bytes"
                
            return True, file_size, None
            
        except Exception as exc:
            return False, None, str(exc)
        
    
    def _cleanup_download_files(self, temp_zip: Path, temp_dir: Path = None):
        """Clean up temporary download files and directories completely"""
        
        try:
            # Remove the entire temporary directory and all its contents
            # This includes the zip file, extracted files, and any subdirectories
            if temp_dir and temp_dir.exists():
                self.logger.debug(f"Removing temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir, ignore_errors=True)
                
                # Double-check that temp_dir is completely gone
                if temp_dir.exists():
                    self.logger.warning(f"Temporary directory still exists after cleanup: {temp_dir}")
            
            # Fallback: if temp_zip was passed separately and still exists
            elif temp_zip and temp_zip.exists():
                self.logger.debug(f"Removing temporary zip file: {temp_zip}")
                temp_zip.unlink()
                
                # Clean up extracted files in the same directory as the zip
                if temp_zip.parent.exists():
                    extract_dir = temp_zip.parent
                    # Remove ncbi_dataset directory that gets created during extraction
                    ncbi_dataset_dir = extract_dir / 'ncbi_dataset'
                    if ncbi_dataset_dir.exists():
                        self.logger.debug(f"Removing extracted ncbi_dataset directory: {ncbi_dataset_dir}")
                        shutil.rmtree(ncbi_dataset_dir, ignore_errors=True)
                    
            # Remove any ncbi_dataset directory in current working directory (legacy cleanup)
            if os.path.exists('ncbi_dataset'):
                self.logger.debug("Removing legacy ncbi_dataset directory from working directory")
                shutil.rmtree('ncbi_dataset', ignore_errors=True)
                
        except Exception as exc:
            self.logger.warning(f"Error cleaning up temporary files: {exc}")
    
    def _process_bulk_download(self, zip_path: Path, genomes: List[GenomeMetadata]) -> bool:
        """
        Extract and process files from bulk download.
        Returns True if at least some genomes were processed successfully.
        """
        if not zip_path.exists():
            self.logger.error(f"Zip file not found: {zip_path}")
            return False
            
        # Track successes and failures
        successful_accessions = {}
        failed_accessions = {}
        
        try:
            # Verify zip file integrity
            try:
                with zipfile.ZipFile(zip_path) as zip_ref:
                    # Check file list before extraction
                    fasta_files = [f for f in zip_ref.namelist() if f.endswith('_genomic.fna')]
                    
                    if not fasta_files:
                        self.logger.error(f"No FASTA files found in zip file: {zip_path}")
                        return False
                    
                    self.logger.info(f"Found {len(fasta_files)} FASTA files in the zip")
                    
                    # Extract files in parallel for speed
                    extract_dir = zip_path.parent
                    self.logger.debug(f"Extracting to {extract_dir}")
                    zip_ref.extractall(path=extract_dir)
            except zipfile.BadZipFile:
                self.logger.error(f"Corrupt zip file: {zip_path}")
                return False
                
            # Process each genome with parallel workers for speed
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {}
                for fasta_path in fasta_files:
                    # Parse accession from path
                    try:
                        accession = fasta_path.split('/')[2]
                        source_path = extract_dir / fasta_path
                        dest_path = self.output_fasta_directory / f"{accession}.fna"
                        
                        future = executor.submit(
                            self._process_single_genome, 
                            accession, 
                            source_path, 
                            dest_path
                        )
                        futures[future] = accession
                    except Exception as exc:
                        self.logger.error(f"Error parsing file path {fasta_path}: {exc}")
                
                # Process results as they complete
                for future in as_completed(futures):
                    accession = futures[future]
                    try:
                        success, file_size, error = future.result()
                        if success:
                            successful_accessions[accession] = file_size  # Store as dictionary entry with file size
                        else:
                            failed_accessions[accession] = error  # Store error message
                            self.logger.warning(f"Failed to process {accession}: {error}")
                    except Exception as exc:
                        failed_accessions[accession] = str(exc)  # Store exception message
                        self.logger.error(f"Error processing {accession}: {exc}")
            
            # Update genome metadata
            with self.lock:
                for genome in genomes:
                    accession = genome.assembly_accession
                    if accession in successful_accessions:
                        final_fasta = self.output_fasta_directory / f"{accession}.fna"
                        genome.assembly_filename = str(final_fasta.absolute())
                        genome.download_status = 'completed'
                        genome.file_size = successful_accessions[accession]
                        genome.error_message = None
                        genome.download_time = time.time() - genome.download_time if genome.download_time else None
                        self.successful_count += 1
                        self.processed_count += 1
                    elif accession in failed_accessions:
                        genome.download_status = 'failed'
                        genome.error_message = failed_accessions[accession]
                        self.failed_count += 1
                        self.processed_count += 1
                        
        except Exception as exc:
            self.logger.error(f"Error processing bulk download: {exc}")
            return False

        return True    
    
    def write_output(self):
        """Write the updated data to the output assembly metadtada csv file"""
        
        output_data = pd.DataFrame([
            {
                'uuid': genome.uuid,
                'species_taxid': genome.species_taxid,
                'assembly_accession': genome.assembly_accession,
                'species': genome.species,
                'assembly_filename': genome.assembly_filename
            }
            for genome in self.genomes
        ])
        
        output_data.to_csv(self.output_csv, index=False)
        self.logger.debug(f"Wrote {len(self.genomes)} genomes to {self.output_csv}")
        
        downloaded = sum(1 for genome in self.genomes if genome.assembly_filename is not None)
        self.logger.debug(f"Successfully downloaded {downloaded}/{len(self.genomes)} genomes")
    
    def download(self):
        """Main method to process data and generate output file"""
        self.process_genomes()
        self.write_output()