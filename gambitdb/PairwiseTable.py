# A class which will take in a directory of assemblies in FASTA format and generate a table of pairwise distances

import os
import logging
import subprocess
import tempfile

class PairwiseTable:
    """
    A class which will take in a directory of assemblies in FASTA format and generate a table of pairwise distances.
    """
    def __init__(self, assembly_directory, signatures_output_filename,distance_table_output_filename, kmer, kmer_prefix, accessions_to_ignore_file, cpus, verbose):
        """
    Initializes the PairwiseTable class.
    Args:
      assembly_directory (str): The directory containing the assemblies in FASTA format.
      signatures_output_filename (str): The filename of the output signatures file.
      distance_table_output_filename (str): The filename of the output distance table file.
      kmer (int): The kmer size to use.
      kmer_prefix (str): The kmer prefix to use.
      accessions_to_ignore_file (str): The filename of the accessions to ignore file.
      cpus (int): The number of CPUs to use.
      verbose (bool): Whether to print verbose output.
    """
        self.logger = logging.getLogger(__name__)
        self.assembly_directory = assembly_directory
        self.signatures_output_filename = signatures_output_filename
        self.distance_table_output_filename = distance_table_output_filename
        self.kmer = kmer
        self.kmer_prefix = kmer_prefix
        self.accessions_to_ignore_file = accessions_to_ignore_file

        # FASTA extensions to filter on for assemblies in the assembly directory
        self.valid_extensions = ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz']
        self.cpus   = cpus
        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def generate_sigs_and_pairwise_table(self):
        """
    Generates the genome signatures and the pairwise table.
    Returns:
      None
    """
        self.logger.debug('generate_sigs_and_pairwise_table')
        self.generate_genome_signatures()
        self.generate_pairwise_table()
        self.cleanup()

    def generate_genome_signatures(self):
        """
    Generates the genome signatures.
    Returns:
      None
    """
        self.logger.debug('generate_genome_signatures')
        self.assembly_list_filename = self.create_assembly_list()
        # if the self.signatures_output_filename file exists, delete it
        if os.path.isfile(self.signatures_output_filename):
            os.remove(self.signatures_output_filename)

        self.logger.debug(self.signatures_command(self.assembly_list_filename))
        subprocess.check_call(self.signatures_command(self.assembly_list_filename), shell=True)

    def signatures_command(self, assembly_list_filename):
        """
    Generates the command to generate the genome signatures.
    Args:
      assembly_list_filename (str): The filename of the assembly list.
    Returns:
      str: The command to generate the genome signatures.
    """
        # gambit signatures create -k 11 -p ATGAC -o signatures.h5 -l assemblies.fofn
        return 'gambit signatures create -k %s -p %s -o %s -l %s -c %s ' % (self.kmer, self.kmer_prefix, self.signatures_output_filename, assembly_list_filename, self.cpus)
    
    # read in the accessions_file_to_ignore into a list, checking if it is None first
    def read_accessions_to_ignore(self):
        """
    Reads in the accessions_file_to_ignore into a list.
    Returns:
      list: The list of accessions to ignore.
    """
        self.logger.debug('read_accessions_to_ignore')
        accessions_to_ignore = []
        if self.accessions_to_ignore_file is not None and os.path.isfile(self.accessions_to_ignore_file):
            with open(self.accessions_to_ignore_file) as f:
                accessions_to_ignore = f.read().splitlines()
        self.logger.debug('read_accessions_to_ignore: %s accessions to ignore' % len(accessions_to_ignore))
        return accessions_to_ignore

    # Given a directory of assemblies in FASTA format, create a tempory file with the full path of each assembly, one per line and return the tempory filename  
    def create_assembly_list(self):
        """
    Creates a tempory file with the full path of each assembly, one per line and returns the tempory filename.
    Args:
      self.assembly_directory (str): The directory containing the assemblies in FASTA format.
    Returns:
      str: The tempory filename.
    """
        self.logger.debug('create_assembly_list')
        accessions_to_ignore = self.read_accessions_to_ignore()
        assembly_list = tempfile.NamedTemporaryFile(mode='w', delete=False)
        file_counter = 0
        # iterate over all files in the assembly directory and sort by filename
        for assembly in sorted(os.listdir(self.assembly_directory)):
            # limit filenames to a list of prefixes for FASTA files and can include gz files
            if any(assembly.endswith(ext) for ext in self.valid_extensions):

                for accession in accessions_to_ignore:
                    if accession in assembly:
                        self.logger.debug('Skipping assembly %s as it is in the accessions to ignore list' % assembly)
                        continue

                file_counter += 1
                assembly_list.write('%s/%s\n' % (self.assembly_directory, assembly))
        assembly_list.close()
        self.logger.debug('create_assembly_list: %s assemblies from files' % file_counter)
        return assembly_list.name

    def generate_pairwise_table(self):
        """
    Generates the pairwise table.
    Returns:
      None
    """
        self.logger.debug('generate_pairwise_table')
        self.logger.debug(self.pairwise_table_command())
        # if the pairwise distance table file already exists, delete it
        if os.path.isfile(self.distance_table_output_filename):
            os.remove(self.distance_table_output_filename)

        subprocess.check_call(self.pairwise_table_command(), shell=True)

    def pairwise_table_command(self):
        """
    Generates a command to generate a pairwise distance table.
    Args:
      self (PairwiseTable): The instance of the PairwiseTable class.
    Returns:
      str: The command to generate a pairwise distance table.
    Examples:
      >>> PairwiseTable.pairwise_table_command()
      'gambit dist --qs %s --square -o %s -c %s ' % (self.signatures_output_filename, self.distance_table_output_filename, self.cpus)
    """
        return 'gambit dist --qs %s --square -o %s -c %s ' % (self.signatures_output_filename, self.distance_table_output_filename, self.cpus)

    def cleanup(self):
        """
    Removes the assembly list file.
    Args:
      self (PairwiseTable): The instance of the PairwiseTable class.
    Side Effects:
      Removes the assembly list file.
    Examples:
      >>> PairwiseTable.cleanup()
      'cleanup'
    """
        self.logger.debug('cleanup')
        os.remove(self.assembly_list_filename)
    