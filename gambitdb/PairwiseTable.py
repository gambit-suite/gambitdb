# A class which will take in a directory of assemblies in FASTA format and generate a table of pairwise distances

import os
import logging
import subprocess
import tempfile

class PairwiseTable:
    def __init__(self, options):
        self.logger = logging.getLogger(__name__)
        self.assembly_directory = options.assembly_directory
        self.signatures_output_filename = options.signatures_output_filename
        self.distance_table_output_filename = options.distance_table_output_filename
        self.kmer = options.kmer
        self.kmer_prefix = options.kmer_prefix

        # FASTA extensions to filter on for assemblies in the assembly directory
        self.valid_extensions = ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz']

        self.verbose = options.verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def generate_sigs_and_pairwise_table(self):
        self.logger.debug('generate_sigs_and_pairwise_table')
        self.generate_genome_signatures()
        self.generate_pairwise_table()
        self.cleanup()

    def generate_genome_signatures(self):
        self.logger.debug('generate_genome_signatures')
        self.assembly_list_filename = self.create_assembly_list()

        self.logger.debug(self.signatures_command(self.assembly_list_filename))
        subprocess.check_call(self.signatures_command(self.assembly_list_filename), shell=True)

    def signatures_command(self, assembly_list_filename):
        # gambit signatures create -k 11 -p ATGAC -o signatures.h5 -l assemblies.fofn
        return 'gambit signatures create -k %s -p %s -o %s -l %s' % (self.kmer, self.kmer_prefix, self.signatures_output_filename, assembly_list_filename)
    
    # Given a directory of assemblies in FASTA format, create a tempory file with the full path of each assembly, one per line and return the tempory filename  
    def create_assembly_list(self):
        self.logger.debug('create_assembly_list')
        assembly_list = tempfile.NamedTemporaryFile(mode='w', delete=False)
        for assembly in os.listdir(self.assembly_directory):
            # limit filenames to a list of prefixes for FASTA files and can include gz files
            if any(assembly.endswith(ext) for ext in self.valid_extensions):
                assembly_list.write('%s/%s\n' % (self.assembly_directory, assembly))
        assembly_list.close()
        return assembly_list.name

    def generate_pairwise_table(self):
        self.logger.debug('generate_pairwise_table')
        self.logger.debug(self.pairwise_table_command())

        subprocess.check_call(self.pairwise_table_command(), shell=True)

    def pairwise_table_command(self):
        return 'gambit dist --qs %s --square -o %s' % (self.signatures_output_filename, self.distance_table_output_filename)

    def cleanup(self):
        self.logger.debug('cleanup')
        os.remove(self.assembly_list_filename)
    