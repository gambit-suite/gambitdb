# This class will read in a FASTQ file using BioPython, identify all kmers of length 11 beginning with ATGAC, 
# create a hash of these kmers and their frequencies.
# Then it will filter out all kmers with a frequency of 1, and write the remaining kmers to a FASTA file.
# The reason for it is that low abundance kmers are often background noise and interfer with the accuracy of the taxonomic classification.
# This is an inefficient implementation, ideally it should be done using cython or kmc.

import logging
# import biopython to read in a fastq file
from Bio import SeqIO
import gzip

class FilterFastq:
    def __init__(self, fastq_filename, kmer_prefix, kmer, min_kmer_freq, output_kmer_filename, verbose):
        self.logger = logging.getLogger(__name__)
        self.fastq_filename = fastq_filename
        self.kmer_prefix = kmer_prefix
        self.kmer = kmer
        self.min_kmer_freq = min_kmer_freq
        self.output_kmer_filename = output_kmer_filename
        self.verbose = verbose

        self.kmer_hash = {}
        self.kmer_count = 0
        self.kmer_freq = 0

        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def filter_fastq(self):
        self.read_fastq()
        self.filter_kmers()
        self.write_kmers()

    # Read in a FASTQ file using BioPython
    def read_fastq(self):
        logging.debug("Reading FASTQ file: %s", self.fastq_filename)

        # if the filename is gzipped, open it with gzip.open
        if self.fastq_filename.endswith(".gz"):
            with gzip.open(self.fastq_filename, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    self.count_kmers(str(record.seq))
        else:
            for record in SeqIO.parse(self.fastq_filename, "fastq"):
                self.count_kmers(str(record.seq))

    # Count the kmers in the sequence
    def count_kmers(self, sequence):
        for i in range(len(sequence) - self.kmer + 1):
            kmer = sequence[i:i+self.kmer]
            if kmer.startswith(self.kmer_prefix):
                self.kmer_count += 1
                if kmer in self.kmer_hash:
                    self.kmer_hash[kmer] += 1
                else:
                    self.kmer_hash[kmer] = 1

    # Filter out kmers with a frequency of 1
    def filter_kmers(self):
        logging.debug("Filtering kmers with a frequency of %s", self.min_kmer_freq)
        # Remove all values from the hash with a frequency of < self.min_kmer_freq
        for kmer in list(self.kmer_hash):
            if self.kmer_hash[kmer] < self.min_kmer_freq:
                del self.kmer_hash[kmer]

    # Write the kmers to a FASTA file
    def write_kmers(self):
        logging.debug("Writing kmers to FASTA file: %s", self.output_kmer_filename)
        with open(self.output_kmer_filename, 'w') as output_file:
            for i, kmer in enumerate(self.kmer_hash):
                output_file.write(">" + str(i+1) + " " + str(self.kmer_hash[kmer])+"\n")
                output_file.write("AA" + kmer + "AA" + "\n")
        logging.debug("Total kmers: %s", self.kmer_count)
        logging.debug("Total kmers with a frequency of %s: %s", self.min_kmer_freq, len(self.kmer_hash))