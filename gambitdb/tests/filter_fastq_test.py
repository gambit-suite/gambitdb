# tests for the FilterFastq class
import unittest
import os
import filecmp
from gambitdb.FilterFastq  import FilterFastq

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','filter_fastq')

class TestFilterFastq(unittest.TestCase):
    def test_filter_fastq(self):
        self.cleanup()
        fastq_filename = os.path.join(data_dir, 'test.fastq')
        kmer_prefix = 'ATGAC'
        kmer = 11
        min_kmer_freq = 2
        output_kmer_filename = os.path.join(data_dir, 'test_output_kmers.fa')
        verbose = False

        f = FilterFastq(fastq_filename,
                    kmer_prefix,
                    kmer,
                    min_kmer_freq,
                    output_kmer_filename,
                    verbose)
        f.filter_fastq()

        self.assertTrue(filecmp.cmp(output_kmer_filename, os.path.join(data_dir, 'test_expected_output_kmers.fa'), shallow=False))
        self.cleanup()

    def cleanup(self):
        for f in ['test_output_kmers.fa']:
            if os.path.exists(os.path.join(data_dir, f)):
                os.remove(os.path.join(data_dir, f))  