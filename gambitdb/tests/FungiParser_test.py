import unittest
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
import os
import tempfile
import shutil
from pathlib import Path
from gambitdb.Fungi import FungiParser, GenomeMetadata

class TestFungiParser(unittest.TestCase):
    """
    Test class for FungiParser class.
    """
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        
        #  Lazy sample data
        self.sample_data = """##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
#assembly_accession\tbioproject\tspecies_taxid\torganism_name\tcontig_count
GCF_000002945.2\tPRJNA127\t4896\tSchizosaccharomyces pombe\t3
GCF_003054445.1\tPRJNA545694\t4909\tPichia kudriavzevii\t5
GCF_000243375.1\tPRJNA79345\t4950\tTorulaspora delbrueckii\t8"""
        
        # Write sample data to temporary file
        self.input_file = os.path.join(self.test_dir, "test_input.txt")
        with open(self.input_file, "w") as f:
            f.write(self.sample_data)
            
        # Output filenames
        self.assembly_output = os.path.join(self.test_dir, "assembly_metadata.csv")
        self.taxon_output = os.path.join(self.test_dir, "species_taxon.csv")
        self.filtered_output = os.path.join(self.test_dir, "filtered_genomes.csv")
        self.fasta_dir = os.path.join(self.test_dir, "fasta")
        
        # Mock genome data from NCBI
        self.mock_genome_data = {
            "reports": [{
                "accession": "GCF_000002945.2",
                "assembly_info": {"assembly_name": "ASM294v2"},
                "source_database": "RefSeq",
                "assembly_stats": {"number_of_contigs": 3},
                "organism": {"organism_name": "Schizosaccharomyces pombe"}
            }],
            "total_count": 1
        }
        
        self.mock_taxon_data = {
            "reports": [{
                "taxonomy": {
                    "rank": "species",
                    "classification": {
                        "genus": {"id": 4895}
                    }
                }
            }]
        }

    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.test_dir)

    @patch('gambitdb.Fungi.NCBIDatasetClient')
    def test_initialization(self, mock_client):
        """Test FungiParser initialization"""
        parser = FungiParser(
            fungi_metadata_spreadsheet=self.input_file,
            max_contigs=100,
            minimum_genomes_per_species=2,
            genome_assembly_metadata_output_filename=self.assembly_output,
            output_fasta_directory=self.fasta_dir,
            taxon_output_filename=self.taxon_output,
            filtered_out_genomes_filename=self.filtered_output,
            exclude_atypical=True,
            is_metagenome_derived="metagenome_derived_exclude",
            parent_taxonomic_level="genus",
            verbose=False
        )
        
        self.assertEqual(parser.max_contigs, 100)
        self.assertEqual(parser.minimum_genomes_per_species, 2)
        self.assertTrue(isinstance(parser.df, pd.DataFrame))
        self.assertEqual(len(parser.df), 3)

    @patch('gambitdb.Fungi.NCBIDatasetClient')
    def test_get_parent_taxon(self, mock_client):
        """Test getting parent taxon information"""
        # Setup mock
        mock_instance = mock_client.return_value
        mock_instance.get_json.return_value = self.mock_taxon_data
        
        parser = FungiParser(
            fungi_metadata_spreadsheet=self.input_file,
            max_contigs=100,
            minimum_genomes_per_species=2,
            genome_assembly_metadata_output_filename=self.assembly_output,
            output_fasta_directory=self.fasta_dir,
            taxon_output_filename=self.taxon_output,
            filtered_out_genomes_filename=self.filtered_output,
            exclude_atypical=True,
            is_metagenome_derived="metagenome_derived_exclude",
            parent_taxonomic_level="genus"
        )
        
        parent_taxon, rank = parser.get_parent_taxon("4896")
        self.assertEqual(parent_taxon, 4895)
        self.assertEqual(rank, "species")

    @patch('gambitdb.Fungi.NCBIDatasetClient')
    def test_process_genome_report(self, mock_client):
        """Test processing of genome report data"""
        parser = FungiParser(
            fungi_metadata_spreadsheet=self.input_file,
            max_contigs=100,
            minimum_genomes_per_species=2,
            genome_assembly_metadata_output_filename=self.assembly_output,
            output_fasta_directory=self.fasta_dir,
            taxon_output_filename=self.taxon_output,
            filtered_out_genomes_filename=self.filtered_output,
            exclude_atypical=True,
            is_metagenome_derived="metagenome_derived_exclude",
            parent_taxonomic_level="genus"
        )
        
        report = self.mock_genome_data["reports"][0]
        genome = parser.process_genome_report(report, "4896", 4895, "species")
        
        self.assertEqual(genome.accession, "GCF_000002945.2")
        self.assertEqual(genome.assembly_name, "ASM294v2")
        self.assertEqual(genome.species_taxid, "4896")
        self.assertEqual(genome.parent_taxid, 4895)

    @patch('gambitdb.Fungi.NCBIDatasetClient')
    def test_filter_genome(self, mock_client):
        """Test genome filtering logic"""
        parser = FungiParser(
            fungi_metadata_spreadsheet=self.input_file,
            max_contigs=5,
            minimum_genomes_per_species=2,
            genome_assembly_metadata_output_filename=self.assembly_output,
            output_fasta_directory=self.fasta_dir,
            taxon_output_filename=self.taxon_output,
            filtered_out_genomes_filename=self.filtered_output,
            exclude_atypical=True,
            is_metagenome_derived="metagenome_derived_exclude",
            parent_taxonomic_level="genus"
        )
        
        # Test genome with too many contigs
        genome = GenomeMetadata(
            accession="TEST1",
            assembly_name="test_assembly",
            assembly_source="RefSeq",
            contig_count=10,
            species_taxid="1234",
            organism_name="Test species",
            parent_taxid=123,
            rank="species"
        )
        
        filter_reason = parser._filter_genome(genome)
        self.assertIsNotNone(filter_reason)
        self.assertIn("Contig count", filter_reason)
        
        genome.contig_count = 3
        filter_reason = parser._filter_genome(genome)
        self.assertIsNone(filter_reason)
        
        # Test genome with sp. in name
        genome.organism_name = "Test sp. ABC"
        filter_reason = parser._filter_genome(genome)
        self.assertIsNotNone(filter_reason)
        self.assertIn("sp.", filter_reason)

    def test_handle_assembly_preferences(self):
        """Test assembly source preference handling"""
        parser = FungiParser(
            fungi_metadata_spreadsheet=self.input_file,
            max_contigs=100,
            minimum_genomes_per_species=2,
            genome_assembly_metadata_output_filename=self.assembly_output,
            output_fasta_directory=self.fasta_dir,
            taxon_output_filename=self.taxon_output,
            filtered_out_genomes_filename=self.filtered_output,
            exclude_atypical=True,
            is_metagenome_derived="metagenome_derived_exclude",
            parent_taxonomic_level="genus"
        )
        
        # Create test genomes with different assembly sources
        genomes = [
            GenomeMetadata(
                accession="TEST1",
                assembly_name="assembly1",
                assembly_source="GenBank",
                contig_count=3,
                species_taxid="1234",
                organism_name="Test species",
                parent_taxid=123,
                rank="species"
            ),
            GenomeMetadata(
                accession="TEST2",
                assembly_name="assembly1",
                assembly_source="SOURCE_DATABASE_REFSEQ",
                contig_count=3,
                species_taxid="1234",
                organism_name="Test species",
                parent_taxid=123,
                rank="species"
            )
        ]
        
        preferred_genomes = parser._handle_assembly_preferences(genomes)
        self.assertEqual(len(preferred_genomes), 1)
        self.assertEqual(preferred_genomes[0].accession, "TEST2")

if __name__ == '__main__':
    unittest.main()