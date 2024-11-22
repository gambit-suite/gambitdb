import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

@dataclass
class AnalysisConfig:
    genomes_path: Path
    species_path: Path
    distances_path: Path
    output_dir: Path
    min_distance_threshold: float = 0.7
    debug: bool = False
    verbose: bool = False

class FungiAnalyzer:
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.setup_logging()
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
    def setup_logging(self):
        level = logging.DEBUG if self.config.debug else logging.INFO
        format = '%(asctime)s - %(levelname)s - %(message)s'
        logging.basicConfig(
            level=level,
            format=format,
            handlers=[
                logging.FileHandler('fungi_analysis.log'),
                logging.StreamHandler() if self.config.verbose else logging.NullHandler()
            ]
        )
    
    @property
    def output_dir(self) -> Path:
        return self.config.output_dir
        
    def load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        logging.info("Loading input data...")
        genomes = pd.read_csv(self.config.genomes_path)
        species = pd.read_csv(self.config.species_path).set_index('name')
        dmat = pd.read_csv(self.config.distances_path, index_col=0)
        accs = genomes['assembly_accession']
        dmat = dmat.loc[accs, accs]
        logging.info(f"Loaded distance matrix with {len(dmat)} rows")
        return genomes, species, dmat
    
    def process_taxonomy_data(self, species_df: pd.DataFrame, genomes_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, set]:
        logging.info(f"Processing taxonomy data, initial species count: {len(species_df)}")
        parent_species = set()
        
        # First identify all parent species that have subspecies
        for name in species_df.index:
            if "subspecies" in name:
                base_name = name.split(" subspecies")[0]
                parent_species.add(base_name)
                
        # Filter out parent species from species_df
        filtered_species = species_df[
            species_df.index.str.contains('subspecies') |
            ~species_df.index.isin(parent_species)
        ]
        
        updated_genomes = genomes_df.copy()
        
        # Remove any genomes that are assigned to parent species that have subspecies
        for parent in parent_species:
            # Remove all genomes assigned to parent species
            mask = updated_genomes['species'] == parent
            if mask.any():
                logging.info(f"Removing {mask.sum()} genomes assigned to parent species {parent}")
                updated_genomes = updated_genomes[~mask]
        
        # Write verification files
        filtered_species.to_csv(self.output_dir / 'filtered_species.csv')
        updated_genomes.to_csv(self.output_dir / 'updated_genomes.csv', index=False)
                
        return filtered_species, updated_genomes, parent_species
    
    def calculate_distances(self, genomes: pd.DataFrame, species: pd.DataFrame, dmat: pd.DataFrame) -> Tuple[List[Tuple], pd.DataFrame, pd.DataFrame]:
        logging.info("Calculating species distances...")
        species, genomes, parent_species = self.process_taxonomy_data(species, genomes)
        
        # Get common species between genomes and species df
        genome_species = set(genomes['species'].unique())
        species_names = set(species.index)
        common_species = sorted(genome_species.intersection(species_names))
        
        # Calculate new min_inter_matrix from dmat
        min_inter_matrix = pd.DataFrame(index=common_species, columns=common_species, dtype=float)
        
        # Debug file for distance calculations
        for i, species1 in enumerate(common_species):
            for j, species2 in enumerate(common_species):
                # Diagonal doesn't matter
                if species1 == species2:
                    min_inter_matrix.loc[species1, species2] = 1.0
                    continue
                
                # Get genomes for each species
                species1_genomes = genomes[genomes['species'] == species1]['assembly_accession'].tolist()
                species2_genomes = genomes[genomes['species'] == species2]['assembly_accession'].tolist()
                
                # Get the submatrix for these two species
                species_dmat = dmat.loc[species1_genomes, species2_genomes]
                
                min_distance = float(species_dmat.min().min())
                
                min_inter_matrix.loc[species1, species2] = min_distance
        
        
        overlaps = []
        overlap_details = []
        overlap_counts = {}
        
        for i, species1 in enumerate(common_species):
            for j, species2 in enumerate(common_species[:i]):
                min_dist = float(min_inter_matrix.loc[species1, species2])
                diameter1 = float(species.loc[species1, 'diameter'])
                diameter2 = float(species.loc[species2, 'diameter'])

                # The core condition: if either diameter is greater than min_distance, it's an overlap
                if diameter1 >= min_dist or diameter2 >= min_dist:
                    overlaps.append((i, j))
                    overlaps.append((j, i))
                    overlap_details.append({
                        'species1': species1,
                        'species2': species2,
                        'species1_diameter': diameter1,
                        'species2_diameter': diameter2,
                        'min_distance': min_dist
                    })
                    overlap_counts[species1] = overlap_counts.get(species1, 0) + 1
                    overlap_counts[species2] = overlap_counts.get(species2, 0) + 1
        
        # Save overlap details to CSV, reused later
        if overlap_details:
            pd.DataFrame(overlap_details).to_csv(self.output_dir / 'overlap_details.csv', index=False)
        
        return overlaps, species, min_inter_matrix
    
    def run_analysis(self):
        logging.info("Starting analysis...")
        genomes, species, dmat = self.load_data()
        overlaps, species, min_inter = self.calculate_distances(genomes, species, dmat)
        
        min_inter.to_csv(self.output_dir / 'min-inter.csv')
        np.fill_diagonal(min_inter.values, 1.0)
        min_distances = pd.DataFrame(min_inter.min(axis=1), columns=['min_inter'])
        
        results = species.merge(min_distances, left_index=True, right_index=True, how='left')
        results['overlap'] = results['min_inter'] - results['diameter']
        results.to_csv(self.output_dir / 'species-data.csv')
        
        # Read the actual overlaps from the CSV
        overlap_df = pd.read_csv(self.output_dir / 'overlap_details.csv')
        
        vis_dir = self.output_dir / 'overlaps'
        vis_dir.mkdir(exist_ok=True)
        
        dmat['assembly_accession'] = dmat.index
        
        # Use set to store processed pairs
        processed_pairs = set()
        
        for _, row in overlap_df.iterrows():
            
            pair = frozenset([row['species1'], row['species2']])
            if pair in processed_pairs:
                continue
            processed_pairs.add(pair)
            
            # Get genomes for both species
            genomes_filt = genomes[genomes['species'].isin(pair)]
            accessions = genomes_filt['assembly_accession']
            species_data = dmat.loc[accessions, accessions]
            species_names = list(pair)
            
            species_str = '_vs_'.join(name.replace(' ', '_') for name in species_names)
            
            genomes_filt = genomes_filt[['assembly_accession', 'species_taxid', 'species']]
            merged_data = pd.merge(species_data, genomes_filt, left_index=True, right_on='assembly_accession', how='left')
            merged_data.to_csv(vis_dir / f'{species_str}_pw_dists.csv', index=False)
            
            plt.figure(figsize=(10, 8))
            sns.heatmap(species_data, vmin=0, vmax=1, cmap=sns.cm.rocket_r)
            plt.tight_layout()
            plt.savefig(vis_dir / f'{species_str}_heatmap.png', dpi=300)
            plt.close()
            
            plt.figure(figsize=(12, 8))
            data_mat = species_data.to_numpy()
            dists = squareform(data_mat)
            linkage_matrix = linkage(dists, "average")
            dendrogram(linkage_matrix, color_threshold=self.config.min_distance_threshold)
            plt.tight_layout()
            plt.savefig(vis_dir / f'{species_str}_dendrogram.png', dpi=300)
            plt.close()