import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
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
        
    def load_data(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
        subspecies_mapping = {}
        
        for name in species_df.index:
            if "subspecies" in name:
                base_name = name.split(" subspecies")[0]
                parent_species.add(base_name)
                subspecies_mapping[base_name] = name
        
        filtered_species = species_df[
            species_df.index.str.contains('subspecies') |
            ~species_df.index.isin(parent_species)
        ]
        
        updated_genomes = genomes_df.copy()
        for parent, subspecies in subspecies_mapping.items():
            mask = updated_genomes['species'] == parent
            if mask.any():
                updated_genomes.loc[mask, 'species'] = subspecies
                
        return filtered_species, updated_genomes, parent_species
    
    def calculate_distances(self, genomes: pd.DataFrame, species: pd.DataFrame, dmat: pd.DataFrame) -> Tuple[List[Tuple], pd.DataFrame, pd.DataFrame]:
        logging.info("Calculating species distances...")
        species, genomes, parent_species = self.process_taxonomy_data(species, genomes)
        
        genome_species = set(genomes['species'].unique())
        species_names = set(species.index)
        common_species = sorted(genome_species.intersection(species_names))
        
        logging.info(f"Species in genomes: {len(genome_species)}")
        logging.info(f"Species in species df: {len(species_names)}")
        logging.info(f"Common species: {len(common_species)}")
        
        genomes = genomes[genomes['species'].isin(common_species)]
        species = species.loc[common_species]
        
        gb = genomes.groupby('species')
        species_inds = [gb.indices[sp] for sp in common_species]
        nspecies = len(species_inds)
        
        diameters = np.zeros(nspecies)
        min_inter = np.zeros((nspecies, nspecies))
        
        overlaps = []
        overlap_details = []
        overlap_counts = {}
        
        for i, inds1 in tqdm(enumerate(species_inds)):
            diameters[i] = dmat.values[np.ix_(inds1, inds1)].max()
            
            for j, inds2 in enumerate(species_inds[:i]):
                mi = dmat.values[np.ix_(inds1, inds2)].min()
                min_inter[i, j] = min_inter[j, i] = mi
        
        for i in range(nspecies):
            for j in range(i):
                d = min_inter[i, j]
                if d <= diameters[i]:
                    overlaps.append((i, j))
                    species1 = common_species[i]
                    species2 = common_species[j]
                    overlap_details.append({
                        'species1': species1,
                        'species2': species2,
                        'species1_diameter': diameters[i],
                        'species2_diameter': diameters[j],
                        'min_distance': d
                    })
                    overlap_counts[species1] = overlap_counts.get(species1, 0) + 1
                
                if d <= diameters[j]:
                    overlaps.append((j, i))
                    species1 = common_species[j]
                    species2 = common_species[i]
                    overlap_details.append({
                        'species1': species1,
                        'species2': species2,
                        'species1_diameter': diameters[j],
                        'species2_diameter': diameters[i],
                        'min_distance': d
                    })
                    overlap_counts[species1] = overlap_counts.get(species1, 0) + 1
        
        overlap_df = pd.DataFrame(overlap_details)
        if not overlap_df.empty:
            overlap_df.to_csv(self.output_dir / 'overlap_details.csv', index=False)
        
        overlap_counts_df = pd.DataFrame.from_dict(overlap_counts, orient='index', columns=['overlap_count'])
        overlap_counts_df.index.name = 'species'
        overlap_counts_df.sort_values('overlap_count', ascending=False).to_csv(self.output_dir / 'overlap_counts.csv')
        
        mi_df = pd.DataFrame(min_inter, index=common_species, columns=common_species)
        species['diameter'] = diameters
        species['ngenomes'] = genomes.groupby('species').size()
        
        return overlaps, species, mi_df
    
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
        
        vis_dir = self.output_dir / 'overlaps'
        vis_dir.mkdir(exist_ok=True)
        
        dmat['assembly_accession'] = dmat.index
        
        for i in overlaps:
            filt_species = species.iloc[list(i)]
            genomes_filt = genomes[genomes['species'].isin(filt_species.index)]
            accessions = genomes_filt['assembly_accession']
            species_data = dmat.loc[accessions, accessions]
            
            species_str = '_'.join(str(e) for e in i)
            
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
            linkage_matrix = linkage(dists, "complete")
            dendrogram(linkage_matrix, color_threshold=self.config.min_distance_threshold)
            plt.tight_layout()
            plt.savefig(vis_dir / f'{species_str}_dendrogram.png', dpi=300)
            plt.close()