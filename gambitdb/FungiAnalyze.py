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
        self.species_diameters = {}
        
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
        """Load input data from CSV files.
        """
        logging.info("Loading input data...")
        genomes = pd.read_csv(self.config.genomes_path)
        species = pd.read_csv(self.config.species_path).set_index('name')
        dmat = pd.read_csv(self.config.distances_path, index_col=0)
        
        accs = genomes['assembly_accession']
        dmat = dmat.loc[accs, accs]
        logging.info(f"Loaded distance matrix with {len(dmat)} rows")
        
        return genomes, species, dmat
    
    def process_taxonomy_data(self, species_df: pd.DataFrame, genomes_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, set]:
        """Process taxonomy data and return filtered DataFrames and set of parent species."""
        logging.info(f"Initial species count: {len(species_df)}")
        
        # Find parents to drop and create subspecies mapping
        parent_species_to_drop = set()
        subspecies_mapping = {}
        for name in species_df.index:
            if "subspecies" in name:
                base_name = name.split(" subspecies")[0]
                parent_species_to_drop.add(base_name)
                subspecies_mapping[base_name] = name
        
        # Filter and update DataFrames
        filtered_species = species_df[
            species_df.index.str.contains('subspecies') |
            ~species_df.index.isin(parent_species_to_drop)
        ]
        
        updated_genomes = genomes_df.copy()
        for parent, subspecies in subspecies_mapping.items():
            mask = updated_genomes['species'] == parent
            if mask.any():
                updated_genomes.loc[mask, 'species'] = subspecies
        
        return filtered_species, updated_genomes, parent_species_to_drop

    def calculate_species_diameter(self, species_genomes: np.ndarray, dmat: np.ndarray) -> float:
        """Calculate the maximum distance between any two genomes within a species.
        """
        return dmat.values[np.ix_(species_genomes, species_genomes)].max()

    def calculate_min_interspecies_distance(self, 
                                        species1_genomes: np.ndarray, 
                                        species2_genomes: np.ndarray, 
                                        dmat: np.ndarray) -> float:
        """Calculate the minimum distance between any genomes of two species.
        """
        return dmat.values[np.ix_(species1_genomes, species2_genomes)].min()

    def find_species_overlaps(self, 
                            species1_idx: int, 
                            species2_idx: int,
                            min_distance: float,
                            diameters: np.ndarray) -> List[Tuple[int, int]]:
        """Check if two species overlap based on their diameters and minimum distance.
        """
        overlaps = []
        if min_distance <= diameters[species1_idx]:
            overlaps.append((species1_idx, species2_idx))
        if min_distance <= diameters[species2_idx]:
            overlaps.append((species2_idx, species1_idx))
        return overlaps

    def calculate_distances(self, 
                       genomes: pd.DataFrame, 
                       species: pd.DataFrame, 
                       dmat: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int]], set]:
        """Calculate species distances using approach from original script."""
        logging.info("Calculating species distances...")
        
        # Process taxonomy data first
        species, genomes, parent_species = self.process_taxonomy_data(species, genomes)

        # Create diameter dictionary from original species data
        self.species_diameters = species['diameter'].to_dict()

        # Find common species between Dataframes
        genome_species = set(genomes['species'].unique())
        species_names = set(species.index)
        common_species = sorted(genome_species.intersection(species_names))

        logging.info(f"Species in genomes: {len(genome_species)}")
        logging.info(f"Species in species df: {len(species_names)}")
        logging.info(f"Common species: {len(common_species)}")

        # Filter Datarames using common species after processing
        genomes = genomes[genomes['species'].isin(common_species)]
        species = species.loc[common_species]
        nspecies = len(common_species)

        # Group genomes by species
        gb = genomes.groupby('species')
        species_inds = [gb.indices[sp_name] for sp_name in common_species]

        # Initialize arrays using the diameter dictionary
        diameters = np.array([self.species_diameters[sp] for sp in common_species])
        
        min_inter = np.zeros((nspecies, nspecies))
        
        # Calculate minimum inter-species distances
        for species1_idx, species1_genomes in tqdm(enumerate(species_inds)):
            for species2_idx, species2_genomes in enumerate(species_inds[:species1_idx]):             
                min_distance = dmat.values[np.ix_(species1_genomes, species2_genomes)].min()
                min_inter[species1_idx, species2_idx] = min_distance
                min_inter[species2_idx, species1_idx] = min_distance
        
        overlaps = []
        for species1_idx in range(nspecies):
            species1_name = common_species[species1_idx]
            for species2_idx in range(species1_idx):
                species2_name = common_species[species2_idx]
                distance = min_inter[species1_idx, species2_idx]
                
                logging.info(f"Checking overlap: {species1_name} (d={diameters[species1_idx]}) vs "
                            f"{species2_name} (d={diameters[species2_idx]}) - distance={distance}")
                
                if distance <= diameters[species1_idx]:
                    overlaps.append((species1_idx, species2_idx))
                if distance <= diameters[species2_idx]:
                    overlaps.append((species2_idx, species1_idx))
        
        return diameters, min_inter, overlaps, parent_species
    
    def create_overlap_report(self, 
                         overlaps: List[Tuple[int, int]], 
                         species: pd.DataFrame,
                         min_inter: np.ndarray,
                         parent_species: set) -> pd.DataFrame:
        """Create report of overlapping taxa, excluding parent species."""
        overlap_records = []
        
        for species_idx, compare_species_idx in overlaps:
            species1 = species.iloc[species_idx]
            species2 = species.iloc[compare_species_idx]
            
            # Skip if either species is a parent
            if (species1.name in parent_species or 
                species2.name in parent_species):
                logging.debug(f"Skipping parent species comparison: {species1.name} - {species2.name}")
                continue
            
            distance = min_inter[species_idx, compare_species_idx]
            
            # Use dictionary lookup for diameters instead of array
            species1_diameter = self.species_diameters[species1.name]
            species2_diameter = self.species_diameters[species2.name]
            
            logging.debug(f"Creating overlap record for {species1.name} (d={species1_diameter}) vs "
                        f"{species2.name} (d={species2_diameter})")
            
            overlap_records.append({
                'species1_taxid': species1['species_taxid'],
                'species1_name': species1.name,
                'species1_diameter': species1_diameter, 
                'species2_taxid': species2['species_taxid'],
                'species2_name': species2.name,
                'species2_diameter': species2_diameter,
                'min_distance': distance,
                'overlap_amount': min(species1_diameter - distance, 
                                    species2_diameter - distance)
            })
        
        return pd.DataFrame(overlap_records)

    def generate_visualizations(self, 
                          species_pair: Tuple[int, int],
                          species: pd.DataFrame,
                          genomes: pd.DataFrame,
                          dmat: pd.DataFrame) -> None:
        """Generate visualizations for a pair of overlapping species.
        """
        # Get the two species we're comparing
        filt_species = species.iloc[list(species_pair)]
        species_names = filt_species.index.tolist()
        
        # Filter genomes for these species
        genomes_filt = genomes[genomes['species'].isin(species_names)]
        accessions = genomes_filt['assembly_accession']
        
        # Get distance data for these genomes
        species_data = dmat.loc[accessions, accessions]
        
        # Create directory for visualizations
        species_str = f"{species_names[0]}_{species_names[1]}"
        vis_dir = self.output_dir / 'visualizations'
        vis_dir.mkdir(exist_ok=True)
        
        # Save distance matrix with metadata
        genomes_filt = genomes_filt[['assembly_accession', 'species_taxid', 'species']]
        merged_data = pd.merge(species_data, 
                            genomes_filt, 
                            left_index=True, 
                            right_on='assembly_accession',
                            how='left')
        merged_data.to_csv(vis_dir / f'{species_str}_pw_dists.csv', index=False)
        
        # Generate heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(species_data, 
                    vmin=0, 
                    vmax=1, 
                    cmap=sns.cm.rocket_r,
                    xticklabels=False,
                    yticklabels=False)
        plt.title(f"Distance Heatmap: {filt_species.iloc[0].name} vs {filt_species.iloc[1].name}")
        plt.tight_layout()
        plt.savefig(vis_dir / f'{species_str}_heatmap.png', dpi=300)
        plt.close()
        
        # Generate dendrogram
        plt.figure(figsize=(12, 8))
        data_mat = species_data.to_numpy()
        dists = squareform(data_mat)
        linkage_matrix = linkage(dists, "average")
        dendrogram(linkage_matrix, 
                color_threshold=self.config.min_distance_threshold,
                leaf_rotation=90)
        plt.title(f"Dendrogram: {filt_species.iloc[0].name} vs {filt_species.iloc[1].name}")
        plt.tight_layout()
        plt.savefig(vis_dir / f'{species_str}_dendrogram.png', dpi=300)
        plt.close()

    def run_analysis(self):
        """Run the species analysis.
        """
        logging.info("Starting species analysis...")
    
        genomes, species, dmat = self.load_data()
        diameters, min_inter, overlaps, parent_species = self.calculate_distances(genomes, species, dmat)
        
        # Get species that exist in the genomes file
        valid_species = genomes['species'].unique()
        species = species[species.index.isin(valid_species)]
        
        # Use same subset for species data
        species_subset = species.iloc[:len(diameters)].copy()
        species_subset['diameter'] = diameters
        species_subset['ngenomes'] = genomes[genomes['species_taxid'].isin(species_subset.index)].groupby('species_taxid').size()
        
        overlap_report = self.create_overlap_report(overlaps, species_subset, min_inter, 
                                                    parent_species)
        overlap_report.to_csv(self.output_dir / 'species_overlaps.csv', index=False)
        
        logging.info(f"Generating visualizations for {len(overlaps)} overlapping species pairs")
        vis_dir = self.output_dir / 'visualizations'
        vis_dir.mkdir(exist_ok=True, parents=True)
        
        # Generate visualizations
        for species_pair in overlaps:
            self.generate_visualizations(species_pair, species, genomes, dmat)
            
        logging.info("Analysis complete!")