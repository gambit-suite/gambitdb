# A class which will take in a species with a high diameter and split clusters into sub clusters.
# These then get their own diameters calculated and if they are still too high, they get split again.
# It takes in a species table, the genome assemblies metadata and the pairwise table.
# This probably needs to be refactored to take in objects instead of reading in files, but 
# it does make it easier to have a direct script for the step.

import pandas as pd
import logging
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, leaves_list,to_tree
from scipy.spatial.distance import squareform
import numpy as np

class SplitSpecies:
    def __init__(self, species_taxon_filename, genome_assembly_metadata_filename, pairwise_distances_filename, maximum_diameter, minimum_cluster_size, verbose):
        self.logger = logging.getLogger(__name__)
        self.species_taxon_filename = species_taxon_filename
        self.genome_assembly_metadata_filename = genome_assembly_metadata_filename
        self.pairwise_distances_filename = pairwise_distances_filename
        self.maximum_diameter = maximum_diameter
        self.minimum_cluster_size = minimum_cluster_size

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def split_high_diameter_species(self):
        species, genome_metadata, pairwise_distances = self.read_input_files()
        high_diameter_species = self.filter_high_diameter_species(species)

        for single_species in high_diameter_species.iterrows():
            genomes_no_outliers = self.split_single_high_diameter_species_into_subspecies(single_species, genome_metadata, pairwise_distances)

    def filter_high_diameter_species(self, species):
        # identify all species with a diameter of > 0.7
        high_diameter_species = species[species['diameter'] > self.maximum_diameter]
        self.logger.debug('split_high_diameter_species: high_diameter_species size: %s' % high_diameter_species.shape[0])
        return high_diameter_species
    
    # take a single species, get corresponding genome assembly metadata and a pairwise distance matrix of that species.
    def split_single_high_diameter_species_into_subspecies(self, single_species, genome_metadata, pairwise_distances):
        # get the genome assembly accessions for this species
        genome_accessions = genome_metadata[genome_metadata['species_taxid'] == single_species[0]]
        # get the pairwise distance matrix for this species and make sure the rows and columns are sorted the same
        pairwise_distances_single = pairwise_distances.loc[genome_accessions.index.tolist(),genome_accessions.index.tolist()]
        pairwise_distances_single = pairwise_distances_single.sort_index()
        pairwise_distances_single = pairwise_distances_single.sort_index(axis=1)

        # get the linkage matrix for this species
        linkage_matrix = self.calculate_linkage_matrix(pairwise_distances_single)
        # get the cluster list based on threshold
        cluster_identity = self.get_cluster_identity(linkage_matrix)
        # get the clusters
        clusters = self.get_clusters(cluster_identity, pairwise_distances_single)

        # remove clusters with <2 genomes per cluster
        clusters = self.remove_clusters_with_too_few_genomes(clusters)
        # add genomes to no outliers file only if curation didn't reduce number of genomes to 1 or 0
        if len(clusters) > 1:
            genomes_no_outliers = pd.concat([genomes_no_outliers,clusters])
        return genomes_no_outliers

    def calculate_linkage_matrix(self, pairwise_distances):
        # reindex pw-distance matrix
        data_mat = pairwise_distances.to_numpy()
        dists = squareform(data_mat)
        linkage_matrix = linkage(dists, "average")
        return linkage_matrix

    def get_cluster_identity(self, linkage_matrix):
        # get cluster list based on threshold
        cluster_identity = fcluster(linkage_matrix,t=0.45,criterion='distance')
        return cluster_identity
    
    def get_clusters(self, cluster_identity, pairwise_distances_single):
        # get the clusters
        clusters = pd.DataFrame(np.array(cluster_identity),columns=['cluster_identity'])
        # add in the accession numbers so that you can link the clusters to genomes
        clusters['assembly_accession']= pairwise_distances_single.index.tolist()
        return clusters

    def read_input_files(self):
        species = pd.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')
        self.logger.debug('read_input_files: species size: %s' % species.shape[0])

        # Read in the genome assembly filenames with path
        genome_metadata = pd.read_csv(self.genome_assembly_metadata_filename)
        genome_metadata = genome_metadata.set_index('assembly_accession')
        self.logger.debug('read_input_files: genome_metadata size: %s' % genome_metadata.shape[0])

        # Read in the pairwise distances file
        pairwise_distances = pd.read_csv(self.pairwise_distances_filename, index_col=0)
        # Make sure the columns and rows are sorted
        pairwise_distances = pairwise_distances.sort_index()
        pairwise_distances = pairwise_distances.sort_index(axis=1)
        # sort the pairwise_distances dataframe by the index column
        self.logger.debug('read_input_files: pairwise_distances size: %s' % pairwise_distances.shape[0])

        return species, genome_metadata, pairwise_distances
