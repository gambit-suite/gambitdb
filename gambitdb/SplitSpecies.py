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
    def __init__(self, species, genome_assembly_metadata, pairwise_distances_filename, accessions_removed, maximum_diameter, minimum_cluster_size, verbose):
        self.logger = logging.getLogger(__name__)
        self.species = species
        self.genome_assembly_metadata = genome_assembly_metadata
        self.pairwise_distances_filename = pairwise_distances_filename
        self.accessions_removed = accessions_removed
        self.maximum_diameter = maximum_diameter
        self.minimum_cluster_size = minimum_cluster_size

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    # Main method for this class which drives everything else.
    def split_high_diameter_species(self):
        species, genome_metadata, pairwise_distances = self.read_input_files()
        high_diameter_species = self.filter_high_diameter_species(species)

        for single_species in high_diameter_species.iterrows():
            subspecies, genome_metadata, single_species = self.split_single_high_diameter_species_into_subspecies(single_species, genome_metadata, pairwise_distances)
            # the single species gets updated within the method as does the genome_metadata
            # concat subspecies to the species dataframe
            species = pd.concat([species, subspecies], ignore_index=True)
        return species, genome_metadata, self.accessions_removed
    
    def filter_high_diameter_species(self, species):
        # identify all species with a diameter of > 0.7
        high_diameter_species = species[species['diameter'] > self.maximum_diameter]
        self.logger.debug('filter_high_diameter_species: high_diameter_species size: %s' % high_diameter_species.shape[0])
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
        small_clusters, clusters = self.remove_clusters_with_too_few_genomes(clusters)

        self.save_small_clusters_accessions_removed(small_clusters, single_species)
        subspecies, genome_metadata, single_species = self.create_subspecies_from_clusters(clusters, single_species,genome_metadata)

        # update the diameters for each subspecies
        subspecies = self.calculate_diameters_subspecies(subspecies, genome_metadata, pairwise_distances)
        return subspecies, genome_metadata, single_species 

    def calculate_diameters_subspecies(self, subspecies, genome_metadata, pairwise_distances):
        self.logger.debug('calculate_diameters_subspecies')   
        diameters = np.zeros(subspecies.shape[0])

        for cluster in genome_metadata.groupby('species'):
            assembly_accessions = cluster[1].index.tolist()
            species_name = cluster[0]
            inds1 = pairwise_distances.index.get_indexer(assembly_accessions)
            diameter = pairwise_distances.values[np.ix_(inds1, inds1)].max()
            # update the diameter for the subspecies with the matching species_name
            subspecies.loc[subspecies['name'] == species_name,'diameter'] = diameter

        return subspecies

    # Take in the clusters and create subpecies from them, copying the single_species and adding 'subspecies X' where X is an integer and setting report to 0
    def create_subspecies_from_clusters(self, clusters, single_species, genome_metadata):
        # get the attribute names from the single_species dataframe row and use these to create a new dataframe with the columns populated
        new_subspecies_list = []
        
        for index, cluster in enumerate(clusters.groupby('cluster_identity')):
            # create the name of the subspecies
            subspecies_name = str(single_species[1]['name']) + ' subspecies ' + str(index + 1)

            # get all the assembly_accession in the cluster and update the genome_metadata to change the species it subspecies_name
            # this is a very inefficient way to link for the database lookup later
            genome_metadata.loc[cluster[1]['assembly_accession'], 'species'] = subspecies_name

            new_subspecies_list.append([single_species[0], subspecies_name, 'species', single_species[0], single_species[0], single_species[0], 0, len(cluster[1]), 0])

        #  the single species should be reported by with a diameter of zero
        single_species[1]['report'] = 1
        single_species[1]['diameter'] = 0
        
        subspecies = pd.DataFrame(new_subspecies_list, columns=['species_taxid'] + single_species[1].index.tolist())
        subspecies = subspecies.set_index('species_taxid')
        return subspecies, genome_metadata, single_species
    
    def save_small_clusters_accessions_removed(self, small_clusters, single_species):
        # save the accessions of the small clusters to a file
        small_clusters_accessions = small_clusters['assembly_accession'].tolist()
        self.accessions_removed = self.accessions_removed + small_clusters_accessions

        self.logger.debug('Remove small clusters: '  + str(single_species[0]) + '\t' + str(small_clusters_accessions))

    def remove_clusters_with_too_few_genomes(self, clusters):
        # remove clusters with <=X genomes per cluster
        clusters = clusters.groupby('cluster_identity').filter(lambda x: len(x) > self.minimum_cluster_size)
        small_clusters = clusters.groupby('cluster_identity').filter(lambda x: len(x) <= self.minimum_cluster_size)
        return small_clusters, clusters

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
        # Read in the pairwise distances file
        pairwise_distances = pd.read_csv(self.pairwise_distances_filename, index_col=0)
        # Make sure the columns and rows are sorted
        pairwise_distances = pairwise_distances.sort_index()
        pairwise_distances = pairwise_distances.sort_index(axis=1)
        # sort the pairwise_distances dataframe by the index column
        self.logger.debug('read_input_files: pairwise_distances size: %s' % pairwise_distances.shape[0])

        return self.species, self.genome_assembly_metadata, pairwise_distances
