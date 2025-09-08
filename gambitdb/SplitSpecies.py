# A class which will take in a species with a high diameter and split clusters into sub clusters.
# These then get their own diameters calculated and if they are still too high, they get split again.
# It takes in a species table, the genome assemblies metadata and the pairwise table.
# This probably needs to be refactored to take in objects instead of reading in files, but
# it does make it easier to have a direct script for the step.

import pandas as pd
import logging
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, leaves_list, to_tree
from scipy.spatial.distance import squareform
import numpy as np


class SplitSpecies:
    """
    A class which will take in a species with a high diameter and split clusters into sub clusters.
    Args:
      species (DataFrame): The species table.
      genome_assembly_metadata (DataFrame): The genome assemblies metadata.
      pairwise_distances_filename (str): The pairwise table.
      accessions_removed (list): List of accessions removed.
      maximum_diameter (float): Maximum diameter of a species.
      minimum_cluster_size (int): Minimum size of a cluster.
      linkage_method (str): The linkage method.
      verbose (bool): Verbosity of the logger.
    """

    def __init__(
        self,
        species,
        genome_assembly_metadata,
        pairwise_distances_filename,
        accessions_removed,
        maximum_diameter,
        minimum_cluster_size,
        linkage_method,
        verbose,
    ):
        """
        Initializes the SplitSpecies class.
        Args:
          species (DataFrame): The species table.
          genome_assembly_metadata (DataFrame): The genome assemblies metadata.
          pairwise_distances_filename (str): The pairwise table.
          accessions_removed (list): List of accessions removed.
          maximum_diameter (float): Maximum diameter of a species.
          minimum_cluster_size (int): Minimum size of a cluster.
          linkage_method (str): The linkage method.
          verbose (bool): Verbosity of the logger.
        """
        self.logger = logging.getLogger(__name__)
        self.species = species
        self.genome_assembly_metadata = genome_assembly_metadata
        self.pairwise_distances_filename = pairwise_distances_filename
        self.accessions_removed = accessions_removed
        self.maximum_diameter = maximum_diameter
        self.minimum_cluster_size = minimum_cluster_size
        self.linkage_method = linkage_method

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    # Main method for this class which drives everything else.
    def split_high_diameter_species(self):
        """
        Main method for this class which drives everything else.
        Returns:
          species (DataFrame): The species table.
          genome_metadata (DataFrame): The genome assemblies metadata.
          accessions_removed (list): List of accessions removed.
        """
        species, genome_metadata, pairwise_distances = self.read_input_files()
        high_diameter_species = self.filter_high_diameter_species(species)

        for single_species in high_diameter_species.iterrows():
            subspecies, genome_metadata, single_species = (
                self.split_single_high_diameter_species_into_subspecies(
                    single_species, genome_metadata, pairwise_distances
                )
            )
            # the single species gets updated within the method as does the genome_metadata
            # concat subspecies to the species dataframe

            if subspecies is not None and not subspecies.empty:
                species = pd.concat([species, subspecies], ignore_index=False, sort=False)
                species.loc[
                    species["name"] == single_species[1]["name"], "diameter"
                ] = 0.0  # Changed to 0.0 for consistency
            else:
                # The species gets removed so need to remove the accessions
                genome_accessions = genome_metadata[
                    genome_metadata["species_taxid"] == single_species[0]
                ]
                self.accessions_removed.extend(genome_accessions.index.tolist())

        return species, genome_metadata, self.accessions_removed

    def filter_high_diameter_species(self, species):
        """
        Identifies all species with a diameter of > 0.7.
        Args:
          species (DataFrame): The species table.
        Returns:
          DataFrame: High diameter species.
        """
        # identify all species with a diameter of > 0.7
        high_diameter_species = species[species["diameter"] > self.maximum_diameter]
        self.logger.debug(
            "filter_high_diameter_species: high_diameter_species size: %s"
            % high_diameter_species.shape[0]
        )
        return high_diameter_species

    # take a single species, get corresponding genome assembly metadata and a pairwise distance matrix of that species.
    def split_single_high_diameter_species_into_subspecies(
        self, single_species, genome_metadata, pairwise_distances
    ):
        """
        Takes a single species, gets corresponding genome assembly metadata and a pairwise distance matrix of that species.
        Args:
          single_species (tuple): Single species.
          genome_metadata (DataFrame): The genome assemblies metadata.
          pairwise_distances (DataFrame): The pairwise table.
        Returns:
          subspecies (DataFrame): Subspecies.
          genome_metadata (DataFrame): The genome assemblies metadata.
          single_species (tuple): Single species.
        """
        # get the genome assembly accessions for this species
        genome_accessions = genome_metadata[
            genome_metadata["species_taxid"] == single_species[0]
        ]
        # get the pairwise distance matrix for this species and make sure the rows and columns are sorted the same
        pairwise_distances_single = pairwise_distances.loc[
            genome_accessions.index.tolist(), genome_accessions.index.tolist()
        ]
        pairwise_distances_single = pairwise_distances_single.sort_index()
        pairwise_distances_single = pairwise_distances_single.sort_index(axis=1)

        if len(pairwise_distances_single) <= 1:
            self.logger.debug(
                "split_single_high_diameter_species_into_subspecies: pairwise_distances_single has 1 or less rows/columns"
            )
            return None, genome_metadata, single_species

        # get the linkage matrix for this species
        linkage_matrix = self.calculate_linkage_matrix(pairwise_distances_single)
        # get the cluster list based on threshold
        cluster_identity = self.get_cluster_identity(linkage_matrix)
        # get the clusters
        clusters = self.get_clusters(cluster_identity, pairwise_distances_single)

        # remove clusters with <2 genomes per cluster
        small_clusters, clusters = self.remove_clusters_with_too_few_genomes(clusters)

        # Check if we have at least 2 viable clusters, otherwise keep as single species, after talking to Zach this is more conservative
        num_clusters = clusters['cluster_identity'].nunique()
        if num_clusters < 2:
            self.logger.debug(f"Only {num_clusters} viable cluster(s) found, keeping as single species")
            return None, genome_metadata, single_species

        self.save_small_clusters_accessions_removed(small_clusters, single_species)
        subspecies, genome_metadata, single_species = (
            self.create_subspecies_from_clusters(
                clusters, single_species, genome_metadata
            )
        )

        # update the diameters for each subspecies
        subspecies = self.calculate_diameters_subspecies(
            subspecies, genome_metadata, pairwise_distances
        )
        return subspecies, genome_metadata, single_species

    def calculate_diameters_subspecies(
        self, subspecies, genome_metadata, pairwise_distances
    ):
        """
        Calculates the diameters of subspecies.
        Args:
          subspecies (DataFrame): A DataFrame containing the subspecies.
          genome_metadata (DataFrame): A DataFrame containing the genome metadata.
          pairwise_distances (DataFrame): A DataFrame containing the pairwise distances.
        Returns:
          DataFrame: The updated subspecies DataFrame with diameters and number of genomes calculated.
        Examples:
          >>> calculate_diameters_subspecies(subspecies, genome_metadata, pairwise_distances)
        """
        self.logger.debug("calculate_diameters_subspecies")

        for cluster in genome_metadata.groupby("species"):
            assembly_accessions = cluster[1].index.tolist()
            species_name = cluster[0]
            inds1 = pairwise_distances.index.get_indexer(assembly_accessions)
            diameter = pairwise_distances.values[np.ix_(inds1, inds1)].max()
            ngenomes = len(assembly_accessions)
            # update the diameter for the subspecies with the matching species_name
            # Explicitly convert to float to avoid dtype warnings
            subspecies.loc[subspecies["name"] == species_name, "diameter"] = float(diameter)
            subspecies.loc[subspecies["name"] == species_name, "ngenomes"] = ngenomes

        return subspecies

    # Take in the clusters and create subpecies from them, copying the single_species and adding 'subspecies X' where X is an integer and setting report to 0
    def create_subspecies_from_clusters(
        self, clusters, single_species, genome_metadata
    ):
        """
        Creates subspecies from clusters.
        Args:
          clusters (DataFrame): A DataFrame containing the clusters.
          single_species (DataFrame): A DataFrame containing the single species.
          genome_metadata (DataFrame): A DataFrame containing the genome metadata.
        Returns:
          (DataFrame, DataFrame, DataFrame): A tuple containing the subspecies DataFrame, the updated genome_metadata DataFrame, and the updated single_species DataFrame.
        Examples:
          >>> create_subspecies_from_clusters(clusters, single_species, genome_metadata)
        """
        # get the attribute names from the single_species dataframe row and use these to create a new dataframe with the columns populated
        new_subspecies_list = []

        for index, cluster in enumerate(clusters.groupby("cluster_identity")):
            # create the name of the subspecies
            subspecies_name = (
                str(single_species[1]["name"]) + " subspecies " + str(index + 1)
            )

            # get all the assembly_accession in the cluster and update the genome_metadata to change the species it subspecies_name
            # this is a very inefficient way to link for the database lookup later
            genome_metadata.loc[cluster[1]["assembly_accession"], "species"] = (
                subspecies_name
            )

            new_subspecies_list.append(
                [
                    single_species[0],
                    subspecies_name,
                    "species",
                    single_species[0],
                    single_species[0],
                    single_species[0],
                    0.0,
                    len(cluster[1]),
                    0,
                ]
            )

        #  the single species should be reported by with a diameter of zero
        single_species[1]["report"] = 1
        single_species[1]["diameter"] = 0.0

        subspecies = pd.DataFrame(
            new_subspecies_list,
            columns=["species_taxid"] + single_species[1].index.tolist(),
        )
        subspecies = subspecies.set_index("species_taxid")
        
        # Explicitly ensure diameter column is float dtype
        if "diameter" in subspecies.columns:
            subspecies["diameter"] = subspecies["diameter"].astype(float)
        
        return subspecies, genome_metadata, single_species

    def save_small_clusters_accessions_removed(self, small_clusters, single_species):
        """
        Saves the accessions of small clusters to a file.
        Args:
          small_clusters (DataFrame): A DataFrame containing the small clusters.
          single_species (DataFrame): A DataFrame containing the single species.
        Returns:
          None
        Side Effects:
          Updates the accessions_removed attribute.
        Examples:
          >>> save_small_clusters_accessions_removed(small_clusters, single_species)
        """
        # save the accessions of the small clusters to a file
        small_clusters_accessions = small_clusters["assembly_accession"].tolist()
        self.accessions_removed = self.accessions_removed + small_clusters_accessions

        self.logger.debug(
            "Remove small clusters: "
            + str(single_species[0])
            + "\t"
            + str(small_clusters_accessions)
        )

    def remove_clusters_with_too_few_genomes(self, clusters):
        """
        Removes clusters with too few genomes.
        Args:
          clusters (DataFrame): A DataFrame containing the clusters.
        Returns:
          (DataFrame, DataFrame): A tuple containing the small_clusters DataFrame and the updated clusters DataFrame.
        Examples:
          >>> remove_clusters_with_too_few_genomes(clusters)
        """
        # remove clusters with <=X genomes per cluster
        clusters = clusters.groupby("cluster_identity").filter(
            lambda x: len(x) > self.minimum_cluster_size
        )
        small_clusters = clusters.groupby("cluster_identity").filter(
            lambda x: len(x) <= self.minimum_cluster_size
        )
        return small_clusters, clusters

    def calculate_linkage_matrix(self, pairwise_distances):
        """
        Calculates the linkage matrix.
        Args:
          pairwise_distances (DataFrame): A DataFrame containing the pairwise distances.
        Returns:
          ndarray: The linkage matrix.
        Examples:
          >>> calculate_linkage_matrix(pairwise_distances)
        """
        # reindex pw-distance matrix
        data_mat = pairwise_distances.to_numpy()
        dists = squareform(data_mat)
        linkage_matrix = linkage(dists, self.linkage_method)
        return linkage_matrix

    def get_cluster_identity(self, linkage_matrix):
        """
        Gets the cluster identity.
        Args:
          linkage_matrix (ndarray): The linkage matrix.
        Returns:
          ndarray: The cluster identity.
        Examples:
          >>> get_cluster_identity(linkage_matrix)
        """
        # get cluster list based on threshold
        cluster_identity = fcluster(
            linkage_matrix, t=self.maximum_diameter, criterion="distance"
        )
        return cluster_identity

    def get_clusters(self, cluster_identity, pairwise_distances_single):
        """
        Gets the clusters.
        Args:
          cluster_identity (ndarray): The cluster identity.
          pairwise_distances_single (DataFrame): A DataFrame containing the pairwise distances for a single species.
        Returns:
          DataFrame: The clusters DataFrame.
        Examples:
          >>> get_clusters(cluster_identity, pairwise_distances_single)
        """
        # get the clusters
        clusters = pd.DataFrame(
            np.array(cluster_identity), columns=["cluster_identity"]
        )
        # add in the accession numbers so that you can link the clusters to genomes
        clusters["assembly_accession"] = pairwise_distances_single.index.tolist()
        return clusters

    def read_input_files(self):
        """
        Reads in the input files using memory mapping for pairwise distances.
        Args:
          None
        Returns:
          (DataFrame, DataFrame, DataFrame): A tuple containing the species DataFrame, the genome_assembly_metadata DataFrame, and the pairwise_distances DataFrame.
        Side Effects:
          Updates the logger attribute.
        Examples:
          >>> read_input_files()
        """
        # Read in the pairwise distances index
        self.logger.debug(f"Reading distance matrix index from {self.pairwise_distances_filename}")
        with open(self.pairwise_distances_filename.replace('.npy', '.idx'), 'r') as f:
            dist_matrix_index_labels = [line.strip() for line in f]
        
        # Handle case where index labels are stored as byte string representations due to memmap
        # Convert "b'GCF_014932875.1'" to "GCF_014932875.1"
        processed_labels = []
        for label in dist_matrix_index_labels:
            if label.startswith("b'") and label.endswith("'"):
                # Remove b' prefix and ' suffix
                processed_labels.append(label[2:-1])
            else:
                processed_labels.append(label)
        
        pairwise_distances_index = pd.Index(processed_labels)

        # Memory-map the pairwise distances file instead of loading it
        self.logger.debug(f"Memory-mapping distance matrix from {self.pairwise_distances_filename}")
        pairwise_distances_matrix = np.memmap(self.pairwise_distances_filename, dtype='float32', mode='r')
        
        # The shape of the matrix on disk must be inferred from the index length
        n_genomes = len(pairwise_distances_index)
        if pairwise_distances_matrix.size == n_genomes * n_genomes:
            pairwise_distances_matrix = pairwise_distances_matrix.reshape((n_genomes, n_genomes))
        else:
            raise ValueError("The size of the .npy file does not match the index size for a square matrix.")

        # Create DataFrame from memmap with proper index and columns
        pairwise_distances = pd.DataFrame(pairwise_distances_matrix, 
                                            index=pairwise_distances_index, 
                                            columns=pairwise_distances_index)
        
        self.logger.debug('read_input_files: pairwise_distances size: %s' % pairwise_distances.shape[0])
        return self.species, self.genome_assembly_metadata, pairwise_distances
