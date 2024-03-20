import pandas
import csv
import numpy as np
import logging
import networkx as nx

class CompressClusters:
    """
    Takes in a pairwise file and removes samples which are very closely related to each other, keeping a representative sample. A maximum distance is passed in and we focus on values below (<=) this threshold.
    """
    def __init__(self, pairwise_distances_filename, output_pairwise_filename, max_distance, representative_genomes, verbose):
        self.pairwise_distances_filename = pairwise_distances_filename
        self.output_pairwise_filename = output_pairwise_filename
        self.max_distance = max_distance

        # filename containing a list of representative genome accessions to keep.
        self.representative_genomes = representative_genomes
        self.verbose = verbose

        self.logger = logging.getLogger('CompressClusters')

    def generate_representative_genomes(self):
        if self.representative_genomes is None:
            return []

        # read in the csv file containing the representative genomes
        with open(self.representative_genomes, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip the header
            representative_genomes = [row[0] for row in reader]
        return representative_genomes

    def compress(self):
        representative_accessions = self.generate_representative_genomes()
        return self.identify_samples_to_remove(representative_accessions)

    #  Gives array([[0.1, 0.2, 0.3],
    #         [0.2, 0.1, 0.4],
    #         [0.3, 0.4, 0.1]])
    def read_pairwise_distances(self):
        self.logger.debug('CompressClusters.read_pairwise_distances')
        pairwise_distances = pandas.read_csv(self.pairwise_distances_filename, index_col=0)
        # sort the pairwise_distances dataframe by the index column
        pairwise_distances = pairwise_distances.sort_index()
        return pairwise_distances

    def identify_samples_to_remove(self, representative_accessions):
        pw = self.read_pairwise_distances()
        pw_filtered = pw[(pw > 0) & (pw <= self.max_distance)]
        pw_filtered = pw_filtered.dropna(axis=0, how='all').dropna(axis=1, how='all')

        # count the number of samples (row names) in pw
        num_samples = len(pw.index)

        samples_to_remove = []
    
        # create a graph of paired values
        G = nx.Graph()
        for index, row in pw_filtered.iterrows():
            non_zero_columns = row.dropna(axis=0, how='all').index.tolist()
            for column in non_zero_columns:
                G.add_edge(index, column)

        # find the connected components in the graph
        connected_components = nx.connected_components(G)

        # get the representative samples for each connected component
        for component in connected_components:
            component_samples = list(component)

            # filter the pw pairwise distance matrix to only include the representative samples as rows and columns. Find the sample which has the lowest distance overall to the others
            min_distance = np.inf
            representative_sample = None
            for sample in component_samples.copy():
                sample_distances = pw_filtered.loc[sample, component_samples]
                sample_distances = sample_distances.dropna()
                avg_distance = np.mean(sample_distances)
                if avg_distance < min_distance:
                    min_distance = avg_distance
                    representative_sample = sample

            # filter the component_samples list to remove the representative sample
            component_samples.remove(representative_sample)

            # if the sample is in the list of representative genomes, remove it from the list so its kept
            if len(representative_accessions) > 0:
                component_samples = [sample for sample in component_samples if sample not in representative_accessions]

            # add the component samples to the list of samples to remove
            samples_to_remove.extend(component_samples)
            if self.verbose:
                print('Compressed cluster: ' + representative_sample + ' -> ' + ','.join(component_samples))
            
            if self.verbose:
                print('Samples to compress: '+ str(len(samples_to_remove)))

        return samples_to_remove, num_samples
