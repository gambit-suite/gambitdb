import pandas
import csv
import numpy as np
import logging
import networkx as nx
import os

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
        samples_to_remove, num_samples = self.identify_samples_to_remove(representative_accessions)
        
        # Write the compressed pairwise distances to output file
        self.write_compressed_pairwise_distances(samples_to_remove)
        
        return samples_to_remove, num_samples
    
    def write_compressed_pairwise_distances(self, samples_to_remove):
        """Write the compressed pairwise distances to output file."""
        self.logger.debug('CompressClusters.write_compressed_pairwise_distances')
        
        # Read the full pairwise distances
        pw = self.read_pairwise_distances()
        
        # Remove samples that were identified for removal
        pw_compressed = pw.drop(index=samples_to_remove, errors='ignore')
        pw_compressed = pw_compressed.drop(columns=samples_to_remove, errors='ignore')
        
        # Determine output format based on input format
        if self.pairwise_distances_filename.endswith('.npy'):
            self._write_compressed_npy(pw_compressed)
        else:
            self._write_compressed_csv(pw_compressed)
    
    def _write_compressed_csv(self, pw_compressed):
        """Write compressed pairwise distances to CSV format."""
        self.logger.debug('CompressClusters._write_compressed_csv')
        pw_compressed.to_csv(self.output_pairwise_filename)
        if self.verbose:
            print(f'Compressed pairwise distances written to {self.output_pairwise_filename}')
    
    def _write_compressed_npy(self, pw_compressed):
        """Write compressed pairwise distances to npy format."""
        self.logger.debug('CompressClusters._write_compressed_npy')
        
        # Ensure output filename has .npy extension
        if not self.output_pairwise_filename.endswith('.npy'):
            output_npy_filename = self.output_pairwise_filename.replace('.csv', '.npy')
        else:
            output_npy_filename = self.output_pairwise_filename
        
        # Write the numpy array
        pw_compressed_array = pw_compressed.values.astype('float32')
        np.save(output_npy_filename, pw_compressed_array)
        
        # Write the index file
        index_filename = output_npy_filename.replace('.npy', '.idx')
        with open(index_filename, 'w') as f:
            for label in pw_compressed.index:
                f.write(f'{label}\n')
        
        if self.verbose:
            print(f'Compressed pairwise distances written to {output_npy_filename}')
            print(f'Index written to {index_filename}')

    #  Gives array([[0.1, 0.2, 0.3],
    #         [0.2, 0.1, 0.4],
    #         [0.3, 0.4, 0.1]])
    def read_pairwise_distances(self):
        self.logger.debug('CompressClusters.read_pairwise_distances')
        
        # Check if file is npy format or csv format
        if self.pairwise_distances_filename.endswith('.npy'):
            return self._read_pairwise_distances_npy()
        else:
            return self._read_pairwise_distances_csv()
    
    def _read_pairwise_distances_csv(self):
        """Read pairwise distances from CSV format."""
        self.logger.debug('CompressClusters._read_pairwise_distances_csv')
        pairwise_distances = pandas.read_csv(self.pairwise_distances_filename, index_col=0)
        # sort the pairwise_distances dataframe by the index column
        pairwise_distances = pairwise_distances.sort_index()
        return pairwise_distances
    
    def _read_pairwise_distances_npy(self):
        """Read pairwise distances from npy format using memory mapping."""
        self.logger.debug('CompressClusters._read_pairwise_distances_npy')
        
        # Read in the pairwise distances index
        index_filename = self.pairwise_distances_filename.replace('.npy', '.idx')
        if not os.path.exists(index_filename):
            raise FileNotFoundError(f"Index file not found: {index_filename}")
        
        self.logger.debug(f"Reading distance matrix index from {index_filename}")
        with open(index_filename, 'r') as f:
            dist_matrix_index_labels = [line.strip() for line in f]
        pairwise_distances_index = pandas.Index(dist_matrix_index_labels)

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
        pairwise_distances = pandas.DataFrame(pairwise_distances_matrix, 
                                            index=pairwise_distances_index, 
                                            columns=pairwise_distances_index)
        
        # sort the pairwise_distances dataframe by the index column
        pairwise_distances = pairwise_distances.sort_index()
        self.logger.debug('_read_pairwise_distances_npy: pairwise_distances size: %s' % pairwise_distances.shape[0])
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
