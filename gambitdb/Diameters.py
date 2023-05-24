# A class which will take in the pairwise distance file, accessions, and species taxon data and calculate diameters and minimums
import logging
import numpy
import pandas

class Diameters:
    def __init__(self, genome_assembly_metadata, pairwise_distances_filename, species_taxon_filename, species_taxon_output_filename, min_inter_output_filename, verbose):
        self.logger = logging.getLogger(__name__)

        # input files
        self.genome_assembly_metadata = genome_assembly_metadata
        self.pairwise_distances_filename = pairwise_distances_filename
        self.species_taxon_filename = species_taxon_filename

        # output files
        self.species_taxon_output_filename = species_taxon_output_filename
        self.min_inter_output_filename = min_inter_output_filename

        self.verbose = verbose
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def calculate_diameters(self):
        self.logger.debug('calculate_diameters')

        # Read in the genome assembly filenames with path
        genome_metadata = pandas.read_csv(self.genome_assembly_metadata)
        genomes_grouped_by_species_taxid = genome_metadata.groupby('species_taxid')

        # Read in the species taxon file
        species = pandas.read_csv(self.species_taxon_filename, index_col=False)
        species = species.set_index('species_taxid')
        number_of_species = species.shape[0]

        # Read in the pairwise distances file
        pairwise_distances = pandas.read_csv(self.pairwise_distances_filename, index_col=0)

        # Create a list of indices for each species in the genome metadata DataFrame
        species_inds = [genomes_grouped_by_species_taxid.indices[species_taxid] for species_taxid in species.index]

        diameters, min_inter = self.calculate_thresholds(number_of_species, species_inds, pairwise_distances)

        # Extend the species taxon table and add in the diameters and number of genomes
        species['diameter'] = diameters
        species['ngenomes'] = genome_metadata.groupby('species_taxid').size()
        # Write out the species taxon table to a new file
        species.to_csv(self.species_taxon_output_filename)

        # Take the min-inter values, add to a dataframe and write out to a new file
        mininter_df = pandas.DataFrame(min_inter, index=species.index, columns=species.index)
        mininter_df.to_csv(self.min_inter_output_filename)

    # Calculate the diameters and minimums
    def calculate_thresholds(self, number_of_species, species_inds, pairwise_distances):
        self.logger.debug('calculate_thresholds')   
        # initalise the diameters and minimums, these get returned at the end
        diameters = numpy.zeros(number_of_species)
        min_inter = numpy.zeros((number_of_species, number_of_species))

        for i, inds1 in enumerate(species_inds):
            # Find the maximum diameters for each species. Basically look at the pairwise distances
            # and find the maximum distance between any two genomes in the species
            diameters[i] = pairwise_distances.values[numpy.ix_(inds1, inds1)].max()
    
            for j, inds2 in enumerate(species_inds[:i]):
                mi = pairwise_distances.values[numpy.ix_(inds1, inds2)].min()
                min_inter[i, j] = min_inter[j, i] = mi

        # return the diameters and minimums
        return diameters, min_inter
